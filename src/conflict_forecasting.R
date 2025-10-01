library(arrow)
library(dplyr)
library(tidyr)
library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(future.apply)
library(ggplot2)
sourceCpp("hmm_conformal_prediction.cpp")


my_data <- read_parquet("cm_2024_johan.parquet")

plot(my_data[my_data[,'country_id']==70, 'ged_target'])

print(my_data["month_id"])

##################################################
# Code To Add "State" Variable
##################################################
df_full <- my_data %>%
  complete(month_id = unique(month_id), country_id,
           fill = list(ged_target = 0L))

# 2) Compute previous-month fatalities by country
df_states <- df_full %>%
  group_by(country_id) %>%
  mutate(prev_fatal   = dplyr::lag(ged_target, default = 0L),
         curr_fatal   = ged_target,
         prev_pos     = prev_fatal > 0L,
         curr_pos     = curr_fatal > 0L,
         state = case_when(
           !prev_pos & !curr_pos ~ 1L, # Peaceful
           !prev_pos &  curr_pos ~ 2L, # Escalation
           prev_pos &  curr_pos ~ 3L, # War
           prev_pos & !curr_pos ~ 4L  # De-escalation
         ),
         casualties = case_when(
           ged_target > 0 ~ 1,
           ged_target <= 0 ~ 0
         )) %>%
  ungroup() %>%
  arrange(country_id, month_id)

##################################################
# Code To Subset data
##################################################

cp_subset_data <- function(cid, start_row, num_rows) {
  df_states %>%
    filter(country_id == cid) %>%        # equality test
    ungroup() %>%                        # avoid grouped-slice surprises
    slice(seq.int(from = start_row,
                  to   = min(n(), start_row + num_rows - 1L))) %>% # clamp to n()
    select(state, casualties)
}

subset_data <- cp_subset_data(83, 0, 400)

cp_subset_data_country <- function(cid) {
  df_states %>%
    filter(country_id == cid) %>%        # equality test
    ungroup() %>%
    select(state, casualties)
}

##################################################
# Code To convert data within dataframes 
# to list of integer pairs
##################################################

df_to_pair_list_hidden_markov <- function(df, state_col = "state", obs_col = "casualties", state_zero_index = FALSE, obs_zero_index = TRUE) {
  st <- as.integer(df[[state_col]])
  ob <- as.integer(df[[obs_col]])
  if (!state_zero_index) {
    st <- st - 1L
  }
  if (!obs_zero_index) {
    ob <- ob - 1L
  }
  lapply(seq_len(nrow(df)), function(i) c(st[i], ob[i]))
}

df_to_pair_list_markov <- function(df, state_col = "state", state_zero_index = FALSE) {
  st <- as.integer(df[[state_col]])
  if (anyNA(st)) stop("state column contains NA after integer coercion.")   # <-- added
  if (!state_zero_index) {
    st <- st - 1L
  }
  lapply(seq_len(nrow(df)), function(i) c(st[i], st[i]))  # integer length-2 vectors
}

lst <- df_to_pair_list_markov(subset_data)

##################################################
# Random sampling simulation across multiple
# countries
##################################################

unique_cids <- unique(df_states$country_id)
train_len <- 50
test_len <- 3
alpha <- 0.2
accepts <- 0
k <- 4L
m <- 4L
for (cid in unique_cids) {
  current_subset <- cp_subset_data_country(cid)
  current_pair_list <- df_to_pair_list_markov(current_subset)
  total_len <- train_len + test_len
  maximum_start_index <- length(current_pair_list) - total_len + 1
  start <- sample(1:maximum_start_index, 1)
  
  sequence      <- if (train_len > 0L) current_pair_list[start:(start + train_len - 1)] else list()
  test_sequence <- current_pair_list[(start + train_len):(start + total_len - 1)]
  
  p_val <- get_cp_p_val(sequence, test_sequence, k, m, alpha)
  if (p_val > alpha) {
    accepts <- accepts + 1
  }
  print(paste("FINISHED CID:", cid, "P-value:", p_val))
}
proportion_accepts <- accepts / length(unique_cids)
proportion_accepts

##################################################
# Random sampling simulation with sliding
# window
##################################################

test_cid <- 70
train_len <- 50
test_len <- 3
alpha <- 0.2
accepts <- 0
tests <- 0
current_subset <- cp_subset_data_country(test_cid)
current_pair_list <- df_to_pair_list_markov(current_subset)
k <- 4L
m <- 4L
for (i in 1:(length(current_pair_list) - (test_len + train_len) + 1)) {
  start <- i
  total_len <- train_len + test_len
  
  sequence      <- if (train_len > 0L) current_pair_list[start:(start + train_len)] else list()
  test_sequence <- current_pair_list[(start + train_len + 1L):(start + total_len - 1)]
  
  # ---- Infer k, m consistent with zero-index conversion above ----
  get_states <- function(pairs_list) vapply(pairs_list, function(v) v[[1]], integer(1))
  get_obs    <- function(pairs_list) vapply(pairs_list, function(v) v[[2]], integer(1))
  
  all_states <- get_states(lst)
  all_obs    <- get_obs(lst)
  
  # Because we subtract 1 when state_zero_index = FALSE, use max + 1L here
  k <- max(all_states) + 1L                                          # <-- changed
  m <- max(all_obs)    + 1L                                          # <-- changed
  
  p_val <- get_cp_p_val(sequence, test_sequence, k, m, alpha)
  if (p_val > alpha) {
    accepts <- accepts + 1
  }
  tests <- tests + 1
  print(paste("FINISHED i:", i, "P-value:", p_val))
}
proportion_accepts <- accepts / tests
proportion_accepts


unique(df_states$country_id)
unique(lst$Value)

##################################################
# Random sampling simulation across multiple
# countries
##################################################

get_states <- function(pairs_list) vapply(pairs_list, function(v) v[[1]], integer(1))
get_obs    <- function(pairs_list) vapply(pairs_list, function(v) v[[2]], integer(1))

infer_k_m <- function(pairs_list) {
  st <- get_states(pairs_list)
  ob <- get_obs(pairs_list)
  list(k = max(st) + 1L, m = max(ob) + 1L)
}

# Build pairs list for a given country ID (Markov version as in your code)
pair_list_for_cid <- function(cid) {
  df <- cp_subset_data_country(cid)
  df_to_pair_list_markov(df)
}

# Compute coverage for a single (alpha, test_len) using a rolling window on one country
coverage_rolling <- function(cid, calib_len, test_len, alpha, max_windows = Inf) {
  pairs <- pair_list_for_cid(cid)
  N <- length(pairs)
  total_len <- calib_len + test_len
  if (N < total_len) return(NA_real_)
  k <- 4L
  m <- 4L
  # All valid window starts
  starts <- seq_len(N - total_len + 1L)
  # Optionally subsample windows for speed
  if (length(starts) > max_windows) {
    set.seed(42)
    starts <- sample(starts, max_windows)
  }
  accepts <- 0L
  tests   <- 0L
  for (s in starts) {
    sequence      <- if (calib_len > 0L) pairs[s:(s + calib_len - 1L)] else list()
    test_sequence <- pairs[(s + calib_len):(s + total_len - 1L)]
    p_val <- get_cp_p_val(sequence, test_sequence, k, m, alpha)
    accepts <- accepts + as.integer(p_val > alpha)
    tests   <- tests + 1L
  }
  if (tests == 0L) NA_real_ else accepts / tests
}

# Compute coverage for a single (alpha, test_len) using random country/time samples
coverage_random <- function(calib_len, test_len, alpha) {
  cids <- unique(df_states$country_id)
  
  # Keep only CIDs with enough length
  valid_cids <- cids[vapply(cids, function(cid) {
    length(pair_list_for_cid(cid)) >= (calib_len + test_len)
  }, logical(1))]
  
  if (length(valid_cids) == 0L) return(NA_real_)
  
  set.seed(123)
  accepts <- 0L
  tests   <- 0L
  
  for (cid in valid_cids) {
    pairs <- pair_list_for_cid(cid)
    N <- length(pairs)
    total_len <- calib_len + test_len
    max_start <- N - total_len + 1L
    if (max_start < 1L) next
    
    # one random start per CID
    s <- sample.int(max_start, 1L)
    
    # infer k, m from this CID's pairs (consistent with zero-index scheme)
    k <- 4L
    m <- 4L
    
    sequence      <- if (calib_len > 0L) pairs[s:(s + calib_len - 1L)] else list()
    test_sequence <- pairs[(s + calib_len):(s + total_len - 1L)]
    
    p_val <- get_cp_p_val(sequence, test_sequence, k, m, alpha)
    accepts <- accepts + as.integer(p_val > alpha)
    tests   <- tests + 1L
  }
  
  if (tests == 0L) NA_real_ else accepts / tests
}

# -------- Grid over alphas and test lengths --------
alpha_grid <- seq(0.50, 1.00, by = 0.05)
test_lengths <- 1:6
calib_len <- 50L

# You can adjust these to trade off speed vs. fidelity
ROLLING_COUNTRY <- 70L       # country_id for rolling window panel
MAX_WINDOWS_PER_LEN <- 1000  # cap windows per (alpha, test_len) for speed
N_TRIALS_RANDOM <- 400       # trials per (alpha, test_len) in random panel

results <- list()

# Rolling panel
for (L in test_lengths) {
  for (a in alpha_grid) {
    cov_val <- coverage_rolling(
      cid = ROLLING_COUNTRY,
      calib_len = calib_len,
      test_len = L,
      alpha = a,
      max_windows = MAX_WINDOWS_PER_LEN
    )
    results[[length(results) + 1L]] <- data.frame(
      panel = "Rolling window",
      test_len = L,
      alpha = a,
      coverage = cov_val
    )
    print(paste("FINISHED L:", L, "FINISHED a:", a, "FINISHED cov_val:", cov_val))
  }
}

# Random selection panel
for (L in test_lengths) {
  for (a in alpha_grid) {
    cov_val <- coverage_random(
      calib_len = calib_len,
      test_len = L,
      alpha = a
    )
    results[[length(results) + 1L]] <- data.frame(
      panel = "Random selection",
      test_len = L,
      alpha = a,
      coverage = cov_val
    )
    print(paste("FINISHED L:", L, "FINISHED a:", a, "FINISHED cov_val:", cov_val))
  }
}

coverage_df <- dplyr::bind_rows(results)

# Plot (unchanged)
ggplot(coverage_df, aes(x = alpha, y = coverage, group = test_len)) +
  geom_abline(intercept = 1, slope = -1, linetype = "dashed") +
  geom_line(aes(color = as.factor(test_len))) +
  facet_wrap(~ panel) +
  scale_x_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, by = 0.05)) +
  guides(color = guide_legend(title = "Test length", ncol = 4)) +
  labs(
    title = "Empirical Coverage vs. Alpha",
    subtitle = paste0("Calibration length = ", calib_len, "; Test lengths = 1..6"),
    x = expression(alpha),
    y = "Empirical coverage"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

out_dir <- "outputs/coverage_runs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(results, file.path(out_dir, "results_list.rds"))

# ---- Split (ensure non-empty test window) ----
calib_len <- 50L
N <- length(lst)                                                   # <-- added
if (calib_len >= N) stop("calib_len must be < number of rows.")    # <-- added
sequence      <- if (calib_len > 0L) lst[seq_len(calib_len)] else list()
test_sequence <- lst[(calib_len + 1L):N]

# ---- Infer k, m consistent with zero-index conversion above ----
get_states <- function(pairs_list) vapply(pairs_list, function(v) v[[1]], integer(1))
get_obs    <- function(pairs_list) vapply(pairs_list, function(v) v[[2]], integer(1))

all_states <- get_states(lst)
all_obs    <- get_obs(lst)

# Because we subtract 1 when state_zero_index = FALSE, use max + 1L here
k <- max(all_states) + 1L                                          # <-- changed
m <- max(all_obs)    + 1L                                          # <-- changed

alpha <- 0.1
k <- 4
m <- 2
estimate_transition_probabilities(sequence, k)
estimate_observation_probabilities(sequence, k, m)
p_val <- get_cp_p_val(sequence, test_sequence, k, m, alpha, TRUE)
p_val