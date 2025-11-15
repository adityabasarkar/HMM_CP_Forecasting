library(arrow)
library(dplyr)
library(tidyr)
library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(future.apply)
library(ggplot2)
library(forcats)
library(gridExtra)

setwd("C:/Users/adity_724nfxg/Documents/git clones/HMM_CP_Forecasting")
sourceCpp("src/cpp_files/cpp_files/hmm_conf_pred_final.cpp")

my_data <- read_parquet("raw_data/cm_2024_johan.parquet")

plot(my_data[my_data[,'country_id']==70, 'ged_target'])

print(my_data["month_id"])

##################################################
# Data preprocessing steps
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

cp_subset_data <- function(cid, start_row, num_rows) {
  df_states %>%
    filter(country_id == cid) %>%        # equality test
    ungroup() %>%                        # avoid grouped-slice surprises
    slice(seq.int(from = start_row,
                  to   = min(n(), start_row + num_rows - 1L))) %>% # clamp to n()
    select(state, casualties)
}

cp_subset_data_country <- function(cid) {
  df_states %>%
    filter(country_id == cid) %>%        # equality test
    ungroup() %>%
    select(state, casualties)
}

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

ctr1 <- my_data[my_data$country_id == 1,,drop=F]
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

##################################################
# Random sampling simulation across multiple
# countries + building alpha vs empirical coverage
# graph
##################################################

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
alpha_grid <- seq(0.0, 0.5, by = 0.05)
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
coverage_df$alpha = 1 - coverage_df$alpha


# Plot (unchanged)
# png(filename="plots/reliability_curve.png", width=1000, height=500,
#     units="px", pointsize=12, bg="white", res=NA)
p1 = ggplot(coverage_df, aes(x = alpha, y = coverage, group = test_len)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_line(aes(color = as.factor(test_len))) +
  facet_wrap(~ panel) +
  # scale_x_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, by = 0.05)) +
  # scale_y_continuous(limits = c(0.4, 1)) + 
  guides(color = guide_legend(title = "Test length", ncol = 4)) +
  labs(
    title = "Reliability Curve",
    subtitle = paste0("Calibration length = ", calib_len, "; Test lengths = 1..6"),
    x = expression(alpha),
    y = "Empirical coverage"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )
pdf("plots/reliability_curve.pdf")
grid.arrange(p1, ncol=1, nrow=2)
dev.off()


out_dir <- "outputs/coverage_runs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(out_dir, "results_list.rds"))

##################################################
# Code to display states in selected country
##################################################

plot_country_states <- function(cid) {
  df_states %>%
    filter(country_id == cid) %>%
    count(state) %>%
    mutate(state = factor(state,
                          levels = 1:4,
                          labels = c("Peaceful", "Escalation", "War", "De-escalation"))) %>%
    ggplot(aes(x = "", y = n, fill = state)) +
    geom_bar(stat = "identity", width = 0.6) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Peaceful"     = "#a6cee3",
        "Escalation"   = "#ff7f00",
        "War"          = "#e31a1c",
        "De-escalation"= "#33a02c"
      )
    ) +
    labs(
      title = paste("State distribution for country", cid),
      x = NULL, y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

# Example usage:
plot_country_states(220)

##################################################
# Code to display variable and non-variable countries
##################################################

variable_countries <- df_states %>%
  group_by(country_id) %>%
  summarise(n_states = n_distinct(state), .groups = "drop") %>%
  filter(n_states > 1) %>%
  pull(country_id)

print(variable_countries)

nonvariable_countries <- df_states %>%
  group_by(country_id) %>%
  summarise(n_states = n_distinct(state), .groups = "drop") %>%
  filter(n_states == 1) %>%
  pull(country_id)

print(nonvariable_countries)

##################################################
# Code to build and analyze cp sets.
##################################################

selected_countries <- c(1, 2, 70, 235, 149)
train_len <- 50
test_len <- 6
alpha <- 0.2
accepts <- 0
k <- 4L
m <- 4L


list_of_cp_sets <- setNames(vector("list", length(selected_countries)),
                    as.character(selected_countries))

for (cid in selected_countries) {
  current_subset <- cp_subset_data_country(cid)
  current_pair_list <- df_to_pair_list_markov(current_subset)
  
  total_len <- train_len + test_len
  maximum_start_index <- length(current_pair_list) - total_len + 1
  start <- sample(1:maximum_start_index, 1)
  
  sequence <- if (train_len > 0L) current_pair_list[start:(start + train_len - 1)] else list()
  test_sequence <- current_pair_list[(start + train_len):(start + total_len - 1)]
  
  cp_set <- get_cp_set(sequence, integer(0), k, m, test_len, alpha)
  condensed_set <- Filter(Negate(is.null), cp_set)
  
  # store condensed_set under the country ID
  list_of_cp_sets[[as.character(cid)]] <- condensed_set
}


cp_set_to_counts <- function(cp_set, k = 4L, base0 = TRUE) {
  if (length(cp_set) == 0) return(tibble::tibble())
  # Build an L x N matrix where each column is a candidate sequence
  mat <- do.call(cbind, lapply(cp_set, as.integer))
  if (is.null(dim(mat))) {
    # single candidate; force matrix
    mat <- matrix(mat, ncol = 1L)
  }
  
  # Optionally convert 0..k-1 to 1..k
  if (base0) mat <- mat + 1L
  
  # sanity: all entries should be in 1..k
  if (any(mat < 1L | mat > k, na.rm = TRUE)) {
    stop("cp_set contains states outside 1..k after base adjustment.")
  }
  
  L <- nrow(mat)
  
  counts <- lapply(seq_len(L), function(t) {
    st <- mat[t, ]
    tibble::tibble(
      timepoint = t,
      state = factor(1:k, levels = 1:k),
      count = as.integer(tabulate(st, nbins = k))
    )
  })
  dplyr::bind_rows(counts)
}

# Build a long DF for all countries
cp_density_df <- function(list_of_cp_sets, k = 4L, base0 = TRUE) {
  # keep facet order same as names(list_of_cp_sets)
  cids <- names(list_of_cp_sets)
  if (is.null(cids)) cids <- as.character(seq_along(list_of_cp_sets))
  dfs <- lapply(seq_along(list_of_cp_sets), function(i) {
    cid <- cids[i]
    df  <- cp_set_to_counts(list_of_cp_sets[[i]], k = k, base0 = base0)
    if (nrow(df) == 0) return(tibble())
    df %>% mutate(country_id = factor(cid, levels = cids))
  })
  bind_rows(dfs)
}

plot_cp_density <- function(list_of_cp_sets, k = 4L, base0 = TRUE) {
  df <- cp_density_df(list_of_cp_sets, k = k, base0 = base0)
  if (nrow(df) == 0) {
    message("No CP candidates to plot.")
    return(invisible(NULL))
  }
  
  # Discrete x-axis labels
  df <- df %>%
    dplyr::mutate(
      month_label = factor(
        paste0("Month ", timepoint),
        levels = paste0("Month ", sort(unique(timepoint)))
      ),
      country_label = forcats::fct_inorder(paste0("Country ", as.character(country_id)))
    )
  
  state_labels <- setNames(
    c("Peaceful", "Escalation", "War", "De-escalation")[seq_len(k)],
    as.character(seq_len(k))
  )
  
  ggplot(df, aes(x = month_label, y = count, fill = forcats::fct_inorder(state))) +
    geom_col(position = "fill", width = 0.95) +   # bars touch, no gaps
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +       # snug axis fit
    scale_fill_manual(
      values = c("#a6cee3", "#ff7f00", "#e31a1c", "#33a02c")[seq_len(k)],
      labels = state_labels,
      name = "State"
    ) +
    labs(
      title = "State Composition Across Prediction Horizon",
      subtitle = "Each bar shows relative share of states among CP candidates at each future month",
      x = NULL,
      y = "Percentage of candidates"
    ) +
    facet_wrap(~ country_label, ncol = 1, scales = "free_x") +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = 18, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.key.width = unit(16, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.4, colour = "#dddddd"),
      strip.text = element_text(face = "bold", size = 12, margin = margin(t = 6, b = 6)),
      plot.margin = margin(12, 16, 12, 16)
    )
}

# ---- Example usage ----
# list_of_cp_sets <- results you built earlier, renamed
# Make sure it's named with country IDs in the order you want to display:
# names(list_of_cp_sets) <- as.character(selected_countries)

p <- plot_cp_density(list_of_cp_sets, k = 4L, base0 = TRUE)
print(p)

out_dir <- "outputs/coverage_runs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list_of_cp_sets, file.path(out_dir, "list_of_cp_sets_6_020.rds"))

list_of_cp_sets <- readRDS(file.path(out_dir, "list_of_cp_sets_6_020.rds"))

ctr1 = my_data[my_data$country_id == 1,]
plot(ctr1$ged_target)

##################################################
# function validation experimentation
##################################################

lst <- df_to_pair_list_markov(subset_data)

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