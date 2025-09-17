library(arrow)
library(dplyr)
library(tidyr)
library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
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

cp_subset_data <- function(cid, start_row, num_rows) {
  df_states %>%
    filter(country_id == cid) %>%        # equality test
    ungroup() %>%                        # avoid grouped-slice surprises
    slice(seq.int(from = start_row,
                  to   = min(n(), start_row + num_rows - 1L))) %>% # clamp to n()
    select(state, casualties)
}

subset_data <- cp_subset_data(70, 250, 53)

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

q_val <- get_cp_p_val(sequence, test_sequence, k, m, alpha)
q_val