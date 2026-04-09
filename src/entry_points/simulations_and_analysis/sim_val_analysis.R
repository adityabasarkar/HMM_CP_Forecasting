library(arrow)
library(dplyr)
library(tibble)

library(RcppArmadillo)
library(Rcpp, quietly = TRUE)
library(RcppArmadillo, quietly = TRUE)
library(RcppDist, quietly = TRUE)

set.seed(2569)

sourceCpp("src/cpp_files/hmm_conf_pred_final.cpp")

data_dir <- "outputs/r_objects"
data <- readRDS(file.path(data_dir, "preprocessed_data.rds"))

alpha <- 0.2
calib_len <- 200L
k <- 4L
num_trials <- 100L
test_lengths <- 1:6

allowable_transitions <- matrix(
  c(1, 1, 0, 0,
    0, 0, 1, 1,
    0, 0, 1, 1,
    1, 1, 0, 0),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)

# averaged transition matrix
unique_ctrs <- unique(data$country_id)
average_mat <- matrix(0, nrow = k, ncol = k)
for (ctr in unique_ctrs) {
  est <- estimate_transition_probabilities(as.list(data$state[data$country_id == ctr] - 1), k)
  if (est[4, 4] > 0) {
    print(ctr)
  }
  average_mat <- average_mat + est
}
P <- average_mat / length(unique_ctrs)
P_masked <- P * allowable_transitions
row_sums <- rowSums(P_masked)
P <- P_masked
P[row_sums > 0, ] <- P_masked[row_sums > 0, ] / row_sums[row_sums > 0]

generate_sequence_mm <- function(P, k, n) {
  return_sequence <- list()
  possible_states <- 0:(k - 1)

  first_state <- sample(possible_states, 1)
  return_sequence <- append(return_sequence, list(first_state))

  for (i in 2:n) {
    previous_state <- return_sequence[[i - 1]]
    state <- sample(possible_states, 1, replace = FALSE, prob = P[previous_state + 1, ])
    return_sequence <- append(return_sequence, list(state))
  }

  return(return_sequence)
}

summarize_counts <- function(x) {
  tibble(
    Q1 = as.numeric(stats::quantile(x, 0.25, names = FALSE, type = 7)),
    Mean = mean(x),
    Median = stats::median(x),
    Q3 = as.numeric(stats::quantile(x, 0.75, names = FALSE, type = 7))
  )
}

get_trial_counts <- function(T1) {
  mm_seq <- generate_sequence_mm(P, k, calib_len + T1)
  mm_seq_calib <- mm_seq[1:calib_len]

  cp_res <- get_cp_set(
    mm_seq_calib,
    c(0),
    T1,
    1000,
    4,
    4,
    alpha,
    FALSE,
    "mm",
    allowable_transitions
  )

  eq_list <- as.numeric(cp_res$eq_list)

  transition_probs <- estimate_transition_probabilities(mm_seq_calib, 4)
  start_state <- mm_seq_calib[[length(mm_seq_calib)]]

  grid_df <- expand.grid(
    rep(list(0:(4L - 1L)), T1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  path_prob <- function(seq) {
    states01 <- c(start_state, as.integer(seq))
    from <- states01[-length(states01)] + 1L
    to <- states01[-1] + 1L
    if (any(allowable_transitions[cbind(from, to)] == 0)) {
      return(NA_real_)
    }
    prod(transition_probs[cbind(from, to)])
  }

  probs <- apply(grid_df, 1L, path_prob)
  probs <- probs[!is.na(probs)]

  target <- 1 - alpha
  probs_ord <- sort(probs, decreasing = TRUE)
  cumul <- cumsum(probs_ord)
  cutoff_idx <- which(cumul >= target)[1]
  if (is.na(cutoff_idx)) {
    cutoff_idx <- length(probs_ord)
  }

  boundary_prob <- probs_ord[cutoff_idx]
  lik_count <- sum(probs_ord[cutoff_idx:length(probs_ord)] == boundary_prob)

  list(
    cp_counts = eq_list,
    lik_count = lik_count
  )
}

summary_rows <- list()

for (T1 in test_lengths) {
  cp_counts_all <- numeric()
  lik_counts_all <- numeric(num_trials)

  for (trial_idx in seq_len(num_trials)) {
    trial_counts <- get_trial_counts(T1)
    cp_counts_all <- c(cp_counts_all, trial_counts$cp_counts)
    lik_counts_all[trial_idx] <- trial_counts$lik_count
    cat("Test Length:", T1, "| Trial:", trial_idx, "\n")
  }

  summary_rows[[length(summary_rows) + 1L]] <- bind_cols(
    tibble(Method = paste0("CP (T1 = ", T1, ")")),
    summarize_counts(cp_counts_all)
  )

  summary_rows[[length(summary_rows) + 1L]] <- bind_cols(
    tibble(Method = paste0("Likelihood (T1 = ", T1, ")")),
    summarize_counts(lik_counts_all)
  )
}

sim_val_summary_table <- bind_rows(summary_rows)

saveRDS(
  sim_val_summary_table,
  file.path(data_dir, "sim_val_summary_table.rds")
)
