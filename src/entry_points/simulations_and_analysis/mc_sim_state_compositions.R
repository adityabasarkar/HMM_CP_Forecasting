library(arrow)
library(dplyr)
library(tibble)

library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(future)

set.seed(2569)

sourceCpp("src/cpp_files/hmm_conf_pred_final.cpp")
sourceCpp("src/cpp_files/naive_cp.cpp")
sourceCpp("src/cpp_files/weighted_cp.cpp")

data_dir = "outputs/r_objects"
data = readRDS(file.path(data_dir, "preprocessed_data.rds"))

alpha = 0.2
T1 = 6
calib_len = 200L
split_fraction = 0.7
decay = 0.98
k = 4

allowable_transitions <- matrix(
  c(1,1,0,0,
    0,0,1,1,
    0,0,1,1,
    1,1,0,0),
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




mk_pairs <- function(states_vec) lapply(states_vec, function(s) s-1)
generate_sequence_mm <- function(P, k, n) {
  
  return_sequence = list()
  possible_states = 0:(k-1)
  
  first_state = sample(possible_states, 1)
  
  return_sequence <- append(return_sequence, list(first_state))
  
  for (i in 2:n) {
    previous_state = return_sequence[[i - 1]]
    state = sample(possible_states, 1, replace=FALSE, prob=P[previous_state+1,])
    return_sequence <- append(return_sequence, list(state))
  }
  
  return(return_sequence)
}

get_cp_set_mat = function() {
  cp_set_results = list()
  
  mm_seq = generate_sequence_mm(P, k, calib_len + T1)
  mm_seq_calib <- mm_seq[1:calib_len]
  mm_seq_test <- mm_seq[(calib_len + 1):(calib_len + T1)]
  
  # state at time T+0
  start_state <- mm_seq_calib[[length(mm_seq_calib)]]
  
  #---------------------------------------
  # cp set with no added zero
  cp_set <- get_cp_set(mm_seq_calib, c(0), T1, 1000, 4, 4, alpha, FALSE, "mm", allowable_transitions)
  
  # prepend T+0 to each sequence
  cp_set <- lapply(cp_set, function(seq) c(start_state, seq))
  
  ctr_res <- matrix(0, nrow = 4, ncol = T1 + 1)        # set all entries to 0  
  rownames(ctr_res) <- paste("State", 1:4)             # give row-names  
  true_states = unlist(mm_seq)[(calib_len + 1):(calib_len + T1)]
  colnames(ctr_res) <- c(
    paste0("tp=0, TS=", start_state + 1),
    paste0("tp=", 1:T1, ", TS=", true_states + 1)
  )
  
  for (seq in cp_set) {
    for (t in seq_along(seq)) {
      state <- seq[t]
      ctr_res[state+1, t] <- ctr_res[state+1, t] + 1
    }
  }
  
  ctr_prop <- sweep(ctr_res, 2, colSums(ctr_res), "/")
  cp_set_results[["CP Set"]] <- ctr_prop
  
  #---------------------------------------
  # cp set with likelihood
  transition_probs = estimate_transition_probabilities(mm_seq_calib, 4)
  
  # last observed state becomes the fixed prefix for all candidate futures
  start_state <- mm_seq_calib[[length(mm_seq_calib)]]
  
  # all candidate futures of length T1 over states 0..3
  grid_df <- expand.grid(rep(list(0:(4L - 1L)), T1), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  
  # product of transitions along a full path (start_state, seq_row)
  path_prob <- function(seq) {
    states01 <- c(start_state, as.integer(seq))
    from <- states01[-length(states01)] + 1L
    to   <- states01[-1]               + 1L
    if (any(allowable_transitions[cbind(from, to)] == 0)) {
      return(0)
    }
    prod(transition_probs[cbind(from, to)])
  }
  
  # 1) score all candidate sequences
  probs <- apply(grid_df, 1L, path_prob)
  
  # 2) order by likelihood (desc)
  ord <- order(-probs, runif(length(probs)), method = "radix")
  probs_ord <- probs[ord]
  
  # 3) take top few so that cumulative >= 1 - alpha (prefer overshoot)
  target <- 1 - alpha
  cumul <- cumsum(probs_ord)
  k <- which(cumul >= target)[1]
  if (is.na(k)) k <- length(probs_ord)  # if total < target due to zeros, take all
  
  # selected candidate set, now including T+0
  lik_set <- lapply(seq_len(k), function(j) c(start_state, as.integer(grid_df[ord[j], ])))
  
  ctr_res <- matrix(0, nrow = 4, ncol = T1 + 1)        # set all entries to 0  
  rownames(ctr_res) <- paste("State", 1:4)             # give row-names  
  true_states = unlist(mm_seq)[(calib_len + 1):(calib_len + T1)]
  colnames(ctr_res) <- c(
    paste0("tp=0, TS=", start_state + 1),
    paste0("tp=", 1:T1, ", TS=", true_states + 1)
  )
  
  for (seq in lik_set) {
    for (t in seq_along(seq)) {
      state <- seq[t]
      ctr_res[state+1, t] <- ctr_res[state+1, t] + 1
    }
  }
  
  # normalize to proportions (so each column sums to 1)
  ctr_prop <- sweep(ctr_res, 2, colSums(ctr_res), "/")
  
  cp_set_results[["Likelihood Set"]] <- ctr_prop

  # Both baselines use this same simulated history and held-out future path.
  baseline_sets <- list(
    "Naive CP" = get_naive_cp_set(mm_seq_calib, T1, 4L, split_fraction, alpha, allowable_transitions),
    "Weighted CP" = get_weighted_cp_set(mm_seq_calib, T1, 4L, split_fraction, alpha, decay, allowable_transitions)
  )
  for (set_name in names(baseline_sets)) {
    baseline_sets[[set_name]] <- lapply(baseline_sets[[set_name]], function(seq) c(start_state, seq))
    baseline_set <- baseline_sets[[set_name]]
    ctr_res <- matrix(0, nrow = 4, ncol = T1 + 1)
    rownames(ctr_res) <- paste("State", 1:4)
    colnames(ctr_res) <- colnames(ctr_prop)
    for (seq in baseline_set) {
      for (t in seq_along(seq)) {
        state <- seq[t]
        ctr_res[state+1, t] <- ctr_res[state+1, t] + 1
      }
    }
    cp_set_results[[set_name]] <- if (length(baseline_set) > 0L) sweep(ctr_res, 2, colSums(ctr_res), "/") else ctr_res
  }
  
  print("FINISHED")
  
  return(list(
    cp_set_results = cp_set_results,
    cp_set = cp_set,
    lik_set = lik_set,
    naive_set = baseline_sets[["Naive CP"]],
    weighted_set = baseline_sets[["Weighted CP"]]
  ))
}

res <- get_cp_set_mat() # C1, C2, C3
mc_cp_set_results1 <- res$cp_set_results
cp_set <- res$cp_set
lik_set <- res$lik_set
naive_set <- res$naive_set
weighted_set <- res$weighted_set

true_states <- as.integer(sub(".*TS=", "", colnames(mc_cp_set_results1[[1]])))
true_states <- true_states - 1L
true_states <- as.integer(true_states)

is_true_in_cp <- any(sapply(cp_set, function(seq) identical(seq, true_states)))
is_true_in_cp
is_true_in_lik <- any(sapply(lik_set, function(seq) identical(seq, true_states)))
is_true_in_lik
is_true_in_naive <- any(sapply(naive_set, function(seq) identical(seq, true_states)))
is_true_in_weighted <- any(sapply(weighted_set, function(seq) identical(seq, true_states)))

sum_table_mc_sim <- tibble::tibble(
  set_type = c("CP Set", "Likelihood Set", "Naive CP", "Weighted CP"),
  num      = c(length(cp_set), length(lik_set), length(naive_set), length(weighted_set)),
  true_in  = c(is_true_in_cp, is_true_in_lik, is_true_in_naive, is_true_in_weighted)
)

data_dir = "outputs/r_objects"
saveRDS(mc_cp_set_results1, file.path(data_dir, "sim_compare_compositions.rds"))
saveRDS(sum_table_mc_sim,
        file.path(data_dir, "sim_compare_set_summary_table.rds"))
