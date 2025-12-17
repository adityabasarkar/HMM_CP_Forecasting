library(arrow)
library(dplyr)
library(tibble)

library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)

set.seed(2569)

sourceCpp("src/cpp_files/hmm_conf_pred_final.cpp")

data_dir = "outputs/r_objects"
data = readRDS(file.path(data_dir, "preprocessed_data.rds"))

alpha_grid <- seq(0.0, 0.5, by = 0.05)
test_lengths <- 1:6
calib_len <- 200L
k = 4
N_TRIALS_RANDOM <- 500       # trials per (alpha, test_len) in random panel

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

mc_rel_curve_cp_set <- function() {
  
  rel_curves_results <- tibble(
    test_length        = integer(),
    target_coverage    = numeric(),
    empirical_coverage = numeric(),
  )
  
  for (alpha_val in alpha_grid) {
    for (T1 in test_lengths) {
      
      num_included <- 0L
      num_tests    <- 0L
      
      for (i in 1:N_TRIALS_RANDOM) {
        
        seq = generate_sequence_mm(P, k, calib_len + T1)
        train_seq = seq[1 : calib_len]
        test_seq = seq[(calib_len+1) : (calib_len+T1)]
        
        p_val <- get_p_val(
          train_seq,
          test_seq,
          1000L,
          as.integer(k),
          as.integer(k),
          FALSE,
          "mm"
        )
        
        if (p_val > alpha_val) num_included <- num_included + 1L
        num_tests <- num_tests + 1L
        
      }
      
      print(
        paste(
          "coverage:", num_included/num_tests,
          "Alpha:", alpha_val,
          "Test Length:", T1
        )
      )
      
      empirical_cov <- if (num_tests > 0L) num_included / num_tests else NA_real_
      
      rel_curves_results <- rel_curves_results %>%
        add_row(
          test_length        = T1,
          target_coverage    = 1 - alpha_val,
          empirical_coverage = empirical_cov
        )
    }
  }
  return(rel_curves_results)
}


mc_rel_curve_likelihood <- function() {
  
  rel_curves_results <- tibble(
    test_length        = integer(),
    target_coverage    = numeric(),
    empirical_coverage = numeric(),
  )
  
  for (alpha_val in alpha_grid) {
    for (T1 in test_lengths) {
      
      num_included <- 0L
      num_tests    <- 0L
      
      for (i in seq_len(N_TRIALS_RANDOM)) {
        
        # simulate one sequence: calib part + test part
        mm_seq       <- generate_sequence_mm(P, k, calib_len + T1)
        mm_seq_calib <- mm_seq[1:calib_len]
        mm_seq_test  <- mm_seq[(calib_len + 1):(calib_len + T1)]
        
        # estimate transition matrix (you may be calling into C++ here)
        transition_probs <- estimate_transition_probabilities(mm_seq_calib, 4)
        # transition_probs <- P
        
        # last observed state becomes the fixed prefix for all candidate futures
        start_state <- mm_seq_calib[[length(mm_seq_calib)]][1]
        
        # all candidate futures of length T1 over states 0..3
        grid_df <- expand.grid(
          rep(list(0:(4L - 1L)), T1),
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        )
        
        # product of transitions along a full path (start_state, seq_row)
        path_prob <- function(seq) {
          states01 <- c(start_state, as.integer(seq))
          from <- states01[-length(states01)] + 1L
          to   <- states01[-1] + 1L
          if (any(allowable_transitions[cbind(from, to)] == 0)) {
            return(0)
          }
          prod(transition_probs[cbind(from, to)])
        }
        
        # 1) score all candidate sequences
        probs <- apply(grid_df, 1L, path_prob)
        
        # 2) order by likelihood (desc)
        ord       <- order(probs, decreasing = TRUE)
        probs_ord <- probs[ord]
        
        # 3) take top few so that cumulative >= 1 - alpha
        target <- 1 - alpha_val
        cumul  <- cumsum(probs_ord)
        
        cutoff_idx <- which(cumul >= target)[1]
        if (is.na(cutoff_idx)) cutoff_idx <- length(probs_ord)  # if all zero, take all
        
        # selected candidate set (each is a length-T1 future sequence, without the start state)
        cp_set <- lapply(seq_len(cutoff_idx), function(j) as.integer(grid_df[ord[j], ]))
        
        mm_seq_test_vec <- as.integer(unlist(mm_seq_test))
        
        if (any(vapply(cp_set, identical, logical(1), mm_seq_test_vec))) {
          num_included <- num_included + 1L
        }
        num_tests <- num_tests + 1L
      }
      
      empirical_cov <- if (num_tests > 0L) num_included / num_tests else NA_real_
      
      print(
        paste(
          "coverage:", empirical_cov,
          "Alpha:", alpha_val,
          "Test Length:", T1
        )
      )
      
      rel_curves_results <- rel_curves_results %>%
        add_row(
          test_length        = T1,
          target_coverage    = 1 - alpha_val,
          empirical_coverage = empirical_cov
        )
    }
  }
  
  return(rel_curves_results)
}


# mc_rel_curve_added_zero <- function() {
#   
#   rel_curves_results <- tibble(
#     type = factor("unanchored, anchored"),
#     test_length        = integer(),
#     target_coverage    = numeric(),
#     empirical_coverage = numeric(),
#   )
#   
#   # -- ADD ZERO -- (Unanchored)
#   for (alpha_val in alpha_grid) {
#     for (T1 in test_lengths) {
#       
#       num_included <- 0L
#       num_tests    <- 0L
#       
#       for (i in 1:N_TRIALS_RANDOM) {
#         
#         mm_seq = generate_sequence_mm(P, k, calib_len + T1)
#         mm_seq_calib <- mm_seq[1:calib_len]
#         mm_seq_test <- mm_seq[(calib_len + 1):(calib_len + T1)]
#         
#         p_val <- get_p_val(mm_seq_calib, mm_seq_test, 1000, 4, 4, TRUE, "mm")
#         if (p_val > alpha_val) num_included <- num_included + 1L
#         num_tests <- num_tests + 1L
#         
#       }
#       
#       print(
#         paste(
#           "coverage:", num_included/num_tests,
#           "Alpha:", alpha_val,
#           "Test Length:", T1
#         )
#       )
#       
#       empirical_cov <- if (num_tests > 0L) num_included / num_tests else NA_real_
#       
#       rel_curves_results <- rel_curves_results %>%
#         add_row(
#           type = "Unanchored",
#           test_length        = T1,
#           target_coverage    = 1 - alpha_val,
#           empirical_coverage = empirical_cov
#         )
#     }
#   }
#   
#   # -- ADD ZERO -- (Anchored)
#   for (alpha_val in alpha_grid) {
#     for (T1 in test_lengths) {
#       
#       num_included <- 0L
#       num_tests    <- 0L
#       
#       for (i in 1:N_TRIALS_RANDOM) {
#         
#         mm_seq = generate_sequence_mm(P, k, 500)
#         n_c  <- length(mm_seq)
#         
#         # we need space for calib_len + T1 *plus one more* state after the test seq
#         idx_lim_c <- n_c - (calib_len + T1)   # note: no +1L here now
#         if (idx_lim_c <= 0L) next
#         
#         # all possible starting indices that leave room for the "next" state
#         candidates <- seq_len(idx_lim_c)
#         
#         # next index after test sequence for each candidate
#         next_idx_vec <- candidates + calib_len + T1
#         
#         # keep only those where the next state is 0
#         valid_starts <- candidates[mm_seq[next_idx_vec] == 0L]
#         
#         # if no valid start exists, skip this country
#         if (length(valid_starts) == 0L) next
#         
#         # sample from valid starts
#         start_idx <- sample(valid_starts, 1L)
#         
#         cal_idx   <- start_idx : (start_idx + calib_len - 1L)
#         test_idx  <- (start_idx + calib_len) : (start_idx + calib_len + T1 - 1L)
#         
#         mm_seq_calib <- mm_seq[cal_idx]
#         mm_seq_test  <- mm_seq[test_idx]
#         
#         p_val <- get_p_val(mm_seq_calib, mm_seq_test, 1000, 4, 4, TRUE, "mm")
#         if (p_val > alpha_val) num_included <- num_included + 1L
#         num_tests <- num_tests + 1L
#         
#       }
#       
#       print(
#         paste(
#           "coverage:", num_included/num_tests,
#           "Alpha:", alpha_val,
#           "Test Length:", T1
#         )
#       )
#       
#       empirical_cov <- if (num_tests > 0L) num_included / num_tests else NA_real_
#       
#       rel_curves_results <- rel_curves_results %>%
#         add_row(
#           type = "anchored",
#           test_length        = T1,
#           target_coverage    = 1 - alpha_val,
#           empirical_coverage = empirical_cov
#         )
#     }
#   }
#   
#   return(rel_curves_results)
# }

mc_rel_curve_cp_set_data <- mc_rel_curve_cp_set()
mc_rel_curve_likelihood_data <- mc_rel_curve_likelihood()
# mc_rel_curve_added_zero_data <- mc_rel_curve_added_zero()

saveRDS(mc_rel_curve_cp_set_data, file.path(data_dir, "mc_rel_curve_cp_set_data.rds"))
saveRDS(mc_rel_curve_likelihood_data, file.path(data_dir, "mc_rel_curve_likelihood_data.rds"))
# saveRDS(mc_rel_curve_added_zero_data, file.path(data_dir, "mc_rel_curve_added_zero_data.rds"))