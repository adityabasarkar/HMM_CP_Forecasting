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

target_coverage_grid <- seq(1.0, 0.5, by = -0.05)
test_lengths <- list(1, 3, 5)
calib_len <- 200L
k <- 4L
N_TRIALS_RANDOM <- 100L

allowable_transitions <- matrix(
  c(1,1,0,0,
    0,0,1,1,
    0,0,1,1,
    1,1,0,0),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)

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

mc_avg_set_sizes <- function() {
  
  results <- tibble(
    test_length = integer(),
    target_coverage = numeric(),
    alpha = numeric(),
    avg_set_size_cp = numeric(),
    avg_set_size_lik = numeric()
  )
  
  for (T1 in test_lengths) {
    
    grid_df <- expand.grid(
      rep(list(0:(k - 1L)), T1),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    
    for (cov_target in target_coverage_grid) {
      
      alpha_val <- 1 - cov_target
      
      cp_sizes <- numeric(N_TRIALS_RANDOM)
      lik_sizes <- numeric(N_TRIALS_RANDOM)
      
      for (i in seq_len(N_TRIALS_RANDOM)) {
        
        mm_seq <- generate_sequence_mm(P, k, calib_len + T1)
        mm_seq_calib <- mm_seq[1:calib_len]
        
        cp_set <- get_cp_set(
          mm_seq_calib,
          as.integer(c()),
          as.integer(T1),
          as.integer(1000L),
          as.integer(k),
          as.integer(k),
          alpha_val,
          FALSE,
          "mm",
          allowable_transitions,
          TRUE
        )
        cp_sizes[i] <- length(cp_set)
        
        transition_probs <- estimate_transition_probabilities(mm_seq_calib, 4)
        start_state <- mm_seq_calib[[length(mm_seq_calib)]][1]
        
        path_prob <- function(seq) {
          states01 <- c(start_state, as.integer(seq))
          from <- states01[-length(states01)] + 1L
          to <- states01[-1] + 1L
          if (any(allowable_transitions[cbind(from, to)] == 0)) {
            return(0)
          }
          prod(transition_probs[cbind(from, to)])
        }
        
        probs <- apply(grid_df, 1L, path_prob)
        ord <- order(-probs, runif(length(probs)), method = "radix")
        probs_ord <- probs[ord]
        
        cumul <- cumsum(probs_ord)
        cutoff_idx <- which(cumul >= cov_target)[1]
        if (is.na(cutoff_idx)) cutoff_idx <- length(probs_ord)
        lik_sizes[i] <- cutoff_idx
      }
      
      results <- results %>%
        add_row(
          test_length = T1,
          target_coverage = cov_target,
          alpha = alpha_val,
          avg_set_size_cp = mean(cp_sizes),
          avg_set_size_lik = mean(lik_sizes)
        )
      
      print(paste(
        "T1:", T1,
        "target_coverage:", cov_target,
        "alpha:", alpha_val,
        "avg_cp_size:", results$avg_set_size_cp[nrow(results)],
        "avg_lik_size:", results$avg_set_size_lik[nrow(results)]
      ))
    }
  }
  
  results
}

avg_set_sizes <- mc_avg_set_sizes()

out_path <- file.path(data_dir, "mc_avg_set_sizes.rds")
saveRDS(avg_set_sizes, out_path)
print(paste("Saved:", out_path))
