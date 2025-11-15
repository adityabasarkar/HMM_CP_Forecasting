library(arrow)
library(dplyr)
library(tibble)

library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)

setwd("C:/Users/adity_724nfxg/Documents/git clones/HMM_CP_Forecasting")
sourceCpp("src/cpp_files/hmm_conf_pred_final.cpp")


alpha_grid <- seq(0.0, 0.5, by = 0.05)
test_lengths <- 1:6
calib_len <- 50L

N_TRIALS_RANDOM <- 500       # trials per (alpha, test_len) in random panel

# transition matrix
P = matrix(c(0.1, 0.2, 0.3, 0.4, 
             0.4, 0.1, 0.2, 0.3,
             0.3, 0.4, 0.1, 0.2,
             0.2, 0.3, 0.4, 0.1), nrow=4, ncol=4)
k = 4

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

get_rel_curve_data <- function(add_zero) {
  
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
          add_zero,
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

rel_curves_results1 <- get_rel_curve_data(FALSE)
rel_curves_results2 <- get_rel_curve_data(TRUE)

data_dir = "outputs/r_objects"
saveRDS(rel_curves_results1, file.path(data_dir, "mcsim_reliability_curves_no_added_zero.rds"))
saveRDS(rel_curves_results2, file.path(data_dir, "mcsim_reliability_curves_added_zero.rds"))