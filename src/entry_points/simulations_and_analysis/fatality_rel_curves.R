library(arrow)
library(dplyr)
library(tibble)

library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)

sourceCpp("src/cpp_files/hmm_conf_pred_final.cpp")


# Retrieve preprocessed data from rds
data_dir = "outputs/r_objects"
data = readRDS(file.path(data_dir, "preprocessed_data.rds"))

alpha_grid <- seq(0.0, 0.5, by = 0.05)
test_lengths <- 1:6
calib_len <- 200L
set.seed(2569)


mk_pairs <- function(states_vec) lapply(states_vec, function(s) s-1)
cp_original_rel_curve <- function() {
  
  rel_curves_results <- tibble(
    test_length        = integer(),
    target_coverage    = numeric(),
    empirical_coverage = numeric(),
  )
  
  unique_cids <- unique(data$country_id)
  
  for (alpha_val in alpha_grid) {
    for (T1 in test_lengths) {
      
      num_included <- 0L
      num_tests    <- 0L
      
      for (cid in unique_cids) {
        
        df_c <- data %>% filter(country_id == cid) %>% arrange(month_id)
        n_c  <- nrow(df_c)
        
        idx_lim_c <- n_c - (calib_len + T1) + 1L
        if (idx_lim_c <= 0L) next
        
        start_idx <- sample(seq_len(idx_lim_c), 1)
        cal_idx   <- start_idx : (start_idx + calib_len - 1L)
        test_idx  <- (start_idx + calib_len) : (start_idx + calib_len + T1 - 1L)
        
        s_cal  <- df_c[[ "state" ]][cal_idx]
        s_test <- df_c[[ "state" ]][test_idx]
        
        mm_seq_calib <- mk_pairs(s_cal)
        mm_seq_test  <- mk_pairs(s_test)
        
        p_val <- get_p_val(mm_seq_calib, mm_seq_test, 1000, 4, 4, FALSE, "mm")
        if (p_val > alpha_val) num_included <- num_included + 1L
        num_tests <- num_tests + 1L
        
        print(
          paste(
            "FINISHED CID:", cid, 
            "P-value:", p_val,
            "Alpha:", alpha_val,
            "Test Length:", T1
          )
        )
      }
      
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


likelihood_rel_curve <- function() {
  
  rel_curves_results <- tibble(
    test_length        = integer(),
    target_coverage    = numeric(),
    empirical_coverage = numeric(),
  )
  
  # --- Random sampling panel ---
  unique_cids <- unique(data$country_id)
  
  for (alpha_val in alpha_grid) {
    for (T1 in test_lengths) {
      
      num_included <- 0L
      num_tests    <- 0L
      
      for (cid in unique_cids) {
        
        df_c <- data %>% filter(country_id == cid) %>% arrange(month_id)
        n_c  <- nrow(df_c)
        
        idx_lim_c <- n_c - (calib_len + T1) + 1L
        if (idx_lim_c <= 0L) next
        
        start_idx <- sample(seq_len(idx_lim_c), 1)
        cal_idx   <- start_idx : (start_idx + calib_len - 1L)
        test_idx  <- (start_idx + calib_len) : (start_idx + calib_len + T1 - 1L)
        
        s_cal  <- df_c[[ "state" ]][cal_idx]
        s_test <- df_c[[ "state" ]][test_idx]
        
        mm_seq_calib <- mk_pairs(s_cal)
        mm_seq_test  <- mk_pairs(s_test)
        
        transition_probs = estimate_transition_probabilities(mm_seq_calib, 4)
        
        # last observed state becomes the fixed prefix for all candidate futures
        start_state <- mm_seq_calib[[length(mm_seq_calib)]][1]
        
        # all candidate futures of length T1 over states 0..3
        grid_df <- expand.grid(rep(list(0:(4L - 1L)), T1), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        
        allowable_transitions <- matrix(
          c(1,1,0,0,
            0,0,1,1,
            0,0,1,1,
            1,1,0,0),
          nrow = 4,
          ncol = 4,
          byrow = TRUE
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
        ord <- order(-probs, runif(length(probs)), method = "radix")
        probs_ord <- probs[ord]
        
        # 3) take top few so that cumulative >= 1 - alpha (prefer overshoot)
        target <- 1 - alpha_val
        cumul <- cumsum(probs_ord)
        k <- which(cumul >= target)[1]
        if (is.na(k)) k <- length(probs_ord)  # if total < target due to zeros, take all
        
        # selected candidate set (each is a length-T1 future sequence, without the start state)
        cp_set <- lapply(seq_len(k), function(j) as.integer(grid_df[ord[j], ]))
        
        mm_seq_test_vec <- as.integer(unlist(mm_seq_test))
        
        if (any(sapply(cp_set, identical, mm_seq_test_vec))) num_included <- num_included + 1L
        num_tests <- num_tests + 1L
        
        print(
          paste(
            "FINISHED CID:", cid,
            "Alpha:", alpha_val,
            "Test Length:", T1
          )
        )
      }
      
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


# cp_added_zero_rel_curve <- function() {
#   
#   rel_curves_results <- tibble(
#     type = factor("unanchored, anchored"),
#     test_length        = integer(),
#     target_coverage    = numeric(),
#     empirical_coverage = numeric(),
#   )
#   
#   unique_cids <- unique(data$country_id)
#   
#   # -- ADD ZERO -- (Unanchored)
#   for (alpha_val in alpha_grid) {
#     for (T1 in test_lengths) {
#       
#       num_included <- 0L
#       num_tests    <- 0L
#       
#       for (cid in unique_cids) {
#         
#         df_c <- data %>% filter(country_id == cid) %>% arrange(month_id)
#         n_c  <- nrow(df_c)
#         
#         idx_lim_c <- n_c - (calib_len + T1) + 1L
#         if (idx_lim_c <= 0L) next
#         
#         start_idx <- sample(seq_len(idx_lim_c), 1)
#         cal_idx   <- start_idx : (start_idx + calib_len - 1L)
#         test_idx  <- (start_idx + calib_len) : (start_idx + calib_len + T1 - 1L)
#         
#         s_cal  <- df_c[[ "state" ]][cal_idx]
#         s_test <- df_c[[ "state" ]][test_idx]
#         
#         mm_seq_calib <- mk_pairs(s_cal)
#         mm_seq_test  <- mk_pairs(s_test)
#         
#         p_val <- get_p_val(mm_seq_calib, mm_seq_test, 1000, 4, 4, TRUE, "mm")
#         if (p_val > alpha_val) num_included <- num_included + 1L
#         num_tests <- num_tests + 1L
#         
#         print(
#           paste(
#             "FINISHED CID:", cid, 
#             "P-value:", p_val,
#             "Alpha:", alpha_val,
#             "Test Length:", T1
#           )
#         )
#       }
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
#       for (cid in unique_cids) {
#         
#         df_c <- data %>% filter(country_id == cid) %>% arrange(month_id)
#         n_c  <- nrow(df_c)
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
#         valid_starts <- candidates[df_c[["state"]][next_idx_vec] == 1L]
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
#         s_cal  <- df_c[[ "state" ]][cal_idx]
#         s_test <- df_c[[ "state" ]][test_idx]
#         
#         mm_seq_calib <- mk_pairs(s_cal)
#         mm_seq_test  <- mk_pairs(s_test)
#         
#         p_val <- get_p_val(mm_seq_calib, mm_seq_test, 1000, 4, 4, TRUE, "mm")
#         if (p_val > alpha_val) num_included <- num_included + 1L
#         num_tests <- num_tests + 1L
#         
#         print(
#           paste(
#             "FINISHED CID:", cid, 
#             "P-value:", p_val,
#             "Alpha:", alpha_val,
#             "Test Length:", T1
#           )
#         )
#       }
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


cp_original_rel_data_fatality <- cp_original_rel_curve() # A4
likelihood_rel_data_fatality <- likelihood_rel_curve() # A5
# cp_added_zero_rel_data_fatality <- cp_added_zero_rel_curve() # A6

saveRDS(cp_original_rel_data_fatality, file.path(data_dir, "cp_original_rel_data_fatality.rds"))
saveRDS(likelihood_rel_data_fatality, file.path(data_dir, "likelihood_rel_data_fatality.rds"))
# saveRDS(cp_added_zero_rel_data_fatality, file.path(data_dir, "cp_added_zero_rel_data_fatality.rds"))