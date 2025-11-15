library(arrow)
library(dplyr)
library(tibble)

library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)

setwd("C:/Users/adity_724nfxg/Documents/git clones/HMM_CP_Forecasting")
sourceCpp("src/cpp_files/hmm_conf_pred_final.cpp")


# Retrieve preprocessed data from rds
data_dir = "outputs/r_objects"
data = readRDS(file.path(data_dir, "preprocessed_data.rds"))

alpha_grid <- seq(0.0, 0.5, by = 0.05)
test_lengths <- 1:6
calib_len <- 50L

# You can adjust these to trade off speed vs. fidelity
ROLLING_COUNTRY <- 70L       # country_id for rolling window panel
MAX_WINDOWS_PER_LEN <- 1000  # cap windows per (alpha, test_len) for speed
N_TRIALS_RANDOM <- 400       # trials per (alpha, test_len) in random panel

# rolling window tests
rolling_country_df = data[data$country_id == ROLLING_COUNTRY,]
rolling_country_n = nrow(rolling_country_df)
mk_pairs <- function(states_vec) lapply(states_vec, function(s) s-1)

get_rel_curve_data <- function(add_zero) {
  
  rel_curves_results <- tibble(
    sampling_type      = factor(levels = c("rolling", "random")),
    test_length        = integer(),
    target_coverage    = numeric(),
    empirical_coverage = numeric(),
  )
  
  for (alpha_val in alpha_grid) {
    for (T1 in test_lengths) {
      
      num_included <- 0L
      num_tests    <- 0L
      
      idx_lim <- rolling_country_n - (calib_len + T1) + 1L
      if (idx_lim > 0L) {
        # optional cap
        starts <- seq_len(idx_lim)
        if (length(starts) > MAX_WINDOWS_PER_LEN) {
          starts <- sample(starts, MAX_WINDOWS_PER_LEN)
        }
        
        for (start_idx in starts) {
          cal_idx  <- start_idx : (start_idx + calib_len - 1L)
          test_idx <- (start_idx + calib_len) : (start_idx + calib_len + T1 - 1L)
          
          s_cal  <- rolling_country_df[[ "state" ]][cal_idx]
          s_test <- rolling_country_df[[ "state" ]][test_idx]
          
          mm_seq_calib <- mk_pairs(s_cal)
          mm_seq_test  <- mk_pairs(s_test)
          
          p_val <- get_p_val(mm_seq_calib, mm_seq_test, 1000, 4, 4, add_zero, "mm")
          if (p_val > alpha_val) num_included <- num_included + 1L
          num_tests <- num_tests + 1L
          
          print(
            paste(
              "FINISHED Start Index:", start_idx, 
              "P-value:", p_val,
              "Alpha:", alpha_val,
              "Test Length:", T1
            )
          )
        }
      }
      
      empirical_cov <- if (num_tests > 0L) num_included / num_tests else NA_real_
      
      rel_curves_results <- rel_curves_results %>%
        add_row(
          sampling_type      = "rolling",
          test_length        = T1,
          target_coverage    = 1 - alpha_val,
          empirical_coverage = empirical_cov
        )
    }
  }
  
  # --- Random sampling panel ---
  unique_cids <- unique(data$country_id)
  
  for (alpha_val in alpha_grid) {
    for (T1 in test_lengths) {
      
      num_included <- 0L
      num_tests    <- 0L
      
      # Do exactly N_TRIALS_RANDOM random trials per (alpha, T1)
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
        
        p_val <- get_p_val(mm_seq_calib, mm_seq_test, 1000, 4, 4, add_zero, "mm")
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
          sampling_type      = "random",
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

saveRDS(rel_curves_results1, file.path(data_dir, "reliability_curves_data_no_added_zero.rds"))
saveRDS(rel_curves_results2, file.path(data_dir, "reliability_curves_data_added_zero.rds"))