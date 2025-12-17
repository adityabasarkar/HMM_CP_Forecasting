library(arrow)
library(dplyr)
library(tibble)

library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(future)

sourceCpp("src/cpp_files/hmm_conf_pred_final.cpp")

data_dir = "outputs/r_objects"
data = readRDS(file.path(data_dir, "preprocessed_data.rds"))

countries <- c(70, 235, 149, 57, 245)
end1 = 540 # Dec 2022
end2 = 547
alpha = 0.2
T1 = 6
set.seed(2569)


mk_pairs <- function(states_vec) lapply(states_vec, function(s) s-1)
get_cp_set_mat = function(end_id, add_zero) {
  cp_set_results = list()
  i = 1
  for (country in countries) {
    
    country_df = data[data$country_id == country,]
    country_name = country_df[[1, "name"]]
    
    end_idx <- which(country_df$month_id == end_id)
    cal_idx <- 1:end_idx
    s_cal  <- country_df[[ "state" ]][cal_idx]
    mm_seq_calib <- mk_pairs(s_cal)
    
    # cp set 1
    allowable_transitions <- matrix(
      c(1,1,0,0,
        0,0,1,1,
        0,0,1,1,
        1,1,0,0),
      nrow = 4,
      ncol = 4,
      byrow = TRUE
    )
    cp_set <- get_cp_set(mm_seq_calib, c(0), T1, 1000, 4, 4, alpha, add_zero, "mm", allowable_transitions)
    
    ctr_res <- matrix(0, nrow = 4, ncol = T1)        # set all entries to 0  
    rownames(ctr_res) <- paste("State", 1:4)       # give row-names  
    
    true_idx <- (end_idx + 1L):(end_idx + T1)
    if (max(true_idx) > nrow(country_df)) {
      colnames(ctr_res) <- paste0(
        "tp=", 1:T1
      )
    } else {
      true_states <- country_df[["state"]][true_idx]
      colnames(ctr_res) <- paste0(
        "tp=", 1:T1, ", TS=", true_states
      )
    }
      
    
    for (seq in cp_set) {
      for (t in seq_along(seq)) {
        state <- seq[t]
        ctr_res[state+1, t] <- ctr_res[state+1, t] + 1
      }
    }
    
    ctr_prop <- sweep(ctr_res, 2, colSums(ctr_res), "/")
    cp_set_results[[country_name]] <- ctr_prop
    
    i = i + 1
    
    print(
      paste(
        "FINISHED Country:", country_name
      )
    )
  }
  return(cp_set_results)
}

cp_set_results1 <- get_cp_set_mat(end1, FALSE) # A1
# cp_set_results2 <- get_cp_set_mat(end1, TRUE) # A2
# cp_set_results3 <- get_cp_set_mat(end2, FALSE) # B1
# cp_set_results4 <- get_cp_set_mat(end2, TRUE) # B2

data_dir = "outputs/r_objects"
saveRDS(cp_set_results1, file.path(data_dir, "cp_set_results_no_added_zero_fatality_is.rds"))
# saveRDS(cp_set_results2, file.path(data_dir, "cp_set_results_added_zero_fatality_is.rds"))
# saveRDS(cp_set_results3, file.path(data_dir, "cp_set_results_no_added_zero_fatality_oos.rds"))
# saveRDS(cp_set_results4, file.path(data_dir, "cp_set_results_added_zero_fatality_oos.rds"))


# Building set through ordering likelihoods
likelihood_set = function(end_id) {
  i = 1
  likelihood_sum_set_results = list()
  for (country in countries) {
    
    country_df = data[data$country_id == country,]
    country_name = country_df[[1, "name"]]
    
    end_idx <- which(country_df$month_id == end_id)
    cal_idx <- 1:end_idx
    s_cal  <- country_df[[ "state" ]][cal_idx]
    mm_seq_calib <- mk_pairs(s_cal)
    
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
      to   <- states01[-1]               + 1L
      if (any(allowable_transitions[cbind(from, to)] == 0)) {
        return(0)
      }
      prod(transition_probs[cbind(from, to)])
    }
    
    # 1) score all candidate sequences
    probs <- apply(grid_df, 1L, path_prob)
    
    # 2) order by likelihood (desc)
    ord <- order(probs, decreasing = TRUE)
    probs_ord <- probs[ord]
    
    # 3) take top few so that cumulative >= 1 - alpha (prefer overshoot)
    target <- 1 - alpha
    cumul <- cumsum(probs_ord)
    k <- which(cumul >= target)[1]
    if (is.na(k)) k <- length(probs_ord)  # if total < target due to zeros, take all
    
    # selected candidate set (each is a length-T1 future sequence, without the start state)
    cp_set <- lapply(seq_len(k), function(j) as.integer(grid_df[ord[j], ]))
    
    ctr_res <- matrix(0, nrow = 4, ncol = T1)        # set all entries to 0  
    rownames(ctr_res) <- paste("State", 1:4)       # give row-names  
    
    true_idx <- (end_idx + 1L):(end_idx + T1)
    if (max(true_idx) > nrow(country_df)) {
      colnames(ctr_res) <- paste0(
        "tp=", 1:T1
      )
    } else {
      true_states <- country_df[["state"]][true_idx]
      colnames(ctr_res) <- paste0(
        "tp=", 1:T1, ", TS=", true_states
      )
    } 
    
    for (seq in cp_set) {
      for (t in seq_along(seq)) {
        state <- seq[t]
        ctr_res[state+1, t] <- ctr_res[state+1, t] + 1
      }
    }
    
    # normalize to proportions (so each column sums to 1)
    ctr_prop <- sweep(ctr_res, 2, colSums(ctr_res), "/")
    
    likelihood_sum_set_results[[country_name]] <- ctr_prop
    
    i = i + 1
    
    print(
      paste(
        "FINISHED Country:", country_name
      )
    )
    
  }
  return(likelihood_sum_set_results)
}

likelihood_sum_set_results_fatality1 = likelihood_set(end1) # A3
# likelihood_sum_set_results_fatality2 = likelihood_set(end2) # B3

data_dir = "outputs/r_objects"
saveRDS(likelihood_sum_set_results_fatality1, file.path(data_dir, "likelihood_sum_set_results_is.rds"))
# saveRDS(likelihood_sum_set_results_fatality2, file.path(data_dir, "likelihood_sum_set_results_oos.rds"))