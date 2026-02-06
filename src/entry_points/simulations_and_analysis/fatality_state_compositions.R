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
data1 = readRDS(file.path(data_dir, "preprocessed_data.rds"))
data2 = readRDS(file.path(data_dir, "preprocessed_data_full.rds"))

countries <- c(70, 235, 149, 57, 245)
end1 = 540 # Dec 2022
end2 = 547
alpha = 0.2
T1 = 6
set.seed(2569)


mk_pairs <- function(states_vec) lapply(states_vec, function(s) s-1)
get_cp_set_prop = function(end_id, add_zero, ctr_id, data) {
  
  country_df = data[data$country_id == ctr_id,]
  
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
  
  return(ctr_prop)
  
}

get_likelihood_set_prop = function(end_id, ctr_id, data) {
  
  country_df = data[data$country_id == ctr_id,]
  
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
  
  return(ctr_prop)
}





get_cp_set_mat = function(end_id, add_zero) {
  
  cp_set_results = list()
  
  for (country in countries) {
    
    country_df = data1[data1$country_id == country,]
    country_name = country_df[[1, "name"]]
    
    ctr_prop <- get_cp_set_prop(end_id, add_zero, country, data1)
    cp_set_results[[country_name]] <- ctr_prop
    print(
      paste(
        "FINISHED Country:", country_name
      )
    )
    
  }
  
  return(cp_set_results)
  
}

likelihood_set = function(end_id) {
  likelihood_set_results = list()
  
  for (country in countries) {
    
    country_df = data1[data1$country_id == country,]
    country_name = country_df[[1, "name"]]
    
    ctr_prop <- get_likelihood_set_prop(end_id, country, data1)
    likelihood_set_results[[country_name]] <- ctr_prop
    print(
      paste(
        "FINISHED Country:", country_name
      )
    )
    
  }
  
  return(likelihood_set_results)
}

compare_set = function(end_id) {
  
  comp_set_results = list()
  approaches = c("CP Set", "Likelihood Set", "CP Set & Added P")
    
  country_df = data2[data2$country_id == 66,]
  country_name = country_df[[1, "name"]]
  
  prop1 <- get_cp_set_prop(end_id, FALSE, 66, data2)
  comp_set_results[["CP Set"]] <- prop1
  
  prop2 <- get_likelihood_set_prop(end_id, 66, data2)
  comp_set_results[["Likelihood Set"]] <- prop2
  
  prop3 <- get_cp_set_prop(end_id, TRUE, 66, data2)
  comp_set_results[["CP Set & Added P"]] <- prop3

  return(comp_set_results)
  
}

# CP SET RESULTS
cp_set_results_fatality <- get_cp_set_mat(end1, FALSE) # A1
data_dir = "outputs/r_objects"
saveRDS(cp_set_results_fatality, file.path(data_dir, "fatality_cp_compositions.rds"))


# LIKELIHOOD SET RESULTS
likelihood_sum_set_results_fatality = likelihood_set(end1) # A3
data_dir = "outputs/r_objects"
saveRDS(likelihood_sum_set_results_fatality, file.path(data_dir, "fatality_likelihood_compositions.rds"))


# COMPARISON SET RESULTS
comp_set_results_fatality = compare_set(end1)
data_dir = "outputs/r_objects"
saveRDS(comp_set_results_fatality, file.path(data_dir, "fatality_comp_compositions.rds"))
