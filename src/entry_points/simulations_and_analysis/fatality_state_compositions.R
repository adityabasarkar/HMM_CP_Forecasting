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

# countries <- c(70, 235, 149, 57, 245)
countries <- c(70, 149, 57, 28)
# countries <- c(133, 220)

end1 = 540 # Dec 2022
end2 = 547
alpha = 0.2
T1 = 6
set.seed(2569)


mk_pairs <- function(states_vec) lapply(states_vec, function(s) s-1)
get_cp_set_prop = function(end_id, add_zero, ctr_id, data, coinflip = TRUE) {
  
  set.seed(2569)
  country_df = data[data$country_id == ctr_id,]
  
  end_idx <- which(country_df$month_id == end_id)
  cal_idx <- 1:end_idx
  s_cal  <- country_df[[ "state" ]][cal_idx]
  mm_seq_calib <- mk_pairs(s_cal)
  
  # state at time T+0, encoded the same way as cp_set states (0:3)
  start_state <- country_df[["state"]][end_idx] - 1L
  
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
  cp_set <- get_cp_set(mm_seq_calib, c(0), T1, 1000, 4, 4, alpha, add_zero, "mm", allowable_transitions, coinflip)
  
  # prepend T+0 to each sequence
  cp_set <- lapply(cp_set, function(seq) c(start_state, seq))
  
  ctr_res <- matrix(0, nrow = 4, ncol = T1 + 1)        # set all entries to 0  
  rownames(ctr_res) <- paste("State", 1:4)             # give row-names  
  
  true_idx <- (end_idx + 1L):(end_idx + T1)
  if (max(true_idx) > nrow(country_df)) {
    colnames(ctr_res) <- c("tp=0", paste0("tp=", 1:T1))
  } else {
    true_states <- country_df[["state"]][true_idx]
    colnames(ctr_res) <- c(
      paste0("tp=0, TS=", country_df[["state"]][end_idx]),
      paste0("tp=", 1:T1, ", TS=", true_states)
    )
  }
  
  for (seq in cp_set) {
    for (t in seq_along(seq)) {
      state <- seq[t]
      ctr_res[state+1, t] <- ctr_res[state+1, t] + 1
    }
  }
  
  ctr_prop <- sweep(ctr_res, 2, colSums(ctr_res), "/")
  
  n_cp_set <- length(cp_set)
  
  true_in_cp <- FALSE
  true_states <- c(start_state, as.integer(mk_pairs(country_df[["state"]][true_idx])))
  if (!is.null(true_states)) {
    true_in_cp <- any(sapply(cp_set, function(seq) identical(seq, true_states)))
  }
  
  return(list(
    ctr_prop   = ctr_prop,
    n_cp_set   = n_cp_set,
    true_in_cp = true_in_cp
  ))
  
}

get_likelihood_set_prop = function(end_id, ctr_id, data) {
  
  set.seed(2569)
  country_df = data[data$country_id == ctr_id,]
  
  end_idx <- which(country_df$month_id == end_id)
  cal_idx <- 1:end_idx
  s_cal  <- country_df[[ "state" ]][cal_idx]
  mm_seq_calib <- mk_pairs(s_cal)
  
  transition_probs = estimate_transition_probabilities(mm_seq_calib, 4)
  
  # state at time T+0, encoded the same way as cp_set states (0:3)
  start_state <- country_df[["state"]][end_idx] - 1L
  
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
  ord <- order(-probs, runif(length(probs)), method = "radix")
  probs_ord <- probs[ord]
  
  # 3) take top few so that cumulative >= 1 - alpha (prefer overshoot)
  target <- 1 - alpha
  cumul <- cumsum(probs_ord)
  k <- which(cumul >= target)[1]
  if (is.na(k)) k <- length(probs_ord)  # if total < target due to zeros, take all
  
  # selected candidate set, now including T+0 as the first state
  cp_set <- lapply(seq_len(k), function(j) c(start_state, as.integer(grid_df[ord[j], ])))
  
  ctr_res <- matrix(0, nrow = 4, ncol = T1 + 1)        # set all entries to 0  
  rownames(ctr_res) <- paste("State", 1:4)             # give row-names  
  
  true_idx <- (end_idx + 1L):(end_idx + T1)
  if (max(true_idx) > nrow(country_df)) {
    colnames(ctr_res) <- c("tp=0", paste0("tp=", 1:T1))
  } else {
    true_states <- country_df[["state"]][true_idx]
    colnames(ctr_res) <- c(
      paste0("tp=0, TS=", country_df[["state"]][end_idx]),
      paste0("tp=", 1:T1, ", TS=", true_states)
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
  
  n_lik_set <- length(cp_set)
  
  true_in_lik <- FALSE
  true_states <- c(start_state, as.integer(mk_pairs(country_df[["state"]][true_idx])))
  if (!is.null(true_states)) {
    true_in_lik <- any(sapply(cp_set, function(seq) identical(seq, true_states)))
  }
  
  return(list(
    ctr_prop    = ctr_prop,
    n_lik_set   = n_lik_set,
    true_in_lik = true_in_lik
  ))
}



get_cp_set_mat = function(end_id, add_zero) {
  
  cp_set_results = list()
  summary_rows <- list()   # will become the per-country summary table
  
  for (country in countries) {
    
    country_df = data1[data1$country_id == country,]
    country_name = country_df[[1, "name"]]
    
    # get CP results + diagnostics
    res_cp <- get_cp_set_prop(end_id, add_zero, country, data1)
    
    # 1) keep the proportions matrix as before
    cp_set_results[[country_name]] <- res_cp$ctr_prop
    
    # 2) store per-country summary row
    summary_rows[[country_name]] <- tibble::tibble(
      country_id   = country,
      country_name = country_name,
      n_cp_set     = res_cp$n_cp_set,
      true_in_cp   = res_cp$true_in_cp
    )
    
    print(paste("FINISHED Country:", country_name))
  }
  
  # bind all rows into one table
  cp_set_summary <- dplyr::bind_rows(summary_rows)
  
  # return both objects
  return(list(
    cp_set_results = cp_set_results,
    cp_set_summary = cp_set_summary
  ))
}

likelihood_set = function(end_id) {
  
  likelihood_set_results = list()
  summary_rows <- list()
  
  for (country in countries) {
    
    country_df = data1[data1$country_id == country,]
    country_name = country_df[[1, "name"]]
    
    res_lik <- get_likelihood_set_prop(end_id, country, data1)
    
    # keep proportions matrix
    likelihood_set_results[[country_name]] <- res_lik$ctr_prop
    
    # per-country summary row
    summary_rows[[country_name]] <- tibble::tibble(
      country_id   = country,
      country_name = country_name,
      n_lik_set    = res_lik$n_lik_set,
      true_in_lik  = res_lik$true_in_lik
    )
    
    print(paste("FINISHED Country:", country_name))
  }
  
  likelihood_set_summary <- dplyr::bind_rows(summary_rows)
  
  return(list(
    likelihood_set_results = likelihood_set_results,
    likelihood_set_summary = likelihood_set_summary
  ))
}

compare_set = function(end_id) {
  
  comp_set_results = list()
  
  country_df = data2[data2$country_id == 66,]
  country_name = country_df[[1, "name"]]
  
  # CP Set
  res1 <- get_cp_set_prop(end_id, FALSE, 66, data2)
  comp_set_results[["CP Set"]] <- res1$ctr_prop
  
  # CP Set & Added P
  res2 <- get_cp_set_prop(end_id, TRUE, 66, data2)
  comp_set_results[["CP Set & Added 1"]] <- res2$ctr_prop
  
  # Likelihood Set
  res3 <- get_likelihood_set_prop(end_id, 66, data2)
  comp_set_results[["Likelihood Set"]] <- res3$ctr_prop
  
  
  
  # one-row summary table for this country
  comp_set_summary <- tibble::tibble(
    set_type = c("CP Set", "CP Set & Added 1", "Likelihood Set"),
    num      = c(res1$n_cp_set, res2$n_cp_set, res3$n_lik_set),
    true_in  = c(res1$true_in_cp, res2$true_in_cp, res3$true_in_lik)
  )
  
  return(list(
    comp_set_results = comp_set_results,
    comp_set_summary = comp_set_summary
  ))
}

data_dir = "outputs/r_objects"

# -----------------------
# CP SET RESULTS
# -----------------------
cp_out_fatality <- get_cp_set_mat(end2, FALSE) # A1

# save the whole object (results + summary)
saveRDS(cp_out_fatality, file.path(data_dir, "fatality_cp_compositions.rds"))

# ALSO save the summary table separately
saveRDS(cp_out_fatality$cp_set_summary,
        file.path(data_dir, "fatality_cp_summary.rds"))

# (optional) save compositions-only separately
saveRDS(cp_out_fatality$cp_set_results,
        file.path(data_dir, "fatality_cp_compositions_only.rds"))


# -----------------------
# LIKELIHOOD SET RESULTS
# -----------------------
lik_out_fatality <- likelihood_set(end2) # A3

# save the whole object (results + summary)
saveRDS(lik_out_fatality, file.path(data_dir, "fatality_likelihood_compositions.rds"))

# ALSO save the summary table separately
saveRDS(lik_out_fatality$likelihood_set_summary,
        file.path(data_dir, "fatality_likelihood_summary.rds"))

# (optional) save compositions-only separately
saveRDS(lik_out_fatality$likelihood_set_results,
        file.path(data_dir, "fatality_likelihood_compositions_only.rds"))


# -----------------------
# COMPARISON SET RESULTS
# -----------------------
comp_out_fatality <- compare_set(end1)

# save the whole object (results + summary)
saveRDS(comp_out_fatality, file.path(data_dir, "fatality_comp_compositions.rds"))

# ALSO save the summary table separately
saveRDS(comp_out_fatality$comp_set_summary,
        file.path(data_dir, "fatality_comp_summary.rds"))

# (optional) save compositions-only separately
saveRDS(comp_out_fatality$comp_set_results,
        file.path(data_dir, "fatality_comp_compositions_only.rds"))