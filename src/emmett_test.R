library(arrow)
library(dplyr)
library(tibble)

library(ggplot2)
library(tidyr)
library(stringr)
library(scales)

library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(future)

data_clean = readRDS(file.path("outputs/r_objects/preprocessed_data.rds"))
data_full  = readRDS(file.path("outputs/r_objects/preprocessed_data_full.rds"))
data_full  = data_full[!is.na(data_full$name), , drop=F]

sourceCpp("src/cpp_files/hmm_conf_pred_final.cpp")

data_dir = "outputs/r_objects"
data1 = readRDS(file.path(data_dir, "preprocessed_data.rds"))
data2 = readRDS(file.path(data_dir, "preprocessed_data_full.rds"))

countries <- c(70, 235, 149, 57, 245)
end1 = 540 # Dec 2022
end2 = 547
alpha = 0.2
T1 = 6
set.seed(1)

# Functions to run all analyses ------------------------------------------------

mk_pairs <- function(states_vec) lapply(states_vec, function(s) s-1)

get_cp_set_prop = function(end_id, add_zero, ctr_id, data, coinflip = TRUE) {

  country_df = data[data$country_id == ctr_id,]

  end_idx <- which(country_df$month_id == end_id)
  cal_idx <- 1:end_idx
  s_cal  <- country_df[[ "state" ]][cal_idx]
  mm_seq_calib <- mk_pairs(s_cal)

  end_T_T_1 <- which(country_df$month_id == (end_id+T1))
  true_forecast_seq = country_df$state[(end_idx+1):end_T_T_1]
  true_forecast_seq = true_forecast_seq - 1 # reduce state index from 1,2,3,4 to 0,1,2,3

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

  for (seq in cp_set[[1]]) {
    for (t in seq_along(seq)) {
      state <- seq[t]
      ctr_res[state+1, t] <- ctr_res[state+1, t] + 1
    }
  }

  ctr_prop <- sweep(ctr_res, 2, colSums(ctr_res), "/")

  cp_set_matrix = do.call('rbind', cp_set[[1]])
  cp_set_matrix_prop = rbind(apply(cp_set_matrix, 2, function(x) {mean(x == 0)}),
                             apply(cp_set_matrix, 2, function(x) {mean(x == 1)}),
                             apply(cp_set_matrix, 2, function(x) {mean(x == 2)}),
                             apply(cp_set_matrix, 2, function(x) {mean(x == 3)}))
  rownames(cp_set_matrix_prop) <- paste("State", 1:4)       # give row-names

  # Check if true sequence is in prediction set
  contains_true_seq = any(apply(cp_set_matrix, 1, function(row) all(row == true_forecast_seq)))


  cp_set_results = list("ctr_prop" = ctr_prop,
                        "cp_set_matrix_prop" = cp_set_matrix_prop,
                        "cp_set_matrix" = cp_set_matrix,
                        "cardinality_cp_set" = nrow(cp_set_matrix),
                        "contains_true_seq"  = contains_true_seq)

  return(cp_set_results)

}

get_likelihood_set_prop = function(end_id, ctr_id, data) {

  country_df = data[data$country_id == ctr_id,]

  end_idx <- which(country_df$month_id == end_id)
  cal_idx <- 1:end_idx
  s_cal  <- country_df[[ "state" ]][cal_idx]
  mm_seq_calib <- mk_pairs(s_cal)

  end_T_T_1 <- which(country_df$month_id == (end_id+T1))
  true_forecast_seq = country_df$state[(end_idx+1):end_T_T_1]
  true_forecast_seq = true_forecast_seq - 1 # reduce state index from 1,2,3,4 to 0,1,2,3

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

  cp_set_matrix = do.call('rbind', cp_set)
  cp_set_matrix_prop = rbind(apply(cp_set_matrix, 2, function(x) {mean(x == 0)}),
                             apply(cp_set_matrix, 2, function(x) {mean(x == 1)}),
                             apply(cp_set_matrix, 2, function(x) {mean(x == 2)}),
                             apply(cp_set_matrix, 2, function(x) {mean(x == 3)}))
  rownames(cp_set_matrix_prop) <- paste("State", 1:4)       # give row-names

  # Check if true sequence is in prediction set
  contains_true_seq = any(apply(cp_set_matrix, 1, function(row) all(row == true_forecast_seq)))


  cp_set_results = list("ctr_prop" = ctr_prop,
                        "cp_set_matrix_prop" = cp_set_matrix_prop,
                        "cp_set_matrix" = cp_set_matrix,
                        "cardinality_cp_set" = nrow(cp_set_matrix),
                        "contains_true_seq"  = contains_true_seq)

  return(cp_set_results)
}

# ------------------------------------------------------------------------------
# Get all results Section 5 ----------------------------------------------------
# ------------------------------------------------------------------------------
for(coinflip in c(TRUE,FALSE)) {

  print(paste0("Coinflip = ", coinflip))

  # Getting all of the results ---------------------------------------------------

  prop1 <- get_cp_set_prop(end1, FALSE, 66, data2, coinflip); print(paste0("CP full --> ", prop1$cardinality_cp_set))
  print(paste0("CP full X(T+T1) = 0 --> ", sum(prop1$cp_set_matrix[,ncol(prop1$cp_set_matrix)] == 0)))
  prop2 <- get_likelihood_set_prop(end1, 66, data2); print(paste0("Likelihood full --> ", prop2$cardinality_cp_set))
  print(paste0("Likelihood full X(T+T1) = 0 --> ", sum(prop2$cp_set_matrix[,ncol(prop2$cp_set_matrix)] == 0)))
  prop3 <- get_cp_set_prop(end1, TRUE, 66, data2, coinflip); print(paste0("CP+0 full --> ", prop3$cardinality_cp_set))
  print(paste0("CP+0 full X(T+T1) = 0 --> ", sum(prop3$cp_set_matrix[,ncol(prop3$cp_set_matrix)] == 0)))

  # Getting all of the results with less training data ---------------------------

  sweden = data2[data2$country_id == 66, ]
  sweden_small = sweden[sweden$month_id >= 535, ]

  prop1_small <- get_cp_set_prop(end1, FALSE, 66, sweden_small, coinflip); print(paste0("CP small --> ", prop1_small$cardinality_cp_set))
  print(paste0("CP small X(T+T1) = 0 --> ", sum(prop1_small$cp_set_matrix[,ncol(prop1_small$cp_set_matrix)] == 0)))
  prop2_small <- get_likelihood_set_prop(end1, 66, sweden_small); print(paste0("Likelihood small --> ", prop2_small$cardinality_cp_set))
  print(paste0("Likelihood small X(T+T1) = 0 --> ", sum(prop2_small$cp_set_matrix[,ncol(prop2_small$cp_set_matrix)] == 0)))
  prop3_small <- get_cp_set_prop(end1, TRUE, 66, sweden_small, coinflip); print(paste0("CP+0 small --> ", prop3_small$cardinality_cp_set))
  print(paste0("CP+0 small X(T+T1) = 0 --> ", sum(prop3_small$cp_set_matrix[,ncol(prop3_small$cp_set_matrix)] == 0)))
  
  # Plotting the results ---------------------------------------------------------
  
  state_at_T = sweden$state[sweden$month_id == end1]
  prop_at_T = rep(0, 4); prop_at_T[state_at_T] = 1

  df_all_train = data.frame("State" = rep(c(rep("State 1", T1+1), rep("State 2", T1+1), 
                                            rep("State 3", T1+1), rep("State 4", T1+1)), 3),
                            "colname" = c(rep(c(paste0("tp=0, TS=",state_at_T), colnames(prop1$ctr_prop)), 4), 
                                          rep(c(paste0("tp=0, TS=",state_at_T), colnames(prop2$ctr_prop)), 4), 
                                          rep(c(paste0("tp=0, TS=",state_at_T), colnames(prop3$ctr_prop)), 4)),
                            "Proportion" = c(c(rbind(prop_at_T,t(prop1$ctr_prop))), 
                                             c(rbind(prop_at_T,t(prop2$ctr_prop))), 
                                             c(rbind(prop_at_T,t(prop3$ctr_prop)))),
                            "tp" = rep(0:T1, 4 * 3),
                            "TS" = rep(1, (T1+1) * 4 * 3),
                            "Country" = c(rep("Conformal Prediction", (T1+1) * 4), 
                                          rep("Likelihood-Based Prediction", (T1+1) * 4), 
                                          rep("Conformal Prediction + 1", (T1+1)*4)),
                            "tp_factor" = rep(c("T", paste0("T+", 1:T1)), 4 * 3))
  df_all_train$Country <- factor(df_all_train$Country,
                                 levels = unique(df_all_train$Country))

  df_small = data.frame("State" = rep(c(rep("State 1", T1+1), rep("State 2", T1+1), 
                                        rep("State 3", T1+1), rep("State 4", T1+1)), 3),
                        "colname" = c(rep(c(paste0("tp=0, TS=",state_at_T), colnames(prop1_small$ctr_prop)), 4), 
                                      rep(c(paste0("tp=0, TS=",state_at_T), colnames(prop2_small$ctr_prop)), 4), 
                                      rep(c(paste0("tp=0, TS=",state_at_T), colnames(prop3_small$ctr_prop)), 4)),
                        "Proportion" = c(c(rbind(prop_at_T,t(prop1_small$ctr_prop))), 
                                         c(rbind(prop_at_T,t(prop2_small$ctr_prop))), 
                                         c(rbind(prop_at_T,t(prop3_small$ctr_prop)))),
                        "tp" = rep(0:T1, 4 * 3),
                        "TS" = rep(1, (T1+1) * 4 * 3),
                        "Country" = c(rep("Conformal Prediction", (T1+1) * 4), 
                                      rep("Likelihood-Based Prediction", (T1+1) * 4), 
                                      rep("Conformal Prediction + 1", (T1+1)*4)),
                        "tp_factor" = rep(c("T", paste0("T+", 1:T1)), 4 * 3))
  df_small$Country <- factor(df_small$Country,
                             levels = unique(df_small$Country))


  plot_cp_results <- function(df, main_title, subtitle = NULL) {

    # state names + colors
    state_labels <- c(
      "State 1" = "State 1 (Peaceful)",
      "State 2" = "State 2 (Escalation)",
      "State 3" = "State 3 (War)",
      "State 4" = "State 4 (Deescalation)"
    )

    state_colors <- c(
      "State 1" = "skyblue2",
      "State 2" = "darkorange2",
      "State 3" = "firebrick2",
      "State 4" = "seagreen3"
    )

    ts_labels <- df %>%
      dplyr::distinct(Country, tp, tp_factor, TS) %>%
      dplyr::filter(!is.na(TS)) %>%
      dplyr::mutate(
        TS_label = dplyr::case_when(
          TS == 1L ~ "1",
          TS == 2L ~ "2",
          TS == 3L ~ "3",
          TS == 4L ~ "4"
        ),
        TS_state = factor(paste("State", TS),
                          levels = names(state_colors))
      )

    ggplot(df, aes(x = tp_factor, y = Proportion, fill = State)) +
      geom_col(width = 0.8) +
      # TS label just above the bar, bold
      geom_text(
        data = ts_labels,
        aes(x = tp_factor, y = 1.02, label = TS_label, color = TS_state),
        inherit.aes = FALSE,
        vjust = 0,              # anchor at y = 1.02 baseline
        size = 3.0,
        fontface = "bold"
      ) +
      # allow a little extra space above 1.0 for the labels
      scale_y_continuous(
        limits = c(0, 1.1),
        breaks = seq(0, 1, by = 0.25),
        labels = scales::percent_format(accuracy = 1),
        expand = expansion(mult = c(0, 0.02))
      ) +
      # don't clip the TS labels at the top of each panel
      coord_cartesian(clip = "off") +
      scale_fill_manual(
        values = state_colors,
        labels = state_labels,
        name   = "State"
      ) +
      scale_color_manual(
        values = state_colors,
        guide  = "none"
      ) +
      facet_wrap(~ Country, ncol = length(unique(ts_labels$Country))) +
      labs(
        title    = main_title,
        subtitle = subtitle,
        x        = NULL,
        y        = "Percentage"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 11),
        axis.text.x      = element_text(size = 9, margin = margin(t = 4)),
        axis.text.y      = element_text(size = 9),
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold"),
        panel.spacing    = unit(1.1, "lines"),
        plot.margin      = margin(t = 24, r = 12, b = 12, l = 12)
      )
  }

  p4 <- plot_cp_results(
    df_all_train,
    "Prediction Sets of Forecasted Conflict State-Sequences (T = 420)",
    "Conflict State Data, alpha = 0.2"
  )

  p5 <- plot_cp_results(
    df_small,
    "Prediction Sets of Forecasted Conflict State-Sequences (T = 6)",
    "Conflict State Data, alpha = 0.2"
  )

  ggsave(paste0("outputs/plots/limitation_large_", as.numeric(coinflip), ".pdf"), p4, width = 12, height = 4, units = "in", device = "pdf")
  ggsave(paste0("outputs/plots/limitation_small_", as.numeric(coinflip), ".pdf"), p5, width = 12, height = 4, units = "in", device = "pdf")

}

# ------------------------------------------------------------------------------
# comparing CP set sizes -------------------------------------------------------
# ------------------------------------------------------------------------------
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

marg_dist = matrix(c(0,1,0,0),nrow = 1)
error = 0.000000001
error_big = T
while(error_big) {
  curr_dist = marg_dist
  marg_dist = marg_dist %*% P
  
  error_diff = abs(marg_dist - curr_dist)
  if(sum(error_diff < error) == 4) {
    error_big = F
  }
  print(marg_dist)
}

# generate_sequence_mm <- function(P, k, n) {
#   return_sequence <- list()
#   possible_states <- 0:(k - 1)
#   first_state <- sample(possible_states, 1)
#   return_sequence <- append(return_sequence, list(first_state))
#   for (i in 2:n) {
#     previous_state <- return_sequence[[i - 1]]
#     state <- sample(possible_states, 1, replace = FALSE, prob = P[previous_state + 1, ])
#     return_sequence <- append(return_sequence, list(state))
#   }
#   return(return_sequence)
# }
# 
# mc_avg_set_sizes <- function() {
#   
#   results <- tibble(
#     test_length = integer(),
#     target_coverage = numeric(),
#     alpha = numeric(),
#     avg_set_size_cp = numeric(),
#     avg_set_size_lik = numeric()
#   )
#   
#   for (T1 in test_lengths) {
#     
#     grid_df <- expand.grid(
#       rep(list(0:(k - 1L)), T1),
#       KEEP.OUT.ATTRS = FALSE,
#       stringsAsFactors = FALSE
#     )
#     
#     for (cov_target in target_coverage_grid) {
#       
#       alpha_val <- 1 - cov_target
#       
#       cp_sizes <- numeric(N_TRIALS_RANDOM)
#       lik_sizes <- numeric(N_TRIALS_RANDOM)
#       
#       emp_covg_cp = rep(0, N_TRIALS_RANDOM)
#       emp_covg_like = rep(0, N_TRIALS_RANDOM)
#       emp_prob_like = rep(0, N_TRIALS_RANDOM)
#       set.seed(2569)
#       
#       for (i in seq_len(N_TRIALS_RANDOM)) {
#         print(i)
#         
#         mm_seq <- generate_sequence_mm(P, k, calib_len + T1)
#         mm_seq_calib <- mm_seq[1:calib_len]
#         
#         cp_set <- get_cp_set(
#           mm_seq_calib,
#           as.integer(c()),
#           as.integer(T1),
#           as.integer(1000L),
#           as.integer(k),
#           as.integer(k),
#           alpha_val,
#           FALSE,
#           "mm",
#           allowable_transitions,
#           TRUE
#         )
#         cp_sizes[i] <- length(cp_set)
#         
#         # CP: my stuff ---------------------------------------------------------
#         seq_T_T1 = do.call('c', mm_seq)
#         true_forecast_seq = tail(seq_T_T1, T1)
#         cp_set_matrix = do.call('rbind', cp_set)
#         
#         # Check if true sequence is in prediction set
#         if(length(cp_set) > 0) {
#           cp_contains_true_seq = as.numeric(any(apply(cp_set_matrix, 1, function(row) all(row == true_forecast_seq))))
#           emp_covg_cp[i] = cp_contains_true_seq  
#         }
#         # ----------------------------------------------------------------------
#         
#         transition_probs <- estimate_transition_probabilities(mm_seq_calib, 4)
#         start_state <- mm_seq_calib[[length(mm_seq_calib)]][1]
#         
#         path_prob <- function(seq) {
#           states01 <- c(start_state, as.integer(seq))
#           from <- states01[-length(states01)] + 1L
#           to <- states01[-1] + 1L
#           if (any(allowable_transitions[cbind(from, to)] == 0)) {
#             return(0)
#           }
#           prod(transition_probs[cbind(from, to)])
#         }
#         
#         probs <- apply(grid_df, 1L, path_prob)
#         ord <- order(-probs, runif(length(probs)), method = "radix")
#         probs_ord <- probs[ord]
#         
#         cumul <- cumsum(probs_ord)
#         cutoff_idx <- which(cumul >= cov_target)[1]
#         if (is.na(cutoff_idx)) cutoff_idx <- length(probs_ord)
#         lik_sizes[i] <- cutoff_idx
#         
#         # LIKE: my stuff -------------------------------------------------------
#         like_set_matrix_full = grid_df[ord, ,drop=F]
#         like_set_matrix = like_set_matrix_full[1:cutoff_idx, ,drop=F]
#         
#         # Check if true sequence is in prediction set
#         like_contains_true_seq = as.numeric(any(apply(like_set_matrix, 1, function(row) all(row == true_forecast_seq))))
#         emp_covg_like[i] = like_contains_true_seq
#         emp_prob_like[i] = cumul[cutoff_idx]
#         # ----------------------------------------------------------------------
#       }
#       
#       results <- results %>%
#         add_row(
#           test_length = T1,
#           target_coverage = cov_target,
#           alpha = alpha_val,
#           avg_set_size_cp = mean(cp_sizes),
#           avg_set_size_lik = mean(lik_sizes)
#         )
#       
#       print(paste(
#         "T1:", T1,
#         "target_coverage:", cov_target,
#         "alpha:", alpha_val,
#         "avg_cp_size:", results$avg_set_size_cp[nrow(results)],
#         "avg_lik_size:", results$avg_set_size_lik[nrow(results)]
#       ))
#     }
#   }
#   
#   results
# }
# 
# avg_set_sizes <- mc_avg_set_sizes()
# 
# out_path <- file.path(data_dir, "mc_avg_set_sizes.rds")
# saveRDS(avg_set_sizes, out_path)
# print(paste("Saved:", out_path))
# 
# 
# 
# 
# 
