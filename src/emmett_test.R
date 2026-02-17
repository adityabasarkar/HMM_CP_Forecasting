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

sweden = data_full[data_full$name == "Sweden", , drop=F]

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

# Functions to run all analyses ------------------------------------------------

mk_pairs <- function(states_vec) lapply(states_vec, function(s) s-1)

get_cp_set_prop = function(end_id, add_zero, ctr_id, data) {
  
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


# Getting all of the results ---------------------------------------------------

prop1 <- get_cp_set_prop(end1, FALSE, 66, data2); print(paste0("CP full --> ", prop1$cardinality_cp_set))
prop2 <- get_likelihood_set_prop(end1, 66, data2); print(paste0("Likelihood full --> ", prop2$cardinality_cp_set))
prop3 <- get_cp_set_prop(end1, TRUE, 66, data2); print(paste0("CP+0 full --> ", prop3$cardinality_cp_set))

# Getting all of the results with less training data ---------------------------

sweden = data2[data2$country_id == 66, ]
sweden_small = sweden[sweden$month_id >= 535, ]

prop1_small <- get_cp_set_prop(end1, FALSE, 66, sweden_small); print(paste0("CP small --> ", prop1_small$cardinality_cp_set))
prop2_small <- get_likelihood_set_prop(end1, 66, sweden_small); print(paste0("Likelihood small --> ", prop2_small$cardinality_cp_set))
prop3_small <- get_cp_set_prop(end1, TRUE, 66, sweden_small); print(paste0("CP+0 small --> ", prop3_small$cardinality_cp_set))

# Plotting the results ---------------------------------------------------------

df_all_train = data.frame("State" = rep(c(rep("State 1", T1), rep("State 2", T1), rep("State 3", T1), rep("State 4", T1)), 3), 
                          "colname" = c(rep(colnames(prop1$ctr_prop), 4), rep(colnames(prop2$ctr_prop), 4), rep(colnames(prop3$ctr_prop), 4)),
                          "Proportion" = c(c(t(prop1$ctr_prop)), c(t(prop2$ctr_prop)), c(t(prop3$ctr_prop))),
                          "tp" = rep(1:T1, 4 * 3),
                          "TS" = rep(1, T1 * 4 * 3),
                          "Country" = c(rep("Conformal Prediction", T1 * 4), rep("Likelihood-Based Prediction", T1 * 4), rep("Conformal Prediction + 0", T1*4)),
                          "tp_factor" = rep(paste0("T + ", 1:T1), 4 * 3))
df_all_train$Country <- factor(df_all_train$Country, 
                   levels = unique(df_all_train$Country))

df_small = data.frame("State" = rep(c(rep("State 1", T1), rep("State 2", T1), rep("State 3", T1), rep("State 4", T1)), 3), 
                      "colname" = c(rep(colnames(prop1_small$ctr_prop), 4), rep(colnames(prop2_small$ctr_prop), 4), rep(colnames(prop3_small$ctr_prop), 4)),
                      "Proportion" = c(c(t(prop1_small$ctr_prop)), c(t(prop2_small$ctr_prop)), c(t(prop3_small$ctr_prop))),
                      "tp" = rep(1:T1, 4 * 3),
                      "TS" = rep(1, T1 * 4 * 3),
                      "Country" = c(rep("Conformal Prediction", T1 * 4), rep("Likelihood-Based Prediction", T1 * 4), rep("Conformal Prediction + 0", T1*4)),
                      "tp_factor" = rep(paste0("T + ", 1:T1), 4 * 3))
df_small$Country <- factor(df_small$Country, 
                               levels = unique(df_small$Country))


plot_cp_results <- function(df, main_title, subtitle = NULL) {
  
  # state names + colors
  state_labels <- c(
    "State 1" = "Peaceful (1)",
    "State 2" = "Escalation (2)",
    "State 3" = "War (3)",
    "State 4" = "De-escalation (4)"
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
  "Real Conflict Data, alpha = 0.2"
)

p5 <- plot_cp_results(
  df_small,
  "Prediction Sets of Forecasted Conflict State-Sequences (T = 6)",
  "Real Conflict Data, alpha = 0.2"
)

ggsave("outputs/plots/limitation_large.pdf", p4, width = 12, height = 6, units = "in", device = "pdf")
ggsave("outputs/plots/limitation_small.pdf", p5, width = 12, height = 6, units = "in", device = "pdf")
