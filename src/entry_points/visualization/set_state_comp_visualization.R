library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(scales)

data_dir = "outputs/r_objects"
plot_dir = "outputs/plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read RDS objects (these now contain results + summary) ----
cp_obj  <- readRDS(file.path(data_dir, "fatality_cp_compositions.rds"))
lik_obj <- readRDS(file.path(data_dir, "fatality_likelihood_compositions.rds"))
cp_obj_12  <- readRDS(file.path(data_dir, "fatality_cp_compositions_12.rds"))
lik_obj_12 <- readRDS(file.path(data_dir, "fatality_likelihood_compositions_12.rds"))

sim_obj <- readRDS(file.path(data_dir, "sim_compare_compositions.rds"))

T1 = 6

# ---- Extract matrices + summaries safely ----
cp_mats  <- if (is.list(cp_obj)  && "cp_set_results" %in% names(cp_obj))  cp_obj$cp_set_results else cp_obj
cp_sum   <- if (is.list(cp_obj)  && "cp_set_summary" %in% names(cp_obj))  cp_obj$cp_set_summary else NULL

lik_mats <- if (is.list(lik_obj) && "likelihood_set_results" %in% names(lik_obj)) lik_obj$likelihood_set_results else lik_obj
lik_sum  <- if (is.list(lik_obj) && "likelihood_set_summary" %in% names(lik_obj)) lik_obj$likelihood_set_summary else NULL

cp_mats_12  <- if (is.list(cp_obj_12)  && "cp_set_results" %in% names(cp_obj_12))  cp_obj_12$cp_set_results else cp_obj_12
cp_sum_12   <- if (is.list(cp_obj_12)  && "cp_set_summary" %in% names(cp_obj_12))  cp_obj_12$cp_set_summary else NULL

lik_mats_12 <- if (is.list(lik_obj_12) && "likelihood_set_results" %in% names(lik_obj_12)) lik_obj_12$likelihood_set_results else lik_obj_12
lik_sum_12  <- if (is.list(lik_obj_12) && "likelihood_set_summary" %in% names(lik_obj_12)) lik_obj_12$likelihood_set_summary else NULL

# sim may have a separate summary file if you saved it; try to read it if present
sim_sum_path <- file.path(data_dir, "sim_compare_set_summary_table.rds")
sim_sum <- if (file.exists(sim_sum_path)) readRDS(sim_sum_path) else NULL

# If sim_obj is wrapped like others, unwrap it too
sim_mats <- if (is.list(sim_obj) && "cp_set_results" %in% names(sim_obj)) sim_obj$cp_set_results else sim_obj


cp_results_to_df <- function(cp_results) {
  bind_rows(lapply(names(cp_results), function(country_name) {
    mat <- cp_results[[country_name]]
    
    col_info <- tibble(colname = colnames(mat)) %>%
      mutate(
        tp = as.integer(stringr::str_match(colname, "tp=(\\d+)")[, 2]),
        TS = as.integer(stringr::str_match(colname, "TS=(\\d+)")[, 2])
      )
    
    mat_df <- as.data.frame(mat) %>%
      mutate(State = rownames(mat))
    
    mat_long <- mat_df %>%
      tidyr::pivot_longer(
        cols = -State,
        names_to = "colname",
        values_to = "Proportion"
      ) %>%
      left_join(col_info, by = "colname") %>%
      mutate(
        Country = country_name,
        tp_factor = factor(
          tp,
          levels = sort(unique(tp)),
          labels = ifelse(sort(unique(tp)) == 0, "T", paste0("T+", sort(unique(tp))))
        ),
        State = factor(State,
                       levels = c("State 1", "State 2", "State 3", "State 4"))
      )
    
    mat_long
  }))
}

plot_cp_results <- function(df, main_title, subtitle = NULL, summary_df = NULL) {
  
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
  
  # ---- Build per-facet footer text ----
  footer_labels <- NULL
  if (!is.null(summary_df)) {
    if (all(c("country_name", "n_cp_set", "true_in_cp") %in% names(summary_df))) {
      footer_labels <- summary_df %>%
        transmute(Country = country_name,
                  footer = paste0("n = ", n_cp_set, " | true_in = ", true_in_cp))
    } else if (all(c("country_name", "n_lik_set", "true_in_lik") %in% names(summary_df))) {
      footer_labels <- summary_df %>%
        transmute(Country = country_name,
                  footer = paste0("n = ", n_lik_set, " | true_in = ", true_in_lik))
    } else if (all(c("set_type", "num", "true_in") %in% names(summary_df))) {
      footer_labels <- summary_df %>%
        transmute(Country = set_type,
                  footer = paste0("n = ", num, " | true_in = ", true_in))
    }
  }
  
  # ---- Create facet label that includes footer on line 2 ----
  df2 <- df
  if (!is.null(footer_labels)) {
    df2 <- df2 %>%
      left_join(footer_labels, by = "Country") %>%
      mutate(Facet = ifelse(is.na(footer) | footer == "",
                            as.character(Country),
                            paste0(Country, "\n", footer)))
  } else {
    df2 <- df2 %>% mutate(Facet = as.character(Country))
  }
  df2$Facet <- factor(df2$Facet, levels = unique(df2$Facet))
  
  # TS labels need the same Facet key
  ts_labels <- df2 %>%
    dplyr::distinct(Facet, tp, tp_factor, TS) %>%
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
  
  p <- ggplot(df2, aes(x = tp_factor, y = Proportion, fill = State)) +
    geom_col(width = 0.8) +
    geom_text(
      data = ts_labels,
      aes(x = tp_factor, y = 1.02, label = TS_label, color = TS_state),
      inherit.aes = FALSE,
      vjust = 0,
      size = 3.0,
      fontface = "bold"
    ) +
    scale_y_continuous(
      limits = c(0, 1.10),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(
      values = state_colors,
      labels = state_labels,
      name   = "State"
    ) +
    scale_color_manual(values = state_colors, guide = "none") +
    facet_wrap(~ Facet, scales = "free_x", ncol = 5, strip.position = "bottom") +
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
      
      # put strips outside so they sit under x-axis, not on panel
      strip.placement   = "outside",
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = 10, lineheight = 0.95),
      
      axis.text.x      = element_text(size = 9, margin = margin(t = 4)),
      axis.text.y      = element_text(size = 9),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      panel.spacing    = unit(1.1, "lines"),
      plot.margin      = margin(t = 24, r = 12, b = 12, l = 12)
    )
  
  p
}


# ---- Build long dfs from the matrix lists ----
df1 <- cp_results_to_df(cp_mats)
df2 <- cp_results_to_df(lik_mats)
df3 <- cp_results_to_df(cp_mats_12)
df4 <- cp_results_to_df(lik_mats_12)
df5 <- cp_results_to_df(sim_mats)

# ---- Plot with summaries passed in ----

# Fatality Data Visualizations
p1 <- plot_cp_results(
  df1,
  "State Composition Across Prediction Horizon for CP Sets",
  "Conflict State Data, T1 = 6, alpha = 0.2",
  summary_df = NULL
)

p2 <- plot_cp_results(
  df2,
  "State Composition Across Prediction Horizon for Likelihood Sets",
  "Conflict State Data, T1 = 6, alpha = 0.2",
  summary_df = NULL
)

p3 <- plot_cp_results(
  df3,
  "State Composition Across Prediction Horizon for CP Sets",
  "Conflict State Data, T1 = 12, alpha = 0.2",
  summary_df = NULL
)

p4 <- plot_cp_results(
  df4,
  "State Composition Across Prediction Horizon for Likelihood Sets",
  "Conflict State Data, T1 = 12, alpha = 0.2",
  summary_df = NULL   # will show if sim_compare_summary.rds exists, otherwise no footer
)

# Simulated Data Visualizations
p5 <- plot_cp_results(
  df5,
  "Comparison of State Compositions between CP Sets and Likelihood Sets",
  "Simulated Data, alpha = 0.2",
  summary_df = NULL   # will show if sim_compare_summary.rds exists, otherwise no footer
)

# Save to PDF
ggsave(file.path(plot_dir, "fatality_cp_compositions.pdf"),         p1, width = 12, height = 5, units = "in", device = cairo_pdf)
ggsave(file.path(plot_dir, "fatality_likelihood_compositions.pdf"), p2, width = 12, height = 5, units = "in", device = cairo_pdf)
ggsave(file.path(plot_dir, "fatality_cp_compositions_12.pdf"),         p3, width = 12, height = 5, units = "in", device = cairo_pdf)
ggsave(file.path(plot_dir, "fatality_likelihood_compositions_12.pdf"), p4, width = 12, height = 5, units = "in", device = cairo_pdf)
ggsave(file.path(plot_dir, "sim_compare_compositions.pdf"),         p5, width = 10, height = 5, units = "in", device = cairo_pdf)
