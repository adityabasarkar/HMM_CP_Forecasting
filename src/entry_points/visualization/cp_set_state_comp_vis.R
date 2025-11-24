library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(scales)

data_dir = "outputs/r_objects"
cp_set_results1 = readRDS(file.path(data_dir, "cp_set_results_no_added_zero_fatality_is.rds"))
cp_set_results2 = readRDS(file.path(data_dir, "cp_set_results_added_zero_fatality_is.rds"))
cp_set_results3 = readRDS(file.path(data_dir, "cp_set_results_no_added_zero_fatality_oos.rds"))
cp_set_results4 = readRDS(file.path(data_dir, "cp_set_results_added_zero_fatality_oos.rds"))
cp_set_results5 = readRDS(file.path(data_dir, "likelihood_sum_set_results_is.rds"))
cp_set_results6 = readRDS(file.path(data_dir, "likelihood_sum_set_results_oos.rds"))
cp_set_results7 = readRDS(file.path(data_dir, "cp_set_results_mc.rds"))

T1 = 6

cp_results_to_df <- function(cp_results) {
  bind_rows(lapply(names(cp_results), function(country_name) {
    mat <- cp_results[[country_name]]
    
    col_info <- tibble(colname = colnames(mat)) %>%
      mutate(
        # always grab tp if present
        tp = as.integer(stringr::str_match(colname, "tp=(\\d+)")[, 2]),
        # TS may or may not be present; returns NA if not found
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
          labels = paste0("M", sort(unique(tp)))
        ),
        State = factor(State,
                       levels = c("State 1", "State 2", "State 3", "State 4"))
      )
    
    mat_long
  }))
}

plot_cp_results <- function(df, main_title, subtitle = NULL) {
  
  # state names + colors
  state_labels <- c(
    "State 1" = "Peaceful (P)",
    "State 2" = "Escalation (E)",
    "State 3" = "War (W)",
    "State 4" = "De-escalation (D)"
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
        TS == 1L ~ "P",
        TS == 2L ~ "E",
        TS == 3L ~ "W",
        TS == 4L ~ "D"
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
    facet_wrap(~ Country, scales = "free_x", ncol = 3) +
    labs(
      title    = main_title,
      subtitle = subtitle,
      x        = NULL,
      y        = "Percentage of Candidates"
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

df1 <- cp_results_to_df(cp_set_results1)
df2 <- cp_results_to_df(cp_set_results2)
df3 <- cp_results_to_df(cp_set_results3)
df4 <- cp_results_to_df(cp_set_results4)
df5 <- cp_results_to_df(cp_set_results5)
df6 <- cp_results_to_df(cp_set_results6)
df7 <- cp_results_to_df(cp_set_results7)


p1 <- plot_cp_results(
  df1,
  "State Composition Across Prediction Horizon for CP Sets",
  "Fatality Data, No Added Zero, In-sample, alpha = 0.01"
)

p2 <- plot_cp_results(
  df2,
  "State Composition Across Prediction Horizon for CP Sets",
  "Fatality Data, Added Zero, In-sample, alpha = 0.01"
)

p3 <- plot_cp_results(
  df3,
  "State Composition Across Prediction Horizon for CP Sets",
  "Fatality Data, No Added Zero, Out-of-sample, alpha = 0.01"
)

p4 <- plot_cp_results(
  df4,
  "State Composition Across Prediction Horizon for CP Sets",
  "Fatality Data, Added Zero, Out-of-sample, alpha = 0.01"
)

p5 <- plot_cp_results(
  df5,
  "State Composition Across Prediction Horizon for Likelihood Sets",
  "Fatality Data, In-sample, alpha = 0.01"
)

p6 <- plot_cp_results(
  df6,
  "State Composition Across Prediction Horizon for Likelihood Sets",
  "Fatality Data, Out-of-sample, alpha = 0.01"
)
p7 <- plot_cp_results(
  df7,
  "State Composition Across Prediction Horizon for CP Sets",
  "Monte Carlo Simulations, alpha = 0.01"
)

# View in RStudio
p1; p2; p3; p4; p5; p6; p7

# Save to PDF (tweak width/height if you want)
ggsave("outputs/plots/cp_sets_no_added_zero_is.pdf",   p1, width = 10, height = 11, units = "in", device = cairo_pdf)
ggsave("outputs/plots/cp_sets_added_zero_is.pdf",      p2, width = 10, height = 11, units = "in", device = cairo_pdf)
ggsave("outputs/plots/cp_sets_no_added_zero_oos.pdf",  p3, width = 10, height = 11, units = "in", device = cairo_pdf)
ggsave("outputs/plots/cp_sets_added_zero_oos.pdf",     p4, width = 10, height = 11, units = "in", device = cairo_pdf)
ggsave("outputs/plots/likelihood_sets_is.pdf",     p5, width = 10, height = 11, units = "in", device = cairo_pdf)
ggsave("outputs/plots/likelihood_sets_oos.pdf",     p6, width = 10, height = 11, units = "in", device = cairo_pdf)
ggsave("outputs/plots/monte_carlo_cp_sets.pdf",     p7, width = 10, height = 6, units = "in", device = cairo_pdf)
