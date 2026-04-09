library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

data_dir <- "outputs/r_objects"

sim_val_summary_table <- readRDS(file.path(data_dir, "sim_val_summary_table.rds"))

format_num <- function(x) {
  format(round(x, 2), nsmall = 2, trim = TRUE)
}

cp_table <- sim_val_summary_table %>%
  filter(grepl("^CP", Method)) %>%
  mutate(
    Q1 = format_num(Q1),
    Mean = format_num(Mean),
    Median = format_num(Median),
    Q3 = format_num(Q3)
  ) %>%
  select(Method, Q1, Mean, Median, Q3)

lik_table <- sim_val_summary_table %>%
  filter(grepl("^Likelihood", Method)) %>%
  mutate(
    Q1 = format_num(Q1),
    Mean = format_num(Mean),
    Median = format_num(Median),
    Q3 = format_num(Q3)
  ) %>%
  select(Method, Q1, Mean, Median, Q3)

make_table_grob <- function(tbl) {
  grob <- tableGrob(
    tbl,
    rows = NULL,
    theme = ttheme_minimal(
      core = list(
        fg_params = list(fontfamily = "serif", fontsize = 10, hjust = c(0, 1, 1, 1, 1)),
        bg_params = list(fill = "white", col = "white")
      ),
      colhead = list(
        fg_params = list(fontfamily = "serif", fontface = "bold", fontsize = 11),
        bg_params = list(fill = "white", col = "white")
      )
    )
  )
  grob$widths <- unit(c(2.3, 0.9, 1.0, 1.0, 0.9), "in")
  grob
}

cp_grob <- make_table_grob(cp_table)
lik_grob <- make_table_grob(lik_table)

caption_text <- paste(
  "Summary of the number of times u ~ U(0,1) is sampled for the CP approach (left)",
  "compared to the number of times the proposed state-sequences have the same",
  "probability mass for the likelihood-based approach (right).",
  sep = "\n"
)

caption_grob <- arrangeGrob(
  textGrob(
    "TABLE 1",
    x = 0.01,
    y = 0.8,
    just = c("left", "center"),
    gp = gpar(fontfamily = "serif", fontsize = 18, fontface = "bold")
  ),
  textGrob(
    caption_text,
    x = 0.13,
    y = 0.8,
    just = c("left", "center"),
    gp = gpar(fontfamily = "serif", fontsize = 10)
  )
)

table_plot <- arrangeGrob(
  cp_grob,
  lik_grob,
  caption_grob,
  nrow = 2,
  ncol = 2,
  heights = c(1, 0.18),
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2), c(3, 3))
)

ggsave(
  filename = "sim_val_analysis_table.pdf",
  plot = table_plot,
  path = "outputs/plots",
  width = 16,
  height = 5.3,
  units = "in",
  device = cairo_pdf
)
