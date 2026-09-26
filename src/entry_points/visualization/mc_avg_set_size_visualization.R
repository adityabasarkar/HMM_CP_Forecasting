library(ggplot2)
library(dplyr)
library(tidyr)

data_dir <- "outputs/r_objects"
plot_dir <- "outputs/plots"

in_path <- file.path(data_dir, "mc_avg_set_sizes.rds")
avg_set_sizes <- readRDS(in_path)

plot_df <- avg_set_sizes %>%
  mutate(
    test_length = factor(test_length,
                         levels = sort(unique(test_length)),
                         labels = paste0("Prediction length ", sort(unique(test_length))))
  ) %>%
  select(test_length, alpha, avg_set_size_cp, avg_set_size_lik, avg_set_size_naive, avg_set_size_weighted) %>%
  pivot_longer(
    cols = c(avg_set_size_cp, avg_set_size_lik, avg_set_size_naive, avg_set_size_weighted),
    names_to = "method",
    values_to = "avg_set_size"
  ) %>%
  mutate(
    method = case_when(
      method == "avg_set_size_cp"  ~ "CP set",
      method == "avg_set_size_lik" ~ "Likelihood set",
      method == "avg_set_size_naive" ~ "Naive CP",
      method == "avg_set_size_weighted" ~ "Weighted CP",
      TRUE ~ method
    ),
    method = factor(method, levels = c("CP set", "Likelihood set", "Naive CP", "Weighted CP"))
  )

p <- ggplot(
  plot_df,
  aes(x = alpha, y = avg_set_size, color = method, group = method)
) +
  geom_line() +
  geom_point() +
  facet_wrap(~ test_length, scales = "free_y") +
  labs(
    x = "Alpha",
    y = "Average prediction set size",
    color = NULL,
    title = "Average prediction set size vs alpha",
    subtitle = "Monte Carlo Simulations"
  ) +
  theme_minimal()

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

ggsave(
  filename = file.path(plot_dir, "avg_set_size_mc_curves.pdf"),
  plot = p,
  width = 10,
  height = 5,
  units = "in",
  device = cairo_pdf
)

print("Saved: outputs/plots/avg_set_size_mc_curves.pdf")
