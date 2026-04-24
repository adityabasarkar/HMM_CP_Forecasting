library(ggplot2)
library(dplyr)

data_dir = "outputs/r_objects"

mc_rel_curve_cp_set_data = readRDS(file.path(data_dir, "mc_rel_curve_cp_set_data.rds"))
mc_rel_curve_likelihood_data = readRDS(file.path(data_dir, "mc_rel_curve_likelihood_data.rds"))
# mc_rel_curve_added_zero_data = readRDS(file.path(data_dir, "mc_rel_curve_added_zero_data.rds"))

## Make sure test_length is treated as a factor for coloring
mc_rel_curve_cp_set_data <- mc_rel_curve_cp_set_data %>%
  mutate(test_length = factor(test_length))

mc_rel_curve_likelihood_data <- mc_rel_curve_likelihood_data %>%
  mutate(test_length = factor(test_length))

# mc_rel_curve_added_zero_data <- mc_rel_curve_added_zero_data %>%
#   mutate(
#     test_length = factor(test_length),
#     # Clean up type labels for nicer facets
#     type = case_when(
#       tolower(type) == "anchored"   ~ "Anchored",
#       tolower(type) == "unanchored" ~ "Unanchored",
#       TRUE                          ~ as.character(type)
#     ),
#     type = factor(type, levels = c("Unanchored", "Anchored"))
#   )

########################################
## 1) CP sets reliability plot (single panel)
########################################

p_cp <- ggplot(
  mc_rel_curve_cp_set_data,
  aes(x = target_coverage, y = empirical_coverage,
      color = test_length, group = test_length)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point() +
  labs(
    x = "Target coverage",
    y = "Empirical coverage",
    color = "Test length",
    title = "Reliability curves – CP sets",
    subtitle = "Monte Carlo Simulations"
  ) +
  coord_equal(xlim = c(0.5, 1), ylim = c(0.5, 1)) +
  theme_minimal()

########################################
## 2) Likelihood-based set reliability plot (single panel)
########################################

p_lik <- ggplot(
  mc_rel_curve_likelihood_data,
  aes(x = target_coverage, y = empirical_coverage,
      color = test_length, group = test_length)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point() +
  labs(
    x = "Target coverage",
    y = "Empirical coverage",
    color = "Test length",
    title = "Reliability curves – Likelihood-based sets",
    subtitle = "Monte Carlo Simulations"
  ) +
  coord_equal(xlim = c(0.5, 1), ylim = c(0.5, 1)) +
  theme_minimal()

########################################
## 3) Added-zero CP (anchored vs unanchored) – 2-panel plot
########################################

# p_zero <- ggplot(
#   mc_rel_curve_added_zero_data,
#   aes(x = target_coverage, y = empirical_coverage,
#       color = test_length, group = test_length)
# ) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~ type, nrow = 1) +
#   labs(
#     x = "Target coverage",
#     y = "Empirical coverage",
#     color = "Test length",
#     title = "Reliability curves – Added-zero CP (Anchored vs Unanchored)",
#     subtitle = "Monte Carlo Simulations"
#   ) +
#   coord_equal(xlim = c(0.5, 1), ylim = c(0.5, 1)) +
#   theme_minimal()

ggsave("outputs/plots/rel_curve_cp_mc.pdf",   p_cp, width = 5, height = 5, units = "in", device = cairo_pdf)
ggsave("outputs/plots/rel_curve_lik_mc.pdf",      p_lik, width = 5, height = 5, units = "in", device = cairo_pdf)
# ggsave("outputs/plots/rel_curve_anchor_mc.pdf",  p_zero, width = 10, height = 5, units = "in", device = cairo_pdf)
