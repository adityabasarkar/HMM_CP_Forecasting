library(ggplot2)
library(dplyr)

data_dir = "outputs/r_objects"

cp_original_rel_data_fatality = readRDS(file.path(data_dir, "cp_original_rel_data_fatality.rds"))
likelihood_rel_data_fatality = readRDS(file.path(data_dir, "likelihood_rel_data_fatality.rds"))
naive_rel_data_fatality = readRDS(file.path(data_dir, "naive_rel_data_fatality.rds"))
weighted_rel_data_fatality = readRDS(file.path(data_dir, "weighted_rel_data_fatality.rds"))
# cp_added_zero_rel_data_fatality = readRDS(file.path(data_dir, "cp_added_zero_rel_data_fatality.rds"))

## Make sure test_length is treated as a factor for coloring
cp_original_rel_data_fatality <- cp_original_rel_data_fatality %>%
  mutate(test_length = factor(test_length))

likelihood_rel_data_fatality <- likelihood_rel_data_fatality %>%
  mutate(test_length = factor(test_length))

naive_rel_data_fatality <- naive_rel_data_fatality %>%
  mutate(test_length = factor(test_length))
weighted_rel_data_fatality <- weighted_rel_data_fatality %>%
  mutate(test_length = factor(test_length))

# cp_added_zero_rel_data_fatality <- cp_added_zero_rel_data_fatality %>%
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
  cp_original_rel_data_fatality,
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
    subtitle = "Conflict State Data"
  ) +
  coord_equal(xlim = c(0.5, 1), ylim = c(0.5, 1)) +
  theme_minimal()

########################################
## 2) Likelihood-based set reliability plot (single panel)
########################################

p_lik <- ggplot(
  likelihood_rel_data_fatality,
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
    subtitle = "Conflict State Data"
  ) +
  coord_equal(xlim = c(0.5, 1), ylim = c(0.5, 1)) +
  theme_minimal()

########################################
## 3) Added-zero CP (anchored vs unanchored) – 2-panel plot
########################################

# p_zero <- ggplot(
#   cp_added_zero_rel_data_fatality,
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
#     subtitle = "Fatality Data"
#   ) +
#   coord_equal(xlim = c(0.5, 1), ylim = c(0.5, 1)) +
#   theme_minimal()

ggsave("outputs/plots/rel_curve_cp_fatal.pdf",   p_cp, width = 5, height = 5, units = "in", device = cairo_pdf)
ggsave("outputs/plots/rel_curve_lik_fatal.pdf",      p_lik, width = 5, height = 5, units = "in", device = cairo_pdf)

# Reuse the same plot style, and show undercoverage below 0.5 as well.
p_naive <- p_cp %+% naive_rel_data_fatality +
  labs(title = "Reliability curves - Naive CP") +
  coord_cartesian(xlim = c(0.5, 1), ylim = c(0, 1))
p_weighted <- p_cp %+% weighted_rel_data_fatality +
  labs(title = "Reliability curves - Weighted CP") +
  coord_cartesian(xlim = c(0.5, 1), ylim = c(0, 1))

comparison_data <- bind_rows(
  "HMM CP" = cp_original_rel_data_fatality,
  "Likelihood" = likelihood_rel_data_fatality,
  "Naive CP" = naive_rel_data_fatality,
  "Weighted CP" = weighted_rel_data_fatality,
  .id = "method"
)
comparison_data$method <- factor(comparison_data$method, levels = c("HMM CP", "Likelihood", "Naive CP", "Weighted CP"))
p_compare <- p_cp %+% comparison_data +
  facet_wrap(~ method, ncol = 2) +
  labs(title = "Reliability comparison across all four methods") +
  coord_cartesian(xlim = c(0.5, 1), ylim = c(0, 1))

ggsave("outputs/plots/rel_curve_naive_fatal.pdf", p_naive, width = 5, height = 5, units = "in", device = cairo_pdf)
ggsave("outputs/plots/rel_curve_weighted_fatal.pdf", p_weighted, width = 5, height = 5, units = "in", device = cairo_pdf)
ggsave("outputs/plots/rel_curve_compare_fatal.pdf", p_compare, width = 10, height = 8, units = "in", device = cairo_pdf)
# ggsave("outputs/plots/rel_curve_anchor_fatal.pdf",  p_zero, width = 10, height = 5, units = "in", device = cairo_pdf)
