library(ggplot2)
library(dplyr)

data_dir = "outputs/r_objects"

cp_original_rel_data_fatality = readRDS(file.path(data_dir, "cp_original_rel_data_fatality.rds"))
likelihood_rel_data_fatality = readRDS(file.path(data_dir, "likelihood_rel_data_fatality.rds"))
# cp_added_zero_rel_data_fatality = readRDS(file.path(data_dir, "cp_added_zero_rel_data_fatality.rds"))

## Make sure test_length is treated as a factor for coloring
cp_original_rel_data_fatality <- cp_original_rel_data_fatality %>%
  mutate(test_length = factor(test_length))

likelihood_rel_data_fatality <- likelihood_rel_data_fatality %>%
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
# ggsave("outputs/plots/rel_curve_anchor_fatal.pdf",  p_zero, width = 10, height = 5, units = "in", device = cairo_pdf)
