library(dplyr)
library(tidyr)
library(ggplot2)

setwd("C:/Users/adity_724nfxg/Documents/git clones/HMM_CP_Forecasting")
data_dir = "outputs/r_objects"
rcd_data1 = readRDS(file.path(data_dir, "reliability_curves_data_no_added_zero.rds"))
rcd_data2 = readRDS(file.path(data_dir, "reliability_curves_data_added_zero.rds"))

plot_df1 <- rcd_data1 %>%
  filter(!is.na(empirical_coverage)) %>%
  mutate(
    # Make panel order + nicer labels
    sampling_type = factor(
      sampling_type,
      levels = c("random", "rolling"),
      labels = c("Random Sampling", "Rolling Window")
    ),
    test_length = factor(test_length)  # for a tidy legend
  )

plot_df2 <- rcd_data2 %>%
  filter(!is.na(empirical_coverage)) %>%
  mutate(
    # Make panel order + nicer labels
    sampling_type = factor(
      sampling_type,
      levels = c("random", "rolling"),
      labels = c("Random Sampling", "Rolling Window")
    ),
    test_length = factor(test_length)  # for a tidy legend
  )

reliability_curve_plot1 = ggplot(plot_df1, aes(x = target_coverage, y = empirical_coverage,
                    color = test_length, group = test_length)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(size = 0.7) +
  geom_point(size = 1) +
  facet_wrap(~ sampling_type, nrow = 1) +
  labs(
    title = "Reliability Curves",
    subtitle = "Calibration Length: 50    Test Lengths: 1 ... 6    No Added Zero",
    x = expression("Target coverage ("*1 - alpha*")"),
    y = "Empirical coverage",
    color = "Test length"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

reliability_curve_plot2 = ggplot(plot_df2, aes(x = target_coverage, y = empirical_coverage,
                                               color = test_length, group = test_length)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(size = 0.7) +
  geom_point(size = 1) +
  facet_wrap(~ sampling_type, nrow = 1) +
  labs(
    title = "Reliability Curves",
    subtitle = "Calibration Length: 50    Test Lengths: 1 ... 6    Added Zero",
    x = expression("Target coverage ("*1 - alpha*")"),
    y = "Empirical coverage",
    color = "Test length"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "outputs/plots/reliability_curves_no_added_zero.pdf",  # path + filename
  plot = reliability_curve_plot1,                # or specify your object: plot = my_plot
  width = 15, height = 8,                      # in inches
  units = "in",                                 # "in", "cm", or "mm"
  dpi = 300                                     # dots per inch (affects raster formats, not PDF)
)

ggsave(
  filename = "outputs/plots/reliability_curves_added_zero.pdf",  # path + filename
  plot = reliability_curve_plot2,                # or specify your object: plot = my_plot
  width = 15, height = 8,                      # in inches
  units = "in",                                 # "in", "cm", or "mm"
  dpi = 300                                     # dots per inch (affects raster formats, not PDF)
)