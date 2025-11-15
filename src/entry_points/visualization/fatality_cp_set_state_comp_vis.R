library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(scales)

setwd("C:/Users/adity_724nfxg/Documents/git clones/HMM_CP_Forecasting")
data_dir = "outputs/r_objects"
cp_set_results1 = readRDS(file.path(data_dir, "cp_set_results_no_added_zero_nonextrapolation.rds"))
cp_set_results2 = readRDS(file.path(data_dir, "cp_set_results_added_zero_nonextrapolation.rds"))
cp_set_results3 = readRDS(file.path(data_dir, "cp_set_results_no_added_zero_extrapolation.rds"))
cp_set_results4 = readRDS(file.path(data_dir, "cp_set_results_added_zero_extrapolation.rds"))

T1 = 6

# ---- reshape cp_set_results into a long tidy frame ----
# cp_set_results is a named list: "country 1" -> 4xT1 matrix of proportions
plot_df1 <- lapply(names(cp_set_results1), function(nm) {
  mat <- cp_set_results1[[nm]]              # rows: state 1..4, cols: timepoint 1..T1
  as.data.frame(mat) %>%
    mutate(state_id = rownames(mat)) %>%
    pivot_longer(starts_with("timepoint"),
                 names_to = "timepoint",
                 values_to = "prop") %>%
    mutate(
      country = str_replace(nm, "^country\\s+", "")
    )
}) %>% bind_rows()

plot_df2 <- lapply(names(cp_set_results2), function(nm) {
  mat <- cp_set_results2[[nm]]              # rows: state 1..4, cols: timepoint 1..T1
  as.data.frame(mat) %>%
    mutate(state_id = rownames(mat)) %>%
    pivot_longer(starts_with("timepoint"),
                 names_to = "timepoint",
                 values_to = "prop") %>%
    mutate(
      country = str_replace(nm, "^country\\s+", "")
    )
}) %>% bind_rows()

plot_df3 <- lapply(names(cp_set_results3), function(nm) {
  mat <- cp_set_results3[[nm]]              # rows: state 1..4, cols: timepoint 1..T1
  as.data.frame(mat) %>%
    mutate(state_id = rownames(mat)) %>%
    pivot_longer(starts_with("timepoint"),
                 names_to = "timepoint",
                 values_to = "prop") %>%
    mutate(
      country = str_replace(nm, "^country\\s+", "")
    )
}) %>% bind_rows()

plot_df4 <- lapply(names(cp_set_results4), function(nm) {
  mat <- cp_set_results4[[nm]]              # rows: state 1..4, cols: timepoint 1..T1
  as.data.frame(mat) %>%
    mutate(state_id = rownames(mat)) %>%
    pivot_longer(starts_with("timepoint"),
                 names_to = "timepoint",
                 values_to = "prop") %>%
    mutate(
      country = str_replace(nm, "^country\\s+", "")
    )
}) %>% bind_rows()

# Optional: nicer labels
state_labels <- c(
  "state 1" = "Peaceful",
  "state 2" = "Escalation",
  "state 3" = "War",
  "state 4" = "De-escalation"
)

# order states in legend/stack (top-to-bottom is reverse of this order)
plot_df1 <- plot_df1 %>%
  mutate(
    timepoint = factor(timepoint, levels = paste0("timepoint ", 1:T1),
                       labels = paste("Month", 1:T1)),
    state = factor(state_id, levels = names(state_labels), labels = unname(state_labels)),
    country = factor(country, levels = unique(country))
  )
  
plot_df2 <- plot_df2 %>%
  mutate(
    timepoint = factor(timepoint, levels = paste0("timepoint ", 1:T1),
                       labels = paste("Month", 1:T1)),
    state = factor(state_id, levels = names(state_labels), labels = unname(state_labels)),
    country = factor(country, levels = unique(country))
  )

plot_df3 <- plot_df3 %>%
  mutate(
    timepoint = factor(timepoint, levels = paste0("timepoint ", 1:T1),
                       labels = paste("Month", 1:T1)),
    state = factor(state_id, levels = names(state_labels), labels = unname(state_labels)),
    country = factor(country, levels = unique(country))
  )

plot_df4 <- plot_df4 %>%
  mutate(
    timepoint = factor(timepoint, levels = paste0("timepoint ", 1:T1),
                       labels = paste("Month", 1:T1)),
    state = factor(state_id, levels = names(state_labels), labels = unname(state_labels)),
    country = factor(country, levels = unique(country))
  )

# ---- colors to match your example (feel free to tweak) ----
fill_cols <- c(
  "Peaceful"       = "#9ecae1",  # light blue
  "Escalation"     = "#f28e2b",  # orange
  "War"            = "#e15759",  # red
  "De-escalation"  = "#59a14f"   # green
)

# ---- plot ----
p1 <- ggplot(plot_df1, aes(x = timepoint, y = prop, fill = state)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85, color = NA) +
  facet_wrap(~ country, ncol = 1, labeller = labeller(country = function(x) paste("Country", x))) +
  scale_fill_manual(values = fill_cols, name = "State") +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.25), expand = expansion(c(0,0.02))) +
  labs(
    title = "State Composition Across Prediction Horizon",
    subtitle = "Conformal Prediction Set Results (T = Until month_id 508, T1 = 6, alpha = 0.2, No Added Zero, No Extrapolation)",
    x = NULL,
    y = "Percentage of candidates"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5)),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )
  
p2 <- ggplot(plot_df2, aes(x = timepoint, y = prop, fill = state)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85, color = NA) +
  facet_wrap(~ country, ncol = 1, labeller = labeller(country = function(x) paste("Country", x))) +
  scale_fill_manual(values = fill_cols, name = "State") +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.25), expand = expansion(c(0,0.02))) +
  labs(
    title = "State Composition Across Prediction Horizon",
    subtitle = "Conformal Prediction Set Results (T = Until month_id 508, T1 = 6, alpha = 0.2, Added Zero, No Extrapolation)",
    x = NULL,
    y = "Percentage of candidates"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5)),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

p3 <- ggplot(plot_df3, aes(x = timepoint, y = prop, fill = state)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85, color = NA) +
  facet_wrap(~ country, ncol = 1, labeller = labeller(country = function(x) paste("Country", x))) +
  scale_fill_manual(values = fill_cols, name = "State") +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.25), expand = expansion(c(0,0.02))) +
  labs(
    title = "State Composition Across Prediction Horizon",
    subtitle = "Conformal Prediction Set Results (T = ALL, T1 = 6, alpha = 0.2, No Added Zero, Extrapolation)",
    x = NULL,
    y = "Percentage of candidates"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5)),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

p4 <- ggplot(plot_df4, aes(x = timepoint, y = prop, fill = state)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85, color = NA) +
  facet_wrap(~ country, ncol = 1, labeller = labeller(country = function(x) paste("Country", x))) +
  scale_fill_manual(values = fill_cols, name = "State") +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.25), expand = expansion(c(0,0.02))) +
  labs(
    title = "State Composition Across Prediction Horizon",
    subtitle = "Conformal Prediction Set Results (T = ALL, T1 = 6, alpha = 0.2, Added Zero, Extrapolation)",
    x = NULL,
    y = "Percentage of candidates"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(margin = margin(t = 5)),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

ggsave(
  filename = "outputs/plots/cp_set_results_noext_no_added_zero.pdf",  # path + filename
  plot = p1,                # or specify your object: plot = my_plot
  width = 15, height = 16,                      # in inches
  units = "in",                                 # "in", "cm", or "mm"
  dpi = 300                                     # dots per inch (affects raster formats, not PDF)
)

ggsave(
  filename = "outputs/plots/cp_set_results_noext_added_zero.pdf",  # path + filename
  plot = p2,                # or specify your object: plot = my_plot
  width = 15, height = 16,                      # in inches
  units = "in",                                 # "in", "cm", or "mm"
  dpi = 300                                     # dots per inch (affects raster formats, not PDF)
)

ggsave(
  filename = "outputs/plots/cp_set_results_ext_no_added_zero.pdf",  # path + filename
  plot = p3,                # or specify your object: plot = my_plot
  width = 15, height = 16,                      # in inches
  units = "in",                                 # "in", "cm", or "mm"
  dpi = 300                                     # dots per inch (affects raster formats, not PDF)
)

ggsave(
  filename = "outputs/plots/cp_set_results_ext_added_zero.pdf",  # path + filename
  plot = p4,                # or specify your object: plot = my_plot
  width = 15, height = 16,                      # in inches
  units = "in",                                 # "in", "cm", or "mm"
  dpi = 300                                     # dots per inch (affects raster formats, not PDF)
)