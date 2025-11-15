library(arrow)
library(dplyr)
library(tidyr)

# Retrieve raw data from parquet
my_data <- read_parquet("data/raw/cm_2024_johan.parquet")

# Filters non-contiguous countries.
preprocessed <- my_data %>%
  filter(!is.na(month_id), !is.na(ged_target)) %>%  # combine filters
  group_by(country_id) %>%
  mutate(
    n_unique   = n_distinct(month_id),
    span       = max(month_id) - min(month_id) + 1,
    contiguous = n_unique == span
  ) %>%
  filter(contiguous) %>%
  ungroup() %>%
  select(-n_unique, -span, -contiguous)

# determine the states for each country
preprocessed <- preprocessed %>%
  group_by(country_id) %>%
  mutate(prev_fatal   = lag(ged_target, default = 0L),
         curr_fatal   = ged_target,
         prev_pos     = prev_fatal > 0L,
         curr_pos     = curr_fatal > 0L,
         state = case_when(
           !prev_pos & !curr_pos ~ 1L, # Peaceful
           !prev_pos &  curr_pos ~ 2L, # Escalation
           prev_pos &  curr_pos ~ 3L, # War
           prev_pos & !curr_pos ~ 4L  # De-escalation
         ),
         casualties = case_when(
           ged_target > 0 ~ 1,
           ged_target <= 0 ~ 0
         )) %>%
  ungroup() %>%
  arrange(country_id, month_id)

preprocessed <- preprocessed %>%
  select(month_id, country_id, ged_target, state)


out_dir <- "outputs/r_objects"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(preprocessed, file.path(out_dir, "preprocessed_data.rds"))
