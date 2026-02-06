library(dplyr)
library(tidyr)

# Retrieve raw data from parquet
my_data <- read.csv("data/raw/markov_data2507.csv")
country_names <- read.csv("data/raw/country_list.csv")

# Filters non-contiguous countries.
preprocessed <- my_data %>%
  filter(!is.na(month_id), !is.na(target_ged)) %>%  # combine filters
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
  mutate(prev_fatal   = lag(target_ged, default = 0L),
         curr_fatal   = target_ged,
         prev_pos     = prev_fatal > 0L,
         curr_pos     = curr_fatal > 0L,
         state = case_when(
           !prev_pos & !curr_pos ~ 1L, # Peaceful
           !prev_pos &  curr_pos ~ 2L, # Escalation
           prev_pos &  curr_pos ~ 3L, # War
           prev_pos & !curr_pos ~ 4L  # De-escalation
         ),
         casualties = case_when(
           target_ged > 0 ~ 1,
           target_ged <= 0 ~ 0
         )) %>%
  ungroup() %>%
  arrange(country_id, month_id)

preprocessed <- preprocessed %>%
  left_join(country_names, by = "country_id") %>%
  select(month_id, country_id, name, target_ged, state, year, month)

out_dir <- "outputs/r_objects"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(preprocessed, file.path(out_dir, "preprocessed_data_full.rds"))

preprocessed <- preprocessed %>%
  group_by(country_id) %>%
  filter(
    mean(state == 1) <= 0.99,
    sum(state != 1) >= 5,
    !if_any(everything(), is.na)
  ) %>%
  ungroup()


out_dir <- "outputs/r_objects"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(preprocessed, file.path(out_dir, "preprocessed_data.rds"))
