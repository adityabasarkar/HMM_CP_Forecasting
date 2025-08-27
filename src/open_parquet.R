library(arrow)
library(dplyr)
library(tidyr)


my_data <- read_parquet("cm_2024_johan.parquet")

plot(my_data[my_data[,'country_id']==70, 'ged_target'])

print(my_data["month_id"])

# month_id
# country_id
# ged_target
# states
# Format data so that states exists. 


# df: data.frame with integer columns: month_id, country_id, ged_target

# 1) Make a full country x month panel; treat missing months as 0 fatalities
df_full <- my_data %>%
  complete(month_id = unique(month_id), country_id,
           fill = list(ged_target = 0L))

# 2) Compute previous-month fatalities by country
df_states <- df_full %>%
  group_by(country_id) %>%
  mutate(prev_fatal   = dplyr::lag(ged_target, default = 0L),
         curr_fatal   = ged_target,
         prev_pos     = prev_fatal > 0L,
         curr_pos     = curr_fatal > 0L,
         state = case_when(
           !prev_pos & !curr_pos ~ 1L, # Peaceful
           !prev_pos &  curr_pos ~ 2L, # Escalation
           prev_pos &  curr_pos ~ 3L, # War
           prev_pos & !curr_pos ~ 4L  # De-escalation
         )) %>%
  ungroup() %>%
  arrange(country_id, month_id)