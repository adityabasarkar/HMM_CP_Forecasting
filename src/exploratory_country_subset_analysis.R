library(arrow)
library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)

setwd("C:/Users/adity_724nfxg/Documents/git clones/HMM_CP_Forecasting")
sourceCpp("src/cpp_files/cpp_files/hmm_conf_pred_final.cpp")

my_data <- read_parquet("raw_data/cm_2024_johan.parquet")

my_countries <- c(1,2,70,149,235)

my_countries_df <- my_data[my_data$country_id %in% my_countries,,drop=F]
my_countries_df <- my_countries_df[,c("month_id","country_id","ged_target"),drop=F]
final_countries_df = NULL
for (i in my_countries) {
  subdf = my_countries_df[my_countries_df$country_id == i,,drop=F]
  subdf = subdf[order(subdf$month_id),,drop=F]
  final_countries_df = rbind(final_countries_df, subdf)
}
final_countries_df = final_countries_df[!is.na(final_countries_df$ged_target),,drop=F]

final_countries_df <- final_countries_df %>%
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
         ),
         casualties = case_when(
           ged_target > 0 ~ 1,
           ged_target <= 0 ~ 0
         )) %>%
  ungroup() %>%
  arrange(country_id, month_id)

train_len <- 50
test_len <- 2
alpha <- 0.2
accepts <- 0
k <- 4L
m <- 4L

list_of_cp_sets <- setNames(vector("list", length(my_countries)),
                            as.character(my_countries))

for (cid in my_countries) {
  current_subset <- cp_subset_data_country(cid)
  current_pair_list <- df_to_pair_list_markov(current_subset)
  
  total_len <- train_len + test_len
  maximum_start_index <- length(current_pair_list) - total_len + 1
  
  # change this later
  start <- 10
  
  sequence <- if (train_len > 0L) current_pair_list[start:(start + train_len - 1)] else list()
  test_sequence <- current_pair_list[(start + train_len):(start + total_len - 1)]
  
  cp_set <- get_cp_set(sequence, integer(0), k, m, test_len, alpha)
  condensed_set <- Filter(Negate(is.null), cp_set)
  
  # store condensed_set under the country ID
  list_of_cp_sets[[as.character(cid)]] <- condensed_set
}

cp_set_og2 = list_of_cp_sets

# find the countries with missing timepoints
cids_with_missing_timepoints = list()
unique_cids = unique(my_data$country_id)
for (cid in unique_cids) {
  subset_data = my_data[my_data$country_id == cid,]
  if (nrow(subset_data) < 412) {
    cids_with_missing_timepoints = append(cids_with_missing_timepoints, cid)
  }
}
cids_with_missing_timepoints

df_states_sample <- my_data[my_data$country_id == 102,]

# plotting stuff
cid = 1

subset_data <- df_states %>% 
  filter(country_id == cid) %>%
  select(c(country_id, month_id, state))

df_to_pair_list_markov <- function(df, state_col = "state", state_zero_index = FALSE) {
  st <- as.integer(df[[state_col]])
  if (anyNA(st)) stop("state column contains NA after integer coercion.")   # <-- added
  if (!state_zero_index) {
    st <- st - 1L
  }
  lapply(seq_len(nrow(df)), function(i) c(st[i], st[i]))  # integer length-2 vectors
}

pair_list <- df_to_pair_list_markov(subset_data)

start_idx = 230
calib_len <- 50L
N <- length(pair_list)                                                   # <-- added
if (calib_len >= N) stop("calib_len must be < number of rows.")    # <-- added
sequence      <- pair_list[1 : 50]
alpha <- 0.2
k <- 4
m <- 4
for (i in 0:3) {
  for (j in 0:3) {
    test_sequence = list(c(i,i), c(j,j))
    p_val <- get_cp_p_val(sequence, test_sequence, k, m, alpha)
    print(paste0("ij: ", i, j, " p_val: ", p_val))
  }
}


estimate_transition_probabilities(sequence, k)
estimate_observation_probabilities(sequence, k, m)
cp_set <- get_cp_p_val(sequence, test_sequence, k, m, alpha, TRUE)
p_val






