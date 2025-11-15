library(arrow)
library(dplyr)
library(tidyr)
library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(future.apply)
library(ggplot2)
library(forcats)
library(gridExtra)

sourceCpp("src/cpp_files/hmm_conf_pred_final.cpp")

####################
# Initial Testing Setup
####################

# number of possible states
k = 2
# number of possible observations
m = 2
# transition matrix
P = matrix(c(0.7, 0.3, 
           0.3, 0.7), nrow=2, ncol=2)
# observation matrix
B = matrix(c(0.4, 0.6, 
             0.6, 0.4), nrow=2, ncol=2)
# number of elements in list
n = 50

generate_sequence_hmm <- function(P, B, k, m, n) {
  
  return_sequence = list()
  possible_states = 0:(k-1)
  possible_obs = 0:(m-1)
  
  first_state = sample(possible_states, 1)
  first_observation = sample(possible_obs, 1, replace=FALSE, prob=B[first_state+1,])

  return_sequence <- append(return_sequence, list(c(first_state, first_observation)))
  
  for (i in 2:n) {
    previous_state = return_sequence[[i - 1]][[1]]
    state = sample(possible_states, 1, replace=FALSE, prob=P[previous_state+1,])
    obs = sample(possible_obs, 1, replace=FALSE, prob=B[state+1,])
    return_sequence <- append(return_sequence, list(c(state, obs)))
  }
  
  return(return_sequence)
}

generate_sequence_mm <- function(P, k, n) {
  
  return_sequence = list()
  possible_states = 0:(k-1)
  
  first_state = sample(possible_states, 1)
  
  return_sequence <- append(return_sequence, list(first_state))
  
  for (i in 2:n) {
    previous_state = return_sequence[[i - 1]]
    state = sample(possible_states, 1, replace=FALSE, prob=P[previous_state+1,])
    return_sequence <- append(return_sequence, list(state))
  }
  
  return(return_sequence)
}

hmm_sequence <- generate_sequence_hmm(P, B, k, m, n)
mm_sequence <- generate_sequence_mm(P, k, n)
test_observations <- c(0, 1, 0)
pred_len = 3
num_perms = 1000
alpha = 0.2
add_zero = FALSE
mode1 = "hmm"
mode2 = "mm"

cp_set_hmm <- get_cp_set(hmm_sequence, test_observations, pred_len, num_perms, k, m, alpha, add_zero, mode1)
cp_set_mm <- get_cp_set(mm_sequence, test_observations, pred_len, num_perms, k, m, alpha, add_zero, mode2)

test_seq = list(c(0L,1L), c(1L,1L), c(1L,0L))
p_val_hmm <- get_p_val(hmm_sequence, test_seq, num_perms, k, m, add_zero, mode1)
test_seq2 = list(0, 0, 0)
p_val_mm <- get_p_val(mm_sequence, test_seq2, num_perms, k, m, add_zero, mode2)

p_val_mm <- get_p_val(mm_seq_calib, list(2, 2, 2, 2, 2, 2), 1000, 4, 4, TRUE, "mm")

ij_block <- ij_blocks(hmm_sequence)
perms <- create_permutes(ij_block, 5, 1000)
p_val <- get_p_val(mm_sequence, test_seq2, 1000, 2, 2, FALSE, "mm")
P_hat <- estimate_transition_probabilities(hmm_sequence, 2)
B_hat = matrix(c(1.0, 1.0, 
                   1.0, 1.0), nrow=2, ncol=2)
sumfilter <- filter_subroutine_sum(P_hat, B_hat, 5, hmm_sequence)
sumfilter
####################
# Test Helper Functions
####################

a = c(1,2,3)
b = c(1,2)
c = c(1,2)
d = c(1,3)
e = c(1,2,3)

same_pair(a,b) # FALSE
same_pair(b,c) # TRUE
same_pair(c,d) # FALSE
same_pair(d,e) # FALSE
same_pair(a,e) # FALSE

indices_test = c(1,20,3,4)
perms <- all_r_perms_list(indices_test, 2)
perms

indices_test = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
perms <- r_permutations(indices_test, 2, 5)

indices_test = c(1,2,3,4)
perms <- r_permutations(indices_test, 2, 9)

indices_test = c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16)
perms <- r_permutations(indices_test, 5, 200)

indices_test = c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)
perms <- r_permutations(indices_test, 5, 200)

num_lst = list(1,2,3,4,34,56,78)
list_to_string(num_lst)

x <- sample(1:5, 3, FALSE)
x
list_to_string(x)
####################
# Test Transition and Observation Matrices
####################

# Should be approximately P and B
est_trans_mat <- estimate_transition_probabilities(hmm_sequence, k)
est_obs_mat <- estimate_observation_probabilities(hmm_sequence, k, m)
est_trans_mat
est_obs_mat

####################
# Test ij_blocks
####################

ij_block_list <- ij_blocks(hmm_sequence)
length(ij_block_list)

####################
# Test create permute blocks
####################
permutables <- create_permutes(ij_block_list, 3, 50)

####################
# Test hmm filter subroutine sum
####################

filter_subroutine_sum(est_trans_mat, est_obs_mat, 3, hmm_sequence)

B = matrix(c(1.0, 0.0, 
             0.0, 1.0), nrow=2, ncol=2)
B_hat = matrix(c(1.0, 1.0, 
                 1.0, 1.0), nrow=2, ncol=2)

sequence1 = generate_sequence_hmm(P, B, k, k, 50)
est_trans_mat <- estimate_transition_probabilities(sequence1, k)
est_obs_mat <- estimate_observation_probabilities(sequence1, k, m)
filter_subroutine_sum(est_trans_mat, B_hat, 3, sequence1)

sequence2 = generate_sequence_hmm(P, B, k, k, 50)
est_trans_mat <- estimate_transition_probabilities(sequence2, k)
est_obs_mat <- estimate_observation_probabilities(sequence2, k, m)
filter_subroutine_sum(est_trans_mat, B_hat, 3, sequence2)

sequence3 = generate_sequence_hmm(P, B, k, k, 50)
est_trans_mat <- estimate_transition_probabilities(sequence3, k)
est_obs_mat <- estimate_observation_probabilities(sequence3, k, m)
filter_subroutine_sum(est_trans_mat, B_hat, 3, sequence3)


####################
# Test get_cp_set and get_p_val
####################


test1_cp_set_hmm = function(P, B, k, m, num_reps, calib_len, test_len, num_perms, alpha, add_zero) {
  
  num_included = 0
  for (i in seq_len(num_reps)) {
    
    cat(sprintf("Running repetition %d of %d ...\n", i, num_reps))
    
    sequence = generate_sequence_hmm(P, B, k, m, calib_len + test_len)
    train_seq = sequence[1 : calib_len]
    test_seq = sequence[(calib_len+1) : (calib_len+test_len)]
    
    states <- vapply(test_seq, `[`, integer(1), 1)
    obs    <- vapply(test_seq, `[`, integer(1), 2)
    
    cp_set = get_cp_set(train_seq, obs, test_len, num_perms, k, m, alpha, add_zero, "hmm")
    in_cp_set = any(vapply(cp_set, function(x) identical(x, states), logical(1)))
    
    if (in_cp_set) {
      num_included = num_included + 1
    }
    
    cat(sprintf(" → completed %d/%d; current inclusion rate: %.3f\n",
                i, num_reps, num_included / i))
    
  }
  
  final_rate <- num_included / num_reps
  cat(sprintf("Final inclusion rate: %.3f\n", final_rate))
  return(final_rate)

}

test1_cp_set_hmm(P, B, k, m, 1000, 50, 3, 1000, 0.1, FALSE)

test2_p_val_hmm = function(P, B, k, m, num_reps, calib_len, test_len, num_perms, alpha, add_zero) {
  
  num_included = 0
  for (i in seq_len(num_reps)) {
    
    cat(sprintf("Running repetition %d of %d ...\n", i, num_reps))
    
    sequence = generate_sequence_hmm(P, B, k, m, calib_len + test_len)
    train_seq = sequence[1 : calib_len]
    test_seq = sequence[(calib_len+1) : (calib_len+test_len)]
    
    p_val = get_p_val(train_seq, test_seq, num_perms, k, m, add_zero, "hmm")
    
    if (p_val > alpha) {
      num_included = num_included + 1
    }
    
    cat(sprintf(" → completed %d/%d; current inclusion rate: %.3f\n",
                i, num_reps, num_included / i))
    
  }
  
  final_rate <- num_included / num_reps
  cat(sprintf("Final inclusion rate: %.3f\n", final_rate))
  return(final_rate)
  
}

test2_p_val_hmm(P, B, k, m, 2000, 50, 3, 1000, 0.2, FALSE)


test3_cp_set_mm = function(P, k, num_reps, calib_len, test_len, num_perms, alpha, add_zero) {
  
  num_included = 0
  for (i in seq_len(num_reps)) {
    
    cat(sprintf("Running repetition %d of %d ...\n", i, num_reps))
    
    sequence = generate_sequence_mm(P, k, calib_len + test_len)
    train_seq = sequence[1 : calib_len]
    test_seq = sequence[(calib_len+1) : (calib_len+test_len)]
    
    states <- vapply(test_seq, `[`, integer(1), 1)
    
    cp_set <- get_cp_set(
      train_seq,
      c(0),
      as.integer(test_len),
      as.integer(num_perms),
      as.integer(k),
      as.integer(k),
      as.numeric(alpha),
      add_zero,
      as.character("mm")
    )
    in_cp_set = any(vapply(cp_set, function(x) identical(x, states), logical(1)))
    
    if (in_cp_set) {
      num_included = num_included + 1
    }
    
    cat(sprintf(" → completed %d/%d; current inclusion rate: %.3f\n",
                i, num_reps, num_included / i))
    
  }
  
  final_rate <- num_included / num_reps
  cat(sprintf("Final inclusion rate: %.3f\n", final_rate))
  return(final_rate)
  
}

v3 = test3_cp_set_mm(P, k, 1000, 50, 3, 1000, 0.2, FALSE)
v3

test4_p_val_mm = function(P, k, num_reps, calib_len, test_len, num_perms, alpha, add_zero) {
  
  num_included = 0
  for (i in seq_len(num_reps)) {
    
    cat(sprintf("Running repetition %d of %d ...\n", i, num_reps))
    
    sequence = generate_sequence_mm(P, k, calib_len + test_len)
    train_seq = sequence[1 : calib_len]
    test_seq = sequence[(calib_len+1) : (calib_len+test_len)]
    
    p_val <- get_p_val(
      train_seq,
      test_seq,
      as.integer(num_perms),
      as.integer(k),
      as.integer(k),
      add_zero,
      "mm"
    )
    
    if (p_val > alpha) {
      num_included = num_included + 1
    }
    
    cat(sprintf(" → completed %d/%d; current inclusion rate: %.3f\n",
                i, num_reps, num_included / i))
    
  }
  
  final_rate <- num_included / num_reps
  cat(sprintf("Final inclusion rate: %.3f\n", final_rate))
  return(final_rate)
  
}

v4 = test4_p_val_mm(P, k, 2000, 50, 3, 1000, 0.2, TRUE)
v4
