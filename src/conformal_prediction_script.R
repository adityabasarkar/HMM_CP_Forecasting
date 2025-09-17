library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("hmm_conformal_prediction.cpp")


# Functions

# testing function
test <- function(num_trials,
                 calib_len,
                 test_len,
                 alpha,
                 p,
                 b,
                 k,
                 m) {
  # ---- basic checks (mirroring test_cpp) ----
  if (num_trials <= 0) stop("num_trials must be > 0.")
  if (calib_len < 0 || test_len <= 0) stop("calib_len >= 0 and test_len > 0 required.")
  if (k != 2 || m != 2) {
    stop("This test implementation assumes k = 2 and m = 2 (binary HMM), like the C++.")
  }
  if (p < 0 || p > 1 || b < 0 || b > 1) {
    stop("p and b must be in [0,1].")
  }
  
  # ---- HMM params (binary), same as C++ ----
  pi <- c(0.5, 0.5)
  P  <- matrix(c(p, 1 - p,
                 1 - p, p), nrow = 2, byrow = TRUE)
  B  <- matrix(c(b, 1 - b,
                 1 - b, b), nrow = 2, byrow = TRUE)
  
  num_correct <- 0L
  
  for (u in seq_len(num_trials) - 1L) {
    # mimic occasional responsiveness; user can interrupt the loop anyway
    if ((u %% 128L) == 0L) utils::flush.console()
    
    sequence_length <- calib_len + test_len
    
    # simulate full sequence: list of length N of (state, obs) integer pairs
    sequence_data <- generate_hmm_sequence_list(sequence_length, pi, P, B)
    
    # calibration split (empty list if calib_len == 0, as in C++)
    calib_data <-
      if (calib_len > 0L) sequence_data[seq_len(calib_len)] else list()
    
    # test window: ground-truth states ("actual") and observations
    actual            <- integer(test_len)
    test_observations <- integer(test_len)
    for (i in seq_len(test_len)) {
      pair_ij <- sequence_data[[calib_len + i]]
      # pair_ij is c(state, obs) with 0-based indices as generated in C++
      actual[i]            <- pair_ij[[1L]]
      test_observations[i] <- pair_ij[[2L]]
    }
    
    # conformal prediction set for the test observations
    cp_set <- get_cp_set(calib_data, test_observations, k, m, alpha)
    
    # membership test: is `actual` present in cp_set?
    found <- FALSE
    for (idx in seq_along(cp_set)) {
      elt <- cp_set[[idx]]
      if (is.null(elt)) next
      if (length(elt) != test_len) next
      # exact elementwise equality
      if (all(as.integer(elt) == actual)) {
        found <- TRUE
        break
      }
    }
    if (found) num_correct <- num_correct + 1L
    
    # match the C++ printing behavior (zero-based trial index)
    cat(sprintf("TRIAL NUM: %d\n", u))
  }
  
  as.numeric(num_correct) / as.numeric(num_trials)
}

# faster testing function
test2 <- function(num_trials,
                 calib_len,
                 test_len,
                 alpha,
                 p,
                 b,
                 k,
                 m) {
  # ---- basic checks (mirroring test_cpp) ----
  if (num_trials <= 0) stop("num_trials must be > 0.")
  if (calib_len < 0 || test_len <= 0) stop("calib_len >= 0 and test_len > 0 required.")
  if (k != 2 || m != 2) {
    stop("This test implementation assumes k = 2 and m = 2 (binary HMM), like the C++.")
  }
  if (p < 0 || p > 1 || b < 0 || b > 1) {
    stop("p and b must be in [0,1].")
  }
  
  # ---- HMM params (binary), same as C++ ----
  pi <- c(0.5, 0.5)
  P  <- matrix(c(p, 1 - p,
                 1 - p, p), nrow = 2, byrow = TRUE)
  B  <- matrix(c(b, 1 - b,
                 1 - b, b), nrow = 2, byrow = TRUE)
  
  num_correct <- 0L
  
  for (u in seq_len(num_trials) - 1L) {
    # mimic occasional responsiveness; user can interrupt the loop anyway
    if ((u %% 128L) == 0L) utils::flush.console()
    
    sequence_length <- calib_len + test_len
    
    # simulate full sequence: list of length N of (state, obs) integer pairs
    sequence_data <- generate_hmm_sequence_list(sequence_length, pi, P, B)
    
    # calibration split (empty list if calib_len == 0, as in C++)
    calib_data <-
      if (calib_len > 0L) sequence_data[seq_len(calib_len)] else list()
    
    # test window: ground-truth states ("actual") and observations
    test_data <-
      if (calib_len < sequence_length) sequence_data[(calib_len + 1):sequence_length] else list()
    
    # conformal prediction set for the test observations
    quantile <- get_cp_p_val(calib_data, test_data, k, m, alpha)
    
    # membership test: is `actual` present in cp_set?
    found <- FALSE
    if (quantile > alpha) {
      found <- TRUE
    }
    if (found) num_correct <- num_correct + 1L
    
    # match the C++ printing behavior (zero-based trial index)
    cat(sprintf("TRIAL NUM: %d\n", u))
  }
  
  as.numeric(num_correct) / as.numeric(num_trials)
}

# test script
score <- test2(200, 100, 3, 0.2, 0.4, 0.3, 2, 2)
print(score)

# 0.8

##################################################
# For testing individual components
##################################################
p <- 0.3
b <- 0.75
pi <- c(0.5, 0.5)
P  <- matrix(c(p, 1-p,
               1-p, p), nrow = 2, byrow = TRUE)
B  <- matrix(c(1, 0,
               0, 1), nrow = 2, byrow = TRUE)
alpha <- 0.2
k <- 2
m <- 2

calib_len <- 100
test_len <- 3
sequence_length <- calib_len + test_len

# simulate full sequence: list of length N of (state, obs) tuples (0-based values)
seq_full <- generate_hmm_sequence_list(sequence_length, pi, P, B)
print(seq_full)
trans_probs <- estimate_transition_probabilities(seq_full, k)
obs_probs <- estimate_observation_probabilities(seq_full, k, m)
print(trans_probs)
print(obs_probs)

blocks <- ij_blocks(seq_full)
print(blocks)

permutes <- create_permutes(blocks, test_len)

train_seq = seq_full[1:100]
test_seq = seq_full[101:103]
test_seq = c(1,0,1)
start <- Sys.time()
cp_set <- get_cp_set(train_seq, test_seq, k, m, alpha)
end <- Sys.time()
time_taken <- end - start
print(time_taken)

# test permutations
nums <- c(1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)
x <- k_permutations(nums, 14)