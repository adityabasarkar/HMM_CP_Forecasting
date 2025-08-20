library(RcppArmadillo)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("hmm_conformal_prediction.cpp")


for (u in 1:200) {
  acc <- test_cpp(num_trials = 100,
                  calib_len  = 30,
                  test_len   = 4,
                  alpha      = 0.2,
                  p          = 0.7,
                  b          = 0.8,
                  k          = 2,
                  m          = 2)
  print(acc)
}


acc <- test_cpp(num_trials = 50,
                calib_len  = 30,
                test_len   = 4,
                alpha      = 0.2,
                p          = 0.7,
                b          = 0.8,
                k          = 2,
                m          = 2)
print(acc)

p <- 0.3
b <- 0.75
pi <- c(0.5, 0.5)
P  <- matrix(c(p, 1 - p,
               1 - p, p), nrow = 2, byrow = TRUE)
B  <- matrix(c(b, 1 - b,
               1 - b, b), nrow = 2, byrow = TRUE)
alpha <- 0.2
k <- 2
m <- 2

calib_len <- 30
test_len <- 4
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
print(permutes)
