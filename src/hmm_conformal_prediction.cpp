
#include <RcppDist.h>
#include <iostream>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]

Rcpp::List generate_hmm_sequence_list(int n,
                                const arma::vec& pi,
                                const arma::mat& P,
                                const arma::mat& B) {
  // prepare storage of (state,obs) pairs
  std::vector< std::pair<arma::uword,arma::uword> > seq;
  seq.reserve(n);
  
  // helper to sample a single index from a probability vector
  auto sample_discrete = [&](const arma::vec& probs) {
    double u = R::runif(0.0, 1.0), accum = 0.0;
    for (arma::uword j = 0; j < probs.n_elem; ++j) {
      accum += probs[j];
      if (u < accum) return j;
    }
    return probs.n_elem - 1;  // fallback
  };
  
  // t = 0
  arma::uword s0 = sample_discrete(pi);
  arma::uword o0 = sample_discrete(B.row(s0).t());
  seq.emplace_back(s0, o0);
  
  // t = 1 … n‑1
  for (int t = 1; t < n; ++t) {
    arma::uword prev = seq[t-1].first;
    arma::uword s = sample_discrete(P.row(prev).t());
    arma::uword o = sample_discrete(B.row(s).t());
    seq.emplace_back(s, o);
  }
  
  // convert to an R List of integer vectors length=2
  Rcpp::List out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = Rcpp::IntegerVector::create(
      static_cast<int>(seq[i].first),
      static_cast<int>(seq[i].second)
    );
  }
  
  return out;
}

// [[Rcpp::export]]
arma::mat estimate_transition_probabilities(Rcpp::List& sequence, int k) {
  int n = sequence.size();
  // initialize count matrix
  arma::mat counts = arma::zeros<arma::mat>(k, k);
  
  // 1) Count transitions
  for (int t = 0; t < n - 1; ++t) {
    Rcpp::IntegerVector pair_t   = sequence[t];
    Rcpp::IntegerVector pair_tp1 = sequence[t + 1];
    int i = pair_t[0];
    int j = pair_tp1[0];
    counts(i, j) += 1.0;
  }
  
  // 2) Normalize each row by its total outgoing count
  arma::vec row_sums = arma::sum(counts, 1);  // sum across columns
  
  arma::mat P = arma::zeros<arma::mat>(k, k);
  for (int i = 0; i < k; ++i) {
    if (row_sums[i] > 0) {
      P.row(i) = counts.row(i) / row_sums[i];
    }
    // if row_sums[i] == 0, leave P(i,*) = 0
  }
  
  return P;
}

//[[Rcpp::export]]
arma::mat estimate_observation_probabilities(Rcpp::List& sequence, int k, int m) {
  
  // get the size of the sequence
  int n = sequence.size();
  // initialize count matrix
  arma::mat counts = arma::zeros<arma::mat>(k, m);
  
  // count transitions
  for (int i = 0; i < n; ++i) {
    Rcpp::IntegerVector pair_i = sequence[i];
    int a = pair_i[0];
    int b = pair_i[1];
    counts(a, b) += 1;
  }
  
  arma::vec row_sums = arma::sum(counts, 1);
  arma::mat B = arma::zeros<arma::mat>(k, m);
  for (int i = 0; i < k; ++i) {
    if (row_sums[i] > 0) {
      B.row(i) = counts.row(i) / row_sums[i];
    }
    // if a state never appears, its row remains zero
  }
  
  return B;
}

// helper: compare two (i, j) tuples stored as R vectors
static inline bool equal_tuple_ij(SEXP a, SEXP b) {
  Rcpp::IntegerVector va(a), vb(b);           // will coerce if numeric
  return va.size() == 2 && vb.size() == 2 && va[0] == vb[0] && va[1] == vb[1];
}

//[[Rcpp::export]]
Rcpp::List ij_blocks(const Rcpp::List& sequence) {
  int n = sequence.size();
  Rcpp::List blocks;                    // outer list: blocks
  if (n == 0) return blocks;            // empty input => empty output

  SEXP ij_tuple = sequence[n - 1];      // last tuple
  Rcpp::List current_block;             // inner list: accumulating tuples

  for (int i = 0; i < n; ++i) {
    if (i == 0) {
      current_block.push_back(sequence[i]);
    } else {
      if (equal_tuple_ij(sequence[i], ij_tuple)) {
        // start a new block at occurrences of the last tuple
        blocks.push_back(current_block);
        current_block = Rcpp::List::create(sequence[i]);
      } else {
        current_block.push_back(sequence[i]);
      }
    }
  }

  // append the final block
  blocks.push_back(current_block);
  return blocks;
}

std::vector<std::vector<int>> k_permutations(std::vector<int> nums, int k) {
  std::vector<std::vector<int>> result;
  if (k < 0 || k > (int)nums.size()) return result;
  if (k == 0) { result.push_back({}); return result; }
  
  // Pair values with original indices, then sort by value to handle duplicates
  std::vector<std::pair<int,int>> a;
  a.reserve(nums.size());
  for (int i = 0; i < (int)nums.size(); ++i) a.emplace_back(nums[i], i);
  std::sort(a.begin(), a.end(), [](auto& x, auto& y){ return x.first < y.first; });
  
  std::vector<int> cur; cur.reserve(k);
  std::vector<char> used(a.size(), 0);
  
  std::function<void()> dfs = [&]() {
    if ((int)cur.size() == k) { result.push_back(cur); return; }
    for (int i = 0; i < (int)a.size(); ++i) {
      if (used[i]) continue;
      // Skip duplicates: if same value as previous and previous not used in this level
      if (i > 0 && a[i].first == a[i-1].first && !used[i-1]) continue;
      
      used[i] = 1;
      cur.push_back(a[i].first);
      dfs();
      cur.pop_back();
      used[i] = 0;
    }
  };
  
  dfs();
  return result;
}

// [[Rcpp::export]]
Rcpp::List create_permutes(const Rcpp::List& blocks, int test_len) {
  const int nb = blocks.size();
  Rcpp::List complete_lists;                 // output
  if (nb == 0) return complete_lists;        // [] -> []
  
  // last block and its first tuple (cut marker: ij_tuple = blocks[-1][0])
  Rcpp::List last_block = blocks[nb - 1];
  if (last_block.size() == 0)
    Rcpp::stop("Last block must contain at least one tuple.");
  SEXP ij_tuple = last_block[0];
  
  std::vector<int> permutable_indices;
  for (int i = 0; i < nb; i++) {
    
    if (i == 0) {
      Rcpp::List first_block = blocks[i];
      SEXP first_tuple = first_block[0];
      if (ij_tuple == first_tuple) {
        permutable_indices.push_back(i);
      }
    }
    if (i != 0 && i != nb - 1) {
      permutable_indices.push_back(i);
    }
  }
  
  std::vector<std::vector<int>> perms = k_permutations(permutable_indices, test_len);
  Rcpp::List return_lst;
  for (std::vector<int> perm: perms) {
    Rcpp::List current_lst;
    for (int ind: perm) {
      Rcpp::List current_block = blocks[ind];
      for (int i = 0; i < current_block.size(); i++) {
        current_lst.push_back(current_block[i]);
      }
    }
    current_lst.push_back(ij_tuple);
    return_lst.push_back(current_lst);
  }
  
  return return_lst;
}



// [[Rcpp::export]]
double hmm_filter_subroutine_sum(
  const arma::mat& P,
  const arma::mat& B,
  int test_len,
  const Rcpp::List& sequence) 
{
  
  const int n = sequence.size();
  const int k = P.n_rows;
  
  // checks
  if (k <= 0) {
    Rcpp::stop("P must be non-empty.");
  }
  if ((int)P.n_cols != k) {
    Rcpp::stop("P must be square k x k.");
  }
  if ((int)B.n_rows != k) {
    Rcpp::stop("B must have k rows.");
  }
  if (test_len < 0) {
    Rcpp::stop("test_len must be >= 0.");
  }
  if (n == 0) {
    return 0.0;
  }
  if (test_len > n - 1) {
    Rcpp::stop("test_len (%d) is too large for sequence length (%d).", test_len, n);
  }
  const int start_idx = n - test_len - 1;  // index of sequence[-test_len-1]
  if (start_idx < 0) {
    Rcpp::stop("Invalid start index computed.");
  }
  
  // ---- first_vec: one-hot at the state of sequence[start_idx] ----
  Rcpp::IntegerVector first_pair = sequence[start_idx];
  if (first_pair.size() != 2) {
    Rcpp::stop("sequence[%d] is not a 2-element tuple.", start_idx + 1);
  }
  int s0 = first_pair[0];
  if (s0 < 0 || s0 >= k) {
    Rcpp::stop("State index out of bounds at start: %d (k=%d).", s0, k);
  }
  
  std::vector<arma::vec> prob_vecs;
  prob_vecs.reserve(test_len + 1);
  
  arma::vec first_vec = arma::zeros<arma::vec>(k);
  first_vec(s0) = 1.0;
  prob_vecs.push_back(first_vec);
  
  // ---- iterate i = start_idx .. n-2 (length == test_len) ----
  const arma::mat Pt = P.t();
  const double eps = 1e-12;
  
  for (int i = start_idx; i <= n - 2; ++i) {
    Rcpp::IntegerVector pair_i = sequence[i];
    if (pair_i.size() != 2) {
      Rcpp::stop("sequence[%d] is not a 2-element tuple.", i + 1);
    }
    int obs = pair_i[1];
    if (obs < 0 || obs >= (int)B.n_cols) {
      Rcpp::stop("Observation index out of bounds at %d: %d (m=%d).",
                 i + 1, obs, (int)B.n_cols);
    }
    
    // temp = diag(B[, obs]) * (P^T * prev)  --> elementwise multiply
    const arma::vec& prev = prob_vecs.back();
    arma::vec temp = B.col(obs) % (Pt * prev);
    
    // stabilize and normalize
    temp += eps;
    double s = arma::accu(temp);
    if (s <= 0.0) {
      Rcpp::stop("Normalization failed (sum<=0) at step %d.", i - start_idx + 1);
    }
    temp /= s;
    
    prob_vecs.push_back(std::move(temp));
  }
  
  // ---- sum over i=0..test_len-1 of prob_vecs[i+1][ true_state_at(start_idx+i) ] ----
  double return_sum = 0.0;
  for (int i = 0; i < test_len; ++i) {
    Rcpp::IntegerVector pair_i = sequence[start_idx + i];
    int st = pair_i[0];
    if (st < 0 || st >= k) {
      Rcpp::stop("State index out of bounds at %d: %d (k=%d).", start_idx + i + 1, st, k);
    }
    return_sum += prob_vecs[i + 1](st);
  }
  
  return return_sum;
}

// [[Rcpp::export]]
Rcpp::List get_cp_set(const Rcpp::List& sequence,
                      const Rcpp::IntegerVector& test_observations,
                      int k, int m, double alpha) {
  const int n = sequence.size();
  const int L = test_observations.size();
  if (L <= 0) Rcpp::stop("test_observations must have positive length.");
  if (k <= 0 || m <= 0) Rcpp::stop("k and m must be positive.");
  if (n == 0) Rcpp::stop("sequence must be non-empty.");
  
  Rcpp::List cp_set;                  // output: list of length k^L (potentially huge)
  
  // Candidate state vector (odometer over base-k with L digits)
  std::vector<int> state(L, 0);
  bool done = false;
  
  while (!done) {
    // ---- Build augmented sequence: sequence + zip(state, test_observations) ----
    Rcpp::List augmented(n + L);
    for (int i = 0; i < n; ++i) augmented[i] = sequence[i];
    for (int i = 0; i < L; ++i) {
      int s = state[i];
      int o = test_observations[i];
      if (s < 0 || s >= k) Rcpp::stop("State out of range in candidate.");
      if (o < 0 || o >= m) Rcpp::stop("Observation out of range in test_observations.");
      augmented[n + i] = Rcpp::IntegerVector::create(s, o);
    }
    
    // ---- Estimate P and B from augmented sequence ----
    arma::mat P_hat  = estimate_transition_probabilities(augmented, k);
    arma::mat B_hat  = estimate_observation_probabilities(augmented, k, m);
    
    // ---- Build permuted lists and compute conformity scores ----
    Rcpp::List blocks        = ij_blocks(augmented);
    Rcpp::List permuted_lists = create_permutes(blocks, L);
    
    const int M = permuted_lists.size();
    if (M <= 0) {
      // No permutations produced (shouldn't happen with your create_permutes).
      cp_set.push_back(R_NilValue);
    } else {
      // score(lst) = 1 - avg_prob; avg_prob uses hmm_filter_subroutine_sum(P,B,L,lst)
      auto score_of = [&](const Rcpp::List& lst) -> double {
        double s = hmm_filter_subroutine_sum(P_hat, B_hat, L, lst); // sum over L steps
        return 1.0 - (s / static_cast<double>(L));
      };
      
      // Baseline is the first score (conformity_scores[0])
      double baseline = score_of(Rcpp::List(permuted_lists[0]));
      int ge_count = 0;
      for (int i = 0; i < M; ++i) {
        double sc = score_of(Rcpp::List(permuted_lists[i]));
        if (sc >= baseline) ++ge_count;
      }
      double q = static_cast<double>(ge_count) / static_cast<double>(M);
      
      if (q > alpha) {
        // Keep this candidate: return as IntegerVector
        Rcpp::IntegerVector keep(L);
        for (int i = 0; i < L; ++i) keep[i] = state[i];
        cp_set.push_back(keep);
      } else {
        // Discard: push NULL to mirror Python's None
        cp_set.push_back(R_NilValue);
      }
    }
    
    // ---- Increment odometer in base-k over L positions ----
    int pos = L - 1;
    while (pos >= 0) {
      state[pos] += 1;
      if (state[pos] < k) break;  // no carry needed
      state[pos] = 0;             // carry
      --pos;
    }
    if (pos < 0) done = true;     // overflow -> finished all k^L candidates
  }
  
  return cp_set;
}


// [[Rcpp::export]]
double test_cpp(int num_trials,
                int calib_len,
                int test_len,
                double alpha,
                double p,
                double b,
                int k,
                int m) {
  // ---- basic checks ----
  if (num_trials <= 0) Rcpp::stop("num_trials must be > 0.");
  if (calib_len < 0 || test_len <= 0) Rcpp::stop("calib_len >= 0 and test_len > 0 required.");
  if (k != 2 || m != 2) {
    Rcpp::stop("This test_cpp implementation assumes k = 2 and m = 2 (binary HMM), like the Python.");
  }
  if (p < 0.0 || p > 1.0 || b < 0.0 || b > 1.0) {
    Rcpp::stop("p and b must be in [0,1].");
  }
  
  // ---- HMM params (binary) ----
  arma::vec pi(2);         pi[0] = 0.5; pi[1] = 0.5;
  arma::mat P(2,2);        P(0,0)=p;   P(0,1)=1-p;  P(1,0)=1-p;  P(1,1)=p;
  arma::mat B(2,2);        B(0,0)=b;   B(0,1)=1-b;  B(1,0)=1-b;  B(1,1)=b;
  
  int num_correct = 0;
  
  for (int u = 0; u < num_trials; ++u) {
    if ((u & 127) == 0) Rcpp::checkUserInterrupt();   // stay responsive
    
    const int sequence_length = calib_len + test_len;
    
    // simulate full sequence: list of length N of (state, obs) tuples
    Rcpp::List sequence_data = generate_hmm_sequence_list(sequence_length, pi, P, B);
    
    // calibration split
    Rcpp::List calib_data(calib_len);
    for (int i = 0; i < calib_len; ++i) calib_data[i] = sequence_data[i];
    
    // test window: ground-truth states ("actual") and observations
    Rcpp::IntegerVector actual(test_len);
    Rcpp::IntegerVector test_observations(test_len);
    for (int i = 0; i < test_len; ++i) {
      Rcpp::IntegerVector pair_ij = sequence_data[calib_len + i];
      actual[i]            = pair_ij[0];   // state
      test_observations[i] = pair_ij[1];   // observation
    }
    
    // conformal prediction set for the test observations
    Rcpp::List cp_set = get_cp_set(calib_data, test_observations, k, m, alpha);
    
    // membership test: is `actual` present in cp_set?
    bool found = false;
    for (int idx = 0; idx < cp_set.size() && !found; ++idx) {
      SEXP elt = cp_set[idx];
      if (Rf_isNull(elt)) continue;
      Rcpp::IntegerVector cand(elt);
      if (cand.size() != test_len) continue;
      bool equal = true;
      for (int j = 0; j < test_len; ++j) {
        if (cand[j] != actual[j]) { equal = false; break; }
      }
      if (equal) found = true;
    }
    if (found) ++num_correct;
    Rprintf("TRIAL NUM: %d\n", u);
  }
  
  return static_cast<double>(num_correct) / static_cast<double>(num_trials);
}