#include <RcppDist.h>
#include <cmath>
#include <vector>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

namespace {

// Each R list element must contain one whole state value between 0 and k - 1.
int read_state(SEXP value, int k) {
  Rcpp::NumericVector x(value);
  if (x.size() != 1 || !R_finite(x[0]) || x[0] < 0 || x[0] >= k ||
      std::floor(x[0]) != x[0]) {
    Rcpp::stop("Each state must be one integer in 0,...,k-1.");
  }
  return static_cast<int>(x[0]);
}

// This structure stores the fitted log probabilities and calibration results.
struct Calibration {
  std::vector<double> log_prob;
  std::vector<double> scores;
};

// An empty mask allows every transition; otherwise, the mask must match the state count.
void validate_transitions(const arma::mat& allowable_transitions, int k) {
  if (allowable_transitions.n_elem == 0) return;
  if (k <= 0 || allowable_transitions.n_rows != static_cast<unsigned int>(k) ||
      allowable_transitions.n_cols != static_cast<unsigned int>(k) ||
      !allowable_transitions.is_finite()) {
    Rcpp::stop("allowable_transitions must be a finite k x k matrix or an empty matrix.");
  }
}

bool allowed_path(const std::vector<int>& path, int last_state,
                  const arma::mat& allowable_transitions) {
  if (allowable_transitions.n_elem == 0) return true;
  // The mask checks the first transition from history and all transitions within the path.
  // Zero entries forbid transitions, but the mask never changes the path score.
  int previous = last_state;
  for (int state : path) {
    if (allowable_transitions(previous, state) == 0.0) return false;
    previous = state;
  }
  return true;
}

// This score treats future states as independent draws from one marginal distribution.
// It sums their negative log probabilities and ignores temporal dependence and ordering.
double path_score(const std::vector<int>& path,
                  const std::vector<double>& log_prob) {
  double score = 0.0;
  for (int state : path) score -= log_prob[state];
  return score;
}

Calibration calibrate(const Rcpp::List& history, int pred_len, int k,
                      double split_fraction) {
  const int n = history.size();
  if (!R_finite(split_fraction) || split_fraction <= 0.0 ||
      split_fraction >= 1.0) {
    Rcpp::stop("split_fraction must be in (0,1).");
  }
  // The first split_fraction of the history fits the model, rounded down to whole points.
  const int fit_len = static_cast<int>(std::floor(split_fraction * n));
  // The split must leave fitting data and enough calibration points for a complete path.
  if (k <= 0 || pred_len <= 0 || fit_len <= 0 || fit_len >= n ||
      n - fit_len < pred_len) {
    Rcpp::stop("Require k and pred_len > 0, at least one fitting point, and at least pred_len calibration points.");
  }

  std::vector<int> states(n);
  for (int i = 0; i < n; ++i) {
    states[i] = read_state(history[i], k);
  }

  // Only the fitting portion contributes to the state counts.
  // Adding one count per state prevents zero probabilities and infinite path scores.
  std::vector<double> counts(k, 1.0);
  for (int i = 0; i < fit_len; ++i) {
    ++counts[states[i]];
  }
  
  Calibration result;
  result.log_prob.resize(k);
  for (int s = 0; s < k; ++s) {
    result.log_prob[s] = std::log(counts[s] / (fit_len + k));
  }

  // Every consecutive path that fits inside the calibration region is scored once.
  // The first path starts at fit_len, and the last path ends at n - 1.
  // Overlapping windows are dependent, so ordinary exchangeable coverage guarantees
  // should not automatically be claimed for this deliberately naive sequential baseline.
  const int num_windows = n - fit_len - pred_len + 1;
  result.scores.reserve(num_windows);
  for (int start = fit_len; start <= n - pred_len; ++start) {
    if (((start - fit_len) & 4095) == 0) Rcpp::checkUserInterrupt();
    std::vector<int> path(pred_len);
    for (int j = 0; j < pred_len; ++j) {
      path[j] = states[start + j];
    }
    result.scores.push_back(path_score(path, result.log_prob));
  }
  return result;
}

double p_value(double candidate_score, const std::vector<double>& scores) {
  // The candidate contributes one count to the numerator and denominator.
  // Each calibration score at least as large as the candidate adds one more count.
  // The tolerance includes numerical ties conservatively.
  double tail_count = 1.0;
  for (double score : scores) {
    if (score >= candidate_score - 1e-12) ++tail_count;
  }
  return tail_count / (scores.size() + 1.0);
}

// This function converts an observed R list into a validated candidate state path.
std::vector<int> read_path(const Rcpp::List& path, int k) {
  if (path.size() == 0) Rcpp::stop("candidate must be non-empty.");
  std::vector<int> out(path.size());
  for (int i = 0; i < path.size(); ++i) {
    out[i] = read_state(path[i], k);
  }
  return out;
}

}  // namespace

// [[Rcpp::export]]
// This function returns the p-value for one observed future path.
// The test sequence is only the candidate and never contributes to fitting or calibration.
double get_naive_p_val(Rcpp::List& sequence, Rcpp::List& test_sequence,
                       int k, double split_fraction,
                       const arma::mat& allowable_transitions) {
  std::vector<int> path = read_path(test_sequence, k);
  validate_transitions(allowable_transitions, k);
  if (sequence.size() == 0) Rcpp::stop("sequence must be non-empty.");
  const int last_state = read_state(sequence[sequence.size() - 1], k);
  // An infeasible test path is rejected before calibration or score evaluation.
  if (!allowed_path(path, last_state, allowable_transitions)) return 0.0;
  Calibration cal = calibrate(sequence, path.size(), k, split_fraction);
  return p_value(path_score(path, cal.log_prob), cal.scores);
}

// [[Rcpp::export]]
// This function returns a list of feasible future paths whose p-values exceed alpha.
Rcpp::List get_naive_cp_set(Rcpp::List& sequence, int pred_len, int k,
                            double split_fraction, double alpha,
                            const arma::mat& allowable_transitions) {
  if (!R_finite(alpha) || alpha < 0.0 || alpha > 1.0) {
    Rcpp::stop("alpha must be in [0,1].");
  }
  validate_transitions(allowable_transitions, k);
  Calibration cal = calibrate(sequence, pred_len, k, split_fraction);
  const int last_state = read_state(sequence[sequence.size() - 1], k);

  // The path limit prevents enumeration from creating an excessively large result.
  const int max_paths = 1000000;
  int num_paths = 1;
  for (int j = 0; j < pred_len; ++j) {
    if (num_paths > max_paths / k) {
      Rcpp::stop("Too many candidate paths to enumerate (limit: 1000000).");
    }
    num_paths *= k;
  }

  Rcpp::List result;
  std::vector<int> path(pred_len);
  for (int i = 0; i < num_paths; ++i) {
    if ((i & 4095) == 0) Rcpp::checkUserInterrupt();
    // Reading this index in base k produces one unique candidate state path.
    int digits = i;
    for (int j = pred_len - 1; j >= 0; --j) {
      path[j] = digits % k;
      digits /= k;
    }
    // The mask rejects infeasible paths before their scores or p-values are evaluated.
    // This check includes the transition from the last observed state into the forecast.
    if (!allowed_path(path, last_state, allowable_transitions)) continue;
    if (p_value(path_score(path, cal.log_prob), cal.scores) > alpha) {
      // Each retained path is returned as an integer vector with pred_len states.
      result.push_back(Rcpp::wrap(path));
    }
  }
  return result;
}
