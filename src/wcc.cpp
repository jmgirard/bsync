#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Calculate Windowed Cross-Correlation Core
// -----------------------------------------------------------------------------
//
// M2 rewrite: NA-aware prefix-sum algorithm.
//
// For each distinct lag tau, we build six prefix-sum arrays over the valid
// (non-NA) positions. This reduces per-(i, tau) work from O(w_max) to O(1)
// after an O(n) preprocessing pass per tau, giving O(n_c * n) total cost vs.
// the prior O(n_c * n_r * w_max).
//
// NA modes (na_rm):
//   true  — pairwise deletion: a position j is valid iff both x[j] and
//            y[j+tau] are non-NA. The count prefix p_n[w] then tells
//            calc_wcc_cpp how many valid pairs are in any window.
//   false — any-NA: a window is NA iff its NA-count prefix is > 0.
//
// Signature is unchanged from the pre-M2 core so RcppExports require no diff.

// [[Rcpp::export]]
NumericVector calc_wcc_cpp(NumericVector x, NumericVector y,
                           IntegerVector i_vals, IntegerVector tau_vals,
                           int w_max, bool na_rm = true) {

  int n_x = x.size();
  int n_y = y.size();
  int n_calcs = i_vals.size();
  NumericVector results(n_calcs);

  // --- Group computation indices by distinct tau values ---------------------

  // Collect unique taus preserving first-seen order
  std::vector<int> unique_taus;
  {
    std::unordered_map<int, int> tau_seen;
    for (int k = 0; k < n_calcs; k++) {
      int tau = tau_vals[k];
      if (tau_seen.find(tau) == tau_seen.end()) {
        tau_seen[tau] = (int)unique_taus.size();
        unique_taus.push_back(tau);
      }
    }
  }

  // Map tau -> list of (result-index k, i_val)
  std::unordered_map<int, std::vector<std::pair<int,int>>> tau_groups;
  for (int k = 0; k < n_calcs; k++) {
    tau_groups[tau_vals[k]].emplace_back(k, i_vals[k]);
  }

  // --- Process each tau block ----------------------------------------------

  for (int tau : unique_taus) {
    const auto& group = tau_groups[tau];

    // Determine the range of x-start positions we'll need.
    // i_vals are 1-based; C++ windows span x[i-1 .. i-1+w_max].
    // We need prefix arrays covering at least 0..max(i_val)-1+w_max.
    int max_i_val = 0;
    for (const auto& p : group) max_i_val = std::max(max_i_val, p.second);
    int arr_len = max_i_val + w_max; // covers indices 0..arr_len-1 (0-based)

    // Clamp to actual array lengths
    if (arr_len > n_x) arr_len = n_x;

    // Build per-position arrays (1-based prefix sums over 0..arr_len-1)
    // p_n   : count of valid pairs
    // p_x   : sum of x at valid positions
    // p_y   : sum of y[j+tau] at valid positions
    // p_x2  : sum of x^2 at valid positions
    // p_y2  : sum of y[j+tau]^2 at valid positions
    // p_xy  : sum of x * y[j+tau] at valid positions
    // p_na  : count of NA pairs (for na_rm=false mode)

    // prefix arrays have size arr_len+1 (index 0 is zero sentinel)
    std::vector<int>    p_n (arr_len + 1, 0);
    std::vector<double> p_x (arr_len + 1, 0.0);
    std::vector<double> p_y (arr_len + 1, 0.0);
    std::vector<double> p_x2(arr_len + 1, 0.0);
    std::vector<double> p_y2(arr_len + 1, 0.0);
    std::vector<double> p_xy(arr_len + 1, 0.0);
    std::vector<int>    p_na(arr_len + 1, 0);

    for (int j = 0; j < arr_len; j++) {
      int jy = j + tau; // y index (0-based)

      bool x_na = (j  < 0 || j  >= n_x) || NumericVector::is_na(x[j]);
      bool y_na = (jy < 0 || jy >= n_y) || NumericVector::is_na(y[jy]);

      int   dn  = 0;
      double dx = 0.0, dy = 0.0, dx2 = 0.0, dy2 = 0.0, dxy = 0.0;
      int   dna = 0;

      if (!x_na && !y_na) {
        double vx = x[j], vy = y[jy];
        dn  = 1;
        dx  = vx;
        dy  = vy;
        dx2 = vx * vx;
        dy2 = vy * vy;
        dxy = vx * vy;
      } else {
        dna = 1;
      }

      p_n [j + 1] = p_n [j] + dn;
      p_x [j + 1] = p_x [j] + dx;
      p_y [j + 1] = p_y [j] + dy;
      p_x2[j + 1] = p_x2[j] + dx2;
      p_y2[j + 1] = p_y2[j] + dy2;
      p_xy[j + 1] = p_xy[j] + dxy;
      p_na[j + 1] = p_na[j] + dna;
    }

    // --- Evaluate each window in this tau group ---------------------------

    for (const auto& kv : group) {
      int k     = kv.first;
      int i_val = kv.second;
      int i     = i_val - 1; // 0-based x start

      // Bounds check: window must be valid in both arrays
      if (i < 0 || i + w_max >= n_x ||
          i + tau < 0 || i + tau + w_max >= n_y) {
        results[k] = NA_REAL;
        continue;
      }

      // Prefix window is [i .. i+w_max] inclusive → prefix[i+w_max+1] - prefix[i]
      int lo = i;
      int hi = i + w_max + 1; // exclusive upper bound in prefix array

      if (hi > arr_len) {
        results[k] = NA_REAL;
        continue;
      }

      // na_rm=false: any NA pair in this window → NA result
      if (!na_rm) {
        int na_count = p_na[hi] - p_na[lo];
        if (na_count > 0) {
          results[k] = NA_REAL;
          continue;
        }
      }

      int    valid_n = p_n [hi] - p_n [lo];
      double sum_x   = p_x [hi] - p_x [lo];
      double sum_y   = p_y [hi] - p_y [lo];
      double sum_x2  = p_x2[hi] - p_x2[lo];
      double sum_y2  = p_y2[hi] - p_y2[lo];
      double sum_xy  = p_xy[hi] - p_xy[lo];

      if (valid_n <= 1) {
        results[k] = NA_REAL;
        continue;
      }

      double var_x = (sum_x2 - (sum_x * sum_x) / valid_n) / (valid_n - 1);
      double var_y = (sum_y2 - (sum_y * sum_y) / valid_n) / (valid_n - 1);
      double cov_xy = (sum_xy - (sum_x * sum_y) / valid_n) / (valid_n - 1);

      if (var_x <= 0.0 || var_y <= 0.0) {
        results[k] = NA_REAL;
      } else {
        results[k] = cov_xy / sqrt(var_x * var_y);
      }
    }
  }

  return results;
}
