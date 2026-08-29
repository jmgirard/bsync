#include <Rcpp.h>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Calculate Windowed Phase Synchrony Core
// -----------------------------------------------------------------------------
//
// Inputs are *instantaneous phase* vectors (radians), extracted in R via
// gsignal::hilbert (Invariant 1: all validation and the analytic-signal step
// live in the R wrapper; the core assumes clean inputs and only does its own
// bounds checks, returning NA out of range).
//
// For each (i, tau) the phase difference dphi_j = phi_x[j] - phi_y[j + tau]
// over the window j = i .. i + w_max yields
//   plv       = |mean(exp(1i * dphi))|          (lachaux1999, p. 195)
//   rel_phase = arg(mean(exp(1i * dphi)))       (circular mean; > 0 = x ahead)
//
// Same prefix-sum shape as the M2 WCC core: per distinct tau, prefix sums of
// cos(dphi) and sin(dphi) make each window O(1) after an O(n) pass.

// [[Rcpp::export]]
List calc_wphase_cpp(NumericVector phi_x, NumericVector phi_y,
                     IntegerVector i_vals, IntegerVector tau_vals,
                     int w_max) {

  int n_x = phi_x.size();
  int n_y = phi_y.size();
  int n_calcs = i_vals.size();
  NumericVector plv(n_calcs);
  NumericVector rel_phase(n_calcs);

  // Group computation indices by distinct tau values (first-seen order).
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

  std::unordered_map<int, std::vector<std::pair<int, int>>> tau_groups;
  for (int k = 0; k < n_calcs; k++) {
    tau_groups[tau_vals[k]].emplace_back(k, i_vals[k]);
  }

  for (int tau : unique_taus) {
    const auto& group = tau_groups[tau];

    int max_i_val = 0;
    for (const auto& p : group) max_i_val = std::max(max_i_val, p.second);
    int arr_len = max_i_val + w_max; // covers 0-based indices 0..arr_len-1
    if (arr_len > n_x) arr_len = n_x;

    // Prefix sums of cos/sin of the phase difference at this tau.
    std::vector<double> p_c(arr_len + 1, 0.0);
    std::vector<double> p_s(arr_len + 1, 0.0);

    for (int j = 0; j < arr_len; j++) {
      int jy = j + tau;
      double dc = 0.0, ds = 0.0;
      if (jy >= 0 && jy < n_y) {
        double d = phi_x[j] - phi_y[jy];
        dc = std::cos(d);
        ds = std::sin(d);
      }
      p_c[j + 1] = p_c[j] + dc;
      p_s[j + 1] = p_s[j] + ds;
    }

    for (const auto& kv : group) {
      int k = kv.first;
      int i = kv.second - 1; // 0-based window start

      if (i < 0 || i + w_max >= n_x ||
          i + tau < 0 || i + tau + w_max >= n_y) {
        plv[k] = NA_REAL;
        rel_phase[k] = NA_REAL;
        continue;
      }

      int lo = i;
      int hi = i + w_max + 1;
      if (hi > arr_len) {
        plv[k] = NA_REAL;
        rel_phase[k] = NA_REAL;
        continue;
      }

      double n_w = (double)(w_max + 1);
      double c_mean = (p_c[hi] - p_c[lo]) / n_w;
      double s_mean = (p_s[hi] - p_s[lo]) / n_w;

      plv[k] = std::sqrt(c_mean * c_mean + s_mean * s_mean);
      rel_phase[k] = std::atan2(s_mean, c_mean);
    }
  }

  return List::create(_["plv"] = plv, _["rel_phase"] = rel_phase);
}
