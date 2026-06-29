#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Calculate Windowed Cross-Correlation Core
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector calc_wcc_cpp(NumericVector x, NumericVector y,
                           IntegerVector i_vals, IntegerVector tau_vals,
                           int w_max, bool na_rm = true) {

  int n_calcs = i_vals.size();
  NumericVector results(n_calcs);

  for(int k = 0; k < n_calcs; k++) {
    int i = i_vals[k] - 1; // 0-based indexing
    int tau = tau_vals[k];

    // Safety check to ensure windows do not go out of bounds
    if (i < 0 || i + w_max >= x.size() || i + tau < 0 || i + tau + w_max >= y.size()) {
      results[k] = NA_REAL;
      continue;
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0;
    double sum_x2 = 0.0, sum_y2 = 0.0;
    int valid_n = 0;
    bool has_na = false;

    for(int w = 0; w <= w_max; w++) {
      double val_x = x[i + w];
      double val_y = y[i + tau + w];

      if (NumericVector::is_na(val_x) || NumericVector::is_na(val_y)) {
        if (!na_rm) {
          has_na = true;
          break;
        }
        continue; // pairwise: skip this pair
      }

      sum_x += val_x;
      sum_y += val_y;
      sum_x2 += val_x * val_x;
      sum_y2 += val_y * val_y;
      sum_xy += val_x * val_y;
      valid_n++;
    }

    if (!na_rm && has_na) {
      results[k] = NA_REAL;
      continue;
    }

    if (valid_n > 1) {
      double var_x = (sum_x2 - (sum_x * sum_x) / valid_n) / (valid_n - 1);
      double var_y = (sum_y2 - (sum_y * sum_y) / valid_n) / (valid_n - 1);
      double cov_xy = (sum_xy - (sum_x * sum_y) / valid_n) / (valid_n - 1);

      if (var_x <= 0 || var_y <= 0) {
        results[k] = NA_REAL;
      } else {
        results[k] = cov_xy / sqrt(var_x * var_y);
      }
    } else {
      results[k] = NA_REAL;
    }
  }

  return results;
}
