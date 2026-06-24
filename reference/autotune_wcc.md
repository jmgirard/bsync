# Auto-Tune WCC Parameters for a Dataset

Automatically determines the optimal Windowed Cross-Correlation
parameters for a multi-dyad dataset by combining Power Spectral Density
(PSD) analysis with a surrogate-driven grid search.

## Usage

``` r
autotune_wcc(
  dyad_list,
  sample_rate,
  n_tune_dyads = 10,
  n_surrogates = 100,
  surrogate_method = c("phase", "circular"),
  trim_odd = FALSE,
  increment_pct = 0.05,
  window_multipliers = c(0.5, 1, 2),
  lag_multipliers = c(0.5, 1, 2),
  progress = TRUE
)
```

## Arguments

- dyad_list:

  A list of data frames, where each data frame represents a dyad and
  contains two numeric columns (the two time series).

- sample_rate:

  A single positive number indicating the sampling rate in Hertz.

- n_tune_dyads:

  Integer. The number of dyads to sample for the tuning phase. Default
  is 10 to provide a robust sample without excessive computation time.
  If the dataset has fewer than this number, all dyads are used.

- n_surrogates:

  Integer. Number of surrogates to generate per test. Default is 100,
  which provides a stable enough standard deviation to calculate
  standardized effect sizes during tuning.

- surrogate_method:

  Character string. "phase" (default) uses phase randomization, which
  preserves the power spectrum and is ideal for continuous physiological
  data. "circular" shifts the time series, which is better for
  preserving local autocorrelation in behavioral data.

- trim_odd:

  Logical. If \`TRUE\` and \`surrogate_method = "phase"\`, automatically
  drops the final observation of any odd-length time series to allow the
  Fourier transform to execute. Default is \`FALSE\`.

- increment_pct:

  Numeric value between 0.01 and 1.0. Determines the step size between
  successive windows as a percentage of the window size. Default is
  0.05.

- window_multipliers:

  A numeric vector. Multipliers applied to the baseline cycle length to
  generate the grid of window sizes. Default is \`c(0.5, 1.0, 2.0)\`.

- lag_multipliers:

  A numeric vector. Multipliers applied to the window size to generate
  the grid of maximum lags. Default is \`c(0.5, 1.0, 2.0)\`.

- progress:

  Logical. If \`TRUE\` (default), displays a dynamic progress bar in the
  console during the grid search.

## Value

A list containing the optimal parameters and the full tuning grid
results.

## Details

\*\*Reproducibility and Parallelization:\*\* This function involves
random sampling (selecting dyads and generating surrogate data). For
reproducible results, call \`set.seed()\` before running this function.
To speed up computation, ensure you have set a parallel backend using
the \`future\` package (e.g., \`future::plan(future::multisession)\`)
prior to execution.
