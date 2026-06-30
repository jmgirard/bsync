# Auto-Tune WCC Parameters for a Multi-Dyad Dataset

Selects Windowed Cross-Correlation hyperparameters that are both
detectable (significant vs. the null) and stable (consistent) across a
collection of dyads. Internally calls \[synchrony_multiverse()\] on each
dyad and applies a gated stability-penalized selection rule via
\[select_specification()\].

## Usage

``` r
autotune_wcc(
  dyad_list,
  sample_rate,
  window_sec,
  lag_sec = NULL,
  increment_pct = 0.1,
  statistic = "mean_abs_z",
  surrogate_method = "phase",
  n_surrogates = 100L,
  n_tune_dyads = 30L,
  sig_pct = 0.5,
  iqr_penalty = 0.5
)
```

## Arguments

- dyad_list:

  A list of data frames or lists. Each element represents one dyad and
  must have at least two numeric columns (or two named list elements
  \`x\` and \`y\`) containing the two time series.

- sample_rate:

  Single positive number; sampling rate in Hz, used to convert
  \`window_sec\` and \`lag_sec\` to samples.

- window_sec:

  Numeric vector; window size(s) in seconds to sweep. Use
  \[suggest_wcc_params()\] on a representative dyad to find a principled
  starting range.

- lag_sec:

  Numeric vector; max lag(s) in seconds. Default \`NULL\` uses
  \`window_sec / 2\` per cell (the SUSY reliability ceiling).

- increment_pct:

  Numeric; window increment as a fraction of window size (e.g., \`0.1\`
  = 10% step). Default is \`0.1\`.

- statistic:

  Character; WCC aggregate statistic. Default \`"mean_abs_z"\`.

- surrogate_method:

  Character; surrogate generator: \`"phase"\` (default) or
  \`"circular"\`.

- n_surrogates:

  Single positive integer; surrogates per cell per dyad. Default
  \`100\`. Increase to \>= 1000 for reporting.

- n_tune_dyads:

  Maximum number of dyads to use. If \`length(dyad_list) \>
  n_tune_dyads\`, a random sample is taken. Default \`30\`.

- sig_pct:

  Detectability gate: minimum proportion of dyads in which a cell must
  be significant (p \< .05). Default \`0.5\`.

- iqr_penalty:

  Penalty weight on cross-dyad IQR of ES. Score = \`median(ES) -
  iqr_penalty \* IQR(ES)\`. Default \`0.5\`.

## Value

A named list with:

- \`window_size\`:

  Selected window size in samples.

- \`lag_max\`:

  Selected max lag in samples.

- \`window_increment\`:

  Selected window increment in samples.

- \`lag_increment\`:

  \`1L\` (standard lag increment).

- \`window_sec\`:

  Selected window size in seconds.

- \`lag_sec\`:

  Selected max lag in seconds.

- \`sig_rate\`:

  Proportion of dyads where selected cell was significant.

- \`median_es\`:

  Median ES across dyads for the selected cell.

- \`iqr_es\`:

  IQR of ES across dyads for the selected cell.

- \`score\`:

  Selection score for the chosen cell.

- \`n_dyads\`:

  Number of dyads used for tuning.

- \`n_cells_gated\`:

  Number of cells that passed the detectability gate.

- \`dyad_multiverses\`:

  List of \`bsync_multiverse\` objects, one per dyad.

## Details

\*\*Why cross-dyad stability?\*\* A parameter set that maximizes raw
synchrony for one dyad may simply match that dyad's autocorrelation
structure. The matched-null surrogate controls for autocorrelation
within a dyad (Invariant 2), but the \*best\* parameters should also
replicate across dyads with structurally different signals – hence the
multi-dyad stability criterion.

\*\*Selection rule.\*\* Cells pass a detectability gate (significant in
at least \`sig_pct\` of dyads). Among passing cells, the score is
\`median(ES) - iqr_penalty \* IQR(ES)\` across dyads, penalizing spread.
If no cell passes the gate, a warning is issued and the
highest-median-ES cell is returned (soft fallback).

\*\*Dyad sampling.\*\* If \`length(dyad_list) \> n_tune_dyads\`, a
random sample of \`n_tune_dyads\` dyads is used for speed; call
\`set.seed()\` beforehand for reproducibility.

## See also

\[synchrony_multiverse()\], \[suggest_wcc_params()\],
\[select_specification()\]
