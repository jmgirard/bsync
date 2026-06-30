# Find Optimum (Peak or Valley) in Windowed Analyses

Find Optimum (Peak or Valley) in Windowed Analyses

## Usage

``` r
pick_optima(
  obj,
  L_size = NULL,
  strict_monotonic = FALSE,
  find_min = NULL,
  search_method = NULL,
  threshold = NULL
)
```

## Arguments

- obj:

  An object of class "wcc_res" or "wdtw_res".

- L_size:

  An odd integer specifying the size of the local search region. Ignored
  if \`search_method = "global"\`. Default is \`NULL\`.

- strict_monotonic:

  Logical indicating whether to strictly enforce monotonic flanks around
  the extremum. Ignored if \`search_method = "global"\`. Default is
  FALSE.

- find_min:

  Logical indicating whether to search for local minima instead of local
  maxima. If \`NULL\` (the default), the function automatically searches
  for maxima (\`FALSE\`) for cross-correlation ("wcc_res") and minima
  (\`TRUE\`) for distance metrics ("wdtw_res").

- search_method:

  Character string specifying "local" or "global" search. "local"
  searches symmetrically outward from lag 0. "global" searches the
  entire window for the absolute extremum. If \`NULL\`, defaults to
  "local" for "wcc_res" and "global" for "wdtw_res".

- threshold:

  A numeric value. For WCC (\`find_min = FALSE\`), optima with an
  absolute value below this threshold are set to NA. For WDTW
  (\`find_min = TRUE\`), optima with a distance above this threshold are
  set to NA. Default is \`NULL\`.

## Value

A data frame of class "wcc_optima" or "wdtw_optima".

## Examples

``` r
wcc_res <- wcc(sim_dyad$x_A, sim_dyad$x_B, window_size = 96, lag_max = 10)
optima <- pick_optima(wcc_res, L_size = 9)
head(optima)
#> 
#> ── WCC Optima Results ──────────────────────────────────────────────────────────
#> Total Windows Analyzed: 6
#> Valid Optima Found: 6 (100%)
#> Search Method: local
#> Search Mode: Peaks (Maxima)
#> Threshold Applied: None
#> Local Search Size: 9
#> Strict Monotonic: FALSE
#> Showing the first 5 results:
#>   i optimum_lag optimum_value
#>  11          -4     0.1454372
#>  12          -4     0.1496750
#>  13          -4     0.1502625
#>  14          -4     0.1434978
#>  15          -4     0.1506591
#> # ... with 1 more row
```
