# Print method for bsync_multiverse objects

Print method for bsync_multiverse objects

## Usage

``` r
# S3 method for class 'bsync_multiverse'
print(x, direction = c("xy", "yx"), ...)
```

## Arguments

- x:

  A `bsync_multiverse` object.

- direction:

  For a Granger result (`estimator = "wgranger"`), the direction to
  report: `"xy"` (x -\> y, the default) or `"yx"` (y -\> x). Other
  estimators have one direction only, so `"yx"` is an error for them.

- ...:

  Additional arguments (not used).

## Value

Returns `x` invisibly.
