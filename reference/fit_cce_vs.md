# Fit an observed two-arm survival comparison

Fits a simple observed survival comparison and returns curves, effect
tables, diagnostics, and metadata on a shared output contract.

## Usage

``` r
fit_cce_vs(
  data,
  arm = "arm",
  time = "time",
  event = "event",
  subgroup = NULL,
  tau,
  times = tau
)
```

## Arguments

- data:

  A data frame.

- arm:

  Arm column name.

- time:

  Follow-up time column name.

- event:

  Event indicator column name.

- subgroup:

  Optional subgroup column name.

- tau:

  Horizon for RMST calculation.

- times:

  Time points used for survival differences.

## Value

A `cce_fit` object.
