# Reverse-solve the hazard ratio needed to hit an RMST target

Reverse-solve the hazard ratio needed to hit an RMST target

## Usage

``` r
required_hr(
  soc_time,
  soc_survival,
  tau,
  target_delta_rmst,
  lower = 0.05,
  upper = 2.5
)
```

## Arguments

- soc_time:

  Time grid for the SOC survival curve.

- soc_survival:

  SOC survival values at `soc_time`.

- tau:

  Truncation horizon used for RMST.

- target_delta_rmst:

  Target RMST gain relative to SOC.

- lower, upper:

  Search interval for the hazard ratio.

## Value

A scalar hazard ratio or `NA_real_` when no root is bracketed.
