# Estimate an assumption-based probability-of-success proxy

Estimate an assumption-based probability-of-success proxy

## Usage

``` r
estimate_pos_proxy(
  soc_time,
  soc_survival,
  tau,
  target_delta_rmst,
  mean_log_hr,
  sd_log_hr,
  draws = 2000L,
  seed = NULL
)
```

## Arguments

- soc_time:

  Time grid for the SOC curve.

- soc_survival:

  SOC survival values at `soc_time`.

- tau:

  Truncation horizon.

- target_delta_rmst:

  Target RMST gain.

- mean_log_hr, sd_log_hr:

  Mean and standard deviation of the log-HR prior.

- draws:

  Number of Monte Carlo draws.

- seed:

  Optional seed.

## Value

A scalar probability.
