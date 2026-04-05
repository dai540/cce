# Project SOC-only curves under proportional hazards scenarios

Project SOC-only curves under proportional hazards scenarios

## Usage

``` r
project_soc_only(
  data,
  time = "time",
  event = "event",
  subgroup = NULL,
  arm = NULL,
  soc_level = "SOC",
  tau = NULL,
  landmark_times = NULL,
  hr_scenarios = c(0.7, 0.85, 1),
  target_delta_rmst = NULL,
  prior_mean_log_hr = NULL,
  prior_sd_log_hr = NULL,
  prior_draws = 2000L,
  bootstrap = 0L,
  seed = 1L,
  n_grid = 100L
)
```

## Arguments

- data:

  Analysis-ready data frame.

- time, event:

  Column names identifying follow-up and event status.

- subgroup:

  Optional subgroup column name.

- arm:

  Optional treatment column. When supplied, only rows matching
  `soc_level` are retained before projection.

- soc_level:

  Reference level used when `arm` is supplied.

- tau:

  RMST truncation horizon.

- landmark_times:

  Survival-difference time points.

- hr_scenarios:

  Numeric vector of proportional-hazard multipliers.

- target_delta_rmst:

  Optional target RMST gain used to reverse-solve `required_hr`.

- prior_mean_log_hr, prior_sd_log_hr:

  Optional normal prior parameters for the log hazard ratio. When both
  are supplied alongside `target_delta_rmst`, a PoS proxy is reported.

- prior_draws:

  Number of Monte Carlo draws for the PoS proxy.

- bootstrap:

  Number of bootstrap resamples used to derive intervals.

- seed:

  Random seed.

- n_grid:

  Number of time points used to summarize curves.

## Value

An object of class `cce_soc_result`.
