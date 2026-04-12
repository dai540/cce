# Project SOC-only survival scenarios

Projects comparator curves from the observed `SOC` survival curve using
a proportional-hazards assumption.

## Usage

``` r
project_soc_only(
  data,
  arm = "arm",
  time = "time",
  event = "event",
  subgroup = NULL,
  soc_level = "SOC",
  tau,
  hr_scenarios = c(0.8, 1),
  target_delta_rmst = NULL
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

- soc_level:

  Label used for standard of care.

- tau:

  Horizon for RMST calculation.

- hr_scenarios:

  Numeric vector of hazard-ratio scenarios.

- target_delta_rmst:

  Optional target RMST gain used to back-solve a required hazard ratio.

## Value

A `cce_fit` object.
