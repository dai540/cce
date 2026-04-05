# Estimate a counterfactual VS comparison

`fit_cce_vs()` fits two complementary estimators for a binary treatment
comparison: a Cox-model-based g-formula standardization, an
IPTW-weighted Kaplan-Meier curve (`iptw_km`), and an IPTW-weighted Cox
standardization (`iptw_cox`). The function returns tidy curves, effects,
diagnostics, and machine-readable metadata.

## Usage

``` r
fit_cce_vs(
  data,
  arm = "arm",
  time = "time",
  event = "event",
  covariates,
  subgroup = NULL,
  tau = NULL,
  landmark_times = NULL,
  n_grid = 100L,
  bootstrap = 0L,
  seed = 1L,
  weight_cap = 50,
  warning_max_weight = 10,
  fail_max_weight = 50,
  warning_smd = 0.1,
  fail_smd = 0.2
)
```

## Arguments

- data:

  Analysis-ready data frame.

- arm, time, event:

  Column names identifying treatment assignment, follow-up time, and
  event indicator.

- covariates:

  Character vector of baseline adjustment covariates.

- subgroup:

  Optional subgroup column name. When supplied, the result includes
  overall and subgroup-specific summaries.

- tau:

  RMST truncation horizon. Defaults to the 90th percentile of observed
  follow-up.

- landmark_times:

  Survival-difference time points.

- n_grid:

  Number of time points used to summarize curves.

- bootstrap:

  Number of bootstrap resamples used to derive intervals.

- seed:

  Random seed used for bootstrap resampling.

- weight_cap:

  Hard upper cap applied to stabilized IPTW weights.

- warning_max_weight, fail_max_weight:

  Diagnostic thresholds for the largest stabilized weight.

- warning_smd, fail_smd:

  Diagnostic thresholds for absolute SMD.

## Value

An object of class `cce_vs_result`.
