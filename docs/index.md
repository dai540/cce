# cce

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-315c86)](https://dai540.github.io/cce/)
[![R-tests](https://github.com/dai540/cce/actions/workflows/R-tests.yaml/badge.svg)](https://github.com/dai540/cce/actions/workflows/R-tests.yaml)
[![GitHub
release](https://img.shields.io/github/v/release/dai540/cce)](https://github.com/dai540/cce/releases)
[![License:
MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://dai540.github.io/cce/LICENSE)

`cce` is a source-first R package for counterfactual comparator work in
survival settings. It focuses on five jobs:

- validating normalized source tables
- building an analysis-ready cohort
- estimating VS-mode curves with `gformula`, `iptw_km`, and `iptw_cox`
- projecting SOC-only curves under hazard-ratio scenarios
- exporting machine-readable outputs for downstream review

The input contract is intentionally strict:

- `patient_baseline`: one row per patient with `index_date` and
  covariates
- `treatment_episodes`: one index-treatment row per patient
- `outcomes`: one endpoint row per patient for the selected endpoint
- `biomarkers`: optional baseline subgroup values

For repeatable analysis, `cce` ships a fixed
[`cce_spec()`](https://dai540.github.io/cce/reference/cce_spec.md)
contract, YAML round-tripping, dataset validation, profiling, and
bundled demo tables through
[`cce_demo_data()`](https://dai540.github.io/cce/reference/cce_demo_data.md).

The package validates inputs, records exclusions, estimates survival
effects, stores analysis metadata in `results.json`, and ships both demo
and real-data tutorials through pkgdown.

## Installation

`cce` is currently intended to be used from a source checkout. The
repository intentionally omits a maintainer email in `DESCRIPTION`, so
it is not set up for standard `install_github()` or CRAN-style
installation.

Clone the repository and load it in place with `pkgload`:

``` r
install.packages(c("pkgload", "testthat", "pkgdown", "rmarkdown"))
pkgload::load_all("path/to/cce", export_all = FALSE)
```

To build the documentation locally:

``` r
pkgdown::build_site("path/to/cce", install = FALSE, new_process = FALSE)
```

Then load the package:

``` r
pkgload::load_all("path/to/cce", export_all = FALSE)
```

## Minimal Example

[`fit_cce_vs()`](https://dai540.github.io/cce/reference/fit_cce_vs.md)
is the main comparative workflow:

``` r
demo <- cce::cce_demo_data(n = 220, seed = 7)

vs_fit <- cce::fit_cce_vs(
  data = demo$analysis_data,
  arm = "arm",
  time = "time",
  event = "event",
  covariates = c("age", "sex", "stage_or_risk", "ps"),
  subgroup = "subgroup",
  tau = 365,
  landmark_times = c(180, 365),
  bootstrap = 10,
  seed = 11
)

head(cce::as_effects_df(vs_fit))
```

For SOC-only planning:

``` r
soc_fit <- cce::project_soc_only(
  data = demo$analysis_data,
  arm = "arm",
  soc_level = "SOC",
  time = "time",
  event = "event",
  subgroup = "subgroup",
  tau = 365,
  hr_scenarios = c(0.65, 0.80, 1.00),
  target_delta_rmst = 30,
  prior_mean_log_hr = log(0.8),
  prior_sd_log_hr = 0.20,
  bootstrap = 10,
  seed = 99
)

head(cce::as_effects_df(soc_fit))
```

## Core Functions

- [`cce_spec()`](https://dai540.github.io/cce/reference/cce_spec.md)
- [`write_cce_spec()`](https://dai540.github.io/cce/reference/write_cce_spec.md)
- [`read_cce_spec()`](https://dai540.github.io/cce/reference/read_cce_spec.md)
- [`validate_cce_tables()`](https://dai540.github.io/cce/reference/validate_cce_tables.md)
- [`build_analysis_dataset()`](https://dai540.github.io/cce/reference/build_analysis_dataset.md)
- [`profile_cce_dataset()`](https://dai540.github.io/cce/reference/profile_cce_dataset.md)
- [`cce_demo_data()`](https://dai540.github.io/cce/reference/cce_demo_data.md)
- [`fit_cce_vs()`](https://dai540.github.io/cce/reference/fit_cce_vs.md)
- [`project_soc_only()`](https://dai540.github.io/cce/reference/project_soc_only.md)
- [`required_hr()`](https://dai540.github.io/cce/reference/required_hr.md)
- [`estimate_pos_proxy()`](https://dai540.github.io/cce/reference/estimate_pos_proxy.md)
- [`write_cce_results()`](https://dai540.github.io/cce/reference/write_cce_results.md)

## Main Outputs

Stable outputs include:

- `results.json`
- `curves.csv`
- `effects.csv`
- `diagnostics.csv`

The JSON payload includes:

- estimators and column mappings
- covariates and thresholds
- subgroup levels and scenario settings
- source spec, exclusions, validation report, and dataset profile

## Documentation

Website: <https://dai540.github.io/cce/>

Articles:

- Demo workflow:
  <https://dai540.github.io/cce/articles/demo-data-workflow.html>
- Real-data workflow:
  <https://dai540.github.io/cce/articles/public-oncology-data.html>

## Citation

``` r
citation("cce")
```
