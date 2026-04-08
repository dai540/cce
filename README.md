# cce

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-315c86)](https://dai540.github.io/cce/)
[![R-tests](https://github.com/dai540/cce/actions/workflows/R-tests.yaml/badge.svg)](https://github.com/dai540/cce/actions/workflows/R-tests.yaml)
[![GitHub release](https://img.shields.io/github/v/release/dai540/cce)](https://github.com/dai540/cce/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

`cce` is an R package for counterfactual comparator analysis in survival
settings. It is designed around two practical workflows:

<https://dai540.github.io/cce/>

- `fit_cce_vs()` for treatment A versus SOC comparisons on a common output contract
- `project_soc_only()` for assumption-based SOC-only scenario planning on the same output contract

The package is intentionally narrow. It does not try to be a full causal
inference framework. Instead, it standardizes the data contract, diagnostic
tables, survival curves, effect tables, and export files needed for
counterfactual comparator work in oncology-style survival analyses.

It is designed for analysts working on:

- external comparator work
- SOC benchmarking
- subgroup survival comparisons
- biomarker-oriented scenario planning
- machine-readable survival reporting

The input contract is intentionally strict:

- `patient_baseline`: one row per patient with `index_date` and baseline covariates
- `treatment_episodes`: one index-treatment row per patient after arm mapping
- `outcomes`: one endpoint row per patient for the selected time-to-event endpoint
- `biomarkers`: optional baseline subgroup values

For repeatable analysis, `cce` ships a fixed `cce_spec()` contract, YAML
round-tripping, validation helpers, dataset profiling, and bundled demo tables
through `cce_demo_data()`.

## Installation

`cce` is currently intended to be used from a source checkout. The repository
intentionally omits a maintainer email in `DESCRIPTION`, so it is not set up for
standard `install_github()` or CRAN-style installation.

Clone the repository and load it in place with `pkgload`:

```r
install.packages(c("pkgload", "testthat", "pkgdown", "rmarkdown"))
pkgload::load_all("path/to/cce", export_all = FALSE)
```

To build the documentation locally:

```r
pkgdown::build_site("path/to/cce", install = FALSE, new_process = FALSE)
```

Then load the package:

```r
pkgload::load_all("path/to/cce", export_all = FALSE)
```

## Core estimands

### VS workflow

`cce` estimates comparative survival summaries on a common output contract:

- standardized survival curves
- landmark survival differences
- restricted mean survival time differences
- diagnostic tables for weighting and overlap

### SOC-only workflow

`cce` projects scenario-based survival summaries from an SOC baseline:

- observed SOC survival curves
- projected comparator curves under hazard-ratio scenarios
- required hazard ratios for target RMST gains
- probability-of-success proxies under a log-HR prior

SOC-only outputs are always labeled as projection-based, not causal estimates.

## What cce does

`cce` does four things.

- Validates normalized source tables before analysis
- Builds a single analysis-ready dataset with exclusions and profiling attached
- Fits VS and SOC-only workflows on a shared output schema
- Exports curves, effects, diagnostics, and metadata for downstream reporting

In practice, the package is doing this:

- `validate_cce_tables()` checks schema completeness and arm mapping
- `build_analysis_dataset()` enforces index-treatment uniqueness and time-zero consistency
- `fit_cce_vs()` returns `gformula`, `iptw_km`, and `iptw_cox` results on one contract
- `project_soc_only()` returns projected curves, required HR values, and PoS proxies
- `write_cce_results()` writes:
  - `results.json`
  - `curves.csv`
  - `effects.csv`
  - `diagnostics.csv`

## Main functions

- `cce_spec()`
- `write_cce_spec()`
- `read_cce_spec()`
- `validate_cce_tables()`
- `build_analysis_dataset()`
- `profile_cce_dataset()`
- `cce_demo_data()`
- `fit_cce_vs()`
- `project_soc_only()`
- `required_hr()`
- `estimate_pos_proxy()`
- `write_cce_results()`

## Example

```r
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

## Built-in tutorial datasets

The package currently ships two tutorial pathways.

### Demo workflow

- `cce_demo_data()`
- fully synthetic normalized source tables plus an analysis dataset

### Public survival workflow

- `survival::veteran`
- compact public time-to-event data used as a real-data tutorial

## Tutorials

The tutorial site is organized around:

- `Get Started`
  - basic package workflow
- `Design`
  - target contract, estimands, diagnostics, and output rules
- `Case Studies`
  - bundled demo data
  - public `survival::veteran` workflow

## What cce cannot do yet

- It does not implement time-varying treatment strategies
- It does not implement dynamic treatment regimes
- It does not implement multi-arm comparisons on one contract
- It does not implement production-grade missing-data workflows beyond the current package interface
- It does not replace a formal target trial protocol or SAP

## Package layout

- `R/spec.R`: analysis specification contract and YAML helpers
- `R/data-build.R`: validation, cohort assembly, exclusions, and profiling
- `R/vs.R`: VS workflow estimators and diagnostics
- `R/soc_projection.R`: SOC-only projection workflow
- `R/outputs.R`: tidy exports and machine-readable output files
- `vignettes/`: getting-started, design notes, and case-study tutorials

## Documentation

Website: <https://dai540.github.io/cce/>

## Citation

```r
citation("cce")
```
