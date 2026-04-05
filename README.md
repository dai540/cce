# cce

`cce` is an R package for counterfactual comparator work in survival settings.
It focuses on a practical workflow:

- validate normalized source tables
- build an analysis-ready cohort
- estimate adjusted VS-mode survival curves with g-formula, `iptw_km`, and `iptw_cox`
- project SOC-only curves under hazard-ratio scenarios
- export machine-readable results for downstream reporting

The package is designed for reproducibility-first development. It keeps inputs
explicit, returns tidy outputs, and ships two end-to-end tutorials:

- a synthetic demo-data workflow
- a public oncology dataset workflow using `survival::veteran`

Documentation site: <https://dai540.github.io/cce/>

The pkgdown site publishes both tutorials as HTML articles and mirrors the
function reference:

- Demo tutorial: <https://dai540.github.io/cce/articles/demo-data-workflow.html>
- Real-data tutorial: <https://dai540.github.io/cce/articles/public-oncology-data.html>

## Installation

From source:

```sh
R CMD INSTALL cce
```

## Quick start

```r
library(cce)

demo <- cce_demo_data(n = 220, seed = 7)
analysis <- demo$analysis_data

vs_fit <- fit_cce_vs(
  data = analysis,
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

head(as_effects_df(vs_fit))

soc_fit <- project_soc_only(
  data = analysis,
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

head(as_effects_df(soc_fit))
```

## Main functions

- `cce_spec()` creates a reusable schema contract
- `validate_cce_tables()` preflights normalized source tables and returns issues
- `build_analysis_dataset()` turns normalized tables into an analysis-ready set
- `profile_cce_dataset()` summarizes sample size, events, subgroup counts, and missingness
- `fit_cce_vs()` estimates g-formula, `iptw_km`, and `iptw_cox` comparator curves
- `project_soc_only()` runs assumption-based PH projections
- `write_cce_results()` writes `results.json`, `curves.csv`, `effects.csv`,
  and `diagnostics.csv`

## Tutorials

- `Demo-data workflow` walks through bundled normalized tables, cohort
  assembly, VS-mode estimation, SOC-only projection, and file export.
- `Public oncology data workflow` shows the same analysis pattern on the
  real patient-level `survival::veteran` dataset.

## Output contract

Both VS and SOC-only results include:

- tidy curve data
- effect summaries with RMST and landmark contrasts
- diagnostics
- machine-readable run metadata including covariates, thresholds, spec, exclusions, and dataset profile

SOC-only outputs are always labeled `Projection (assumption-based)`.
