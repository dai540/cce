# Getting started

## Overview

`cce` keeps one small output contract across two workflows:

- observed two-arm comparison with
  [`fit_cce_vs()`](https://dai540.github.io/cce/reference/fit_cce_vs.md)
- SOC-only projection with
  [`project_soc_only()`](https://dai540.github.io/cce/reference/project_soc_only.md)

``` r
demo <- cce_demo_data(n = 200, seed = 1)
head(demo)
#>   arm time event       subgroup age ps
#> 1 SOC   83     1 Biomarker-high  69  2
#> 2   A  214     0     All-comers  53  1
#> 3 SOC  123     1 Biomarker-high  74  0
#> 4 SOC  405     0 Biomarker-high  64  1
#> 5   A  341     0 Biomarker-high  70  1
#> 6 SOC  280     1  Biomarker-low  73  0
```

## Observed comparison

``` r
vs_fit <- fit_cce_vs(
  data = demo,
  arm = "arm",
  time = "time",
  event = "event",
  subgroup = "subgroup",
  tau = 365,
  times = c(180, 365)
)

vs_fit
#> cce_fit
#>   mode: vs 
#>   label: Observed comparison 
#>   curves: 722 rows
#>   effects: 12 rows
#>   diagnostics: 16 rows
head(as_effects_df(vs_fit))
#>   mode   subgroup contrast         metric time    estimate
#> 1   vs        All  A - SOC delta_survival  180   0.2915253
#> 2   vs        All  A - SOC delta_survival  365   0.1938709
#> 3   vs        All  A - SOC     delta_rmst  365  99.2640724
#> 4   vs All-comers  A - SOC delta_survival  180   0.3017719
#> 5   vs All-comers  A - SOC delta_survival  365   0.2331357
#> 6   vs All-comers  A - SOC     delta_rmst  365 102.4089881
```

## SOC-only projection

``` r
soc_fit <- project_soc_only(
  data = demo,
  arm = "arm",
  time = "time",
  event = "event",
  subgroup = "subgroup",
  tau = 365,
  hr_scenarios = c(0.7, 0.85, 1.0),
  target_delta_rmst = 30
)

soc_fit
#> cce_fit
#>   mode: soc_only 
#>   label: Projection (assumption-based) 
#>   curves: 764 rows
#>   effects: 16 rows
#>   diagnostics: 8 rows
head(as_effects_df(soc_fit))
#>       mode   subgroup          contrast      metric time  estimate
#> 1 soc_only        All     HR_0.70 - SOC  delta_rmst  365 55.889011
#> 2 soc_only        All     HR_0.85 - SOC  delta_rmst  365 25.771121
#> 3 soc_only        All     HR_1.00 - SOC  delta_rmst  365  0.000000
#> 4 soc_only        All required_hr - SOC required_hr  365  0.827455
#> 5 soc_only All-comers     HR_0.70 - SOC  delta_rmst  365 51.605992
#> 6 soc_only All-comers     HR_0.85 - SOC  delta_rmst  365 23.531579
```

## Export outputs

``` r
out_dir <- tempfile("cce-site-")
write_cce_results(vs_fit, out_dir)
list.files(out_dir)
#> [1] "curves.csv"      "diagnostics.csv" "effects.csv"     "results.json"
```

The package always writes the same four files:

- `curves.csv`
- `effects.csv`
- `diagnostics.csv`
- `results.json`
