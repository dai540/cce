# Tutorial: demo data

## Build a compact analysis

``` r
demo <- cce_demo_data(n = 240, seed = 10)
table(demo$arm, demo$subgroup)
#>      
#>       All-comers Biomarker-high Biomarker-low
#>   A           41             45            33
#>   SOC         43             38            40
```

## Run the observed comparison

``` r
vs_fit <- fit_cce_vs(
  data = demo,
  subgroup = "subgroup",
  tau = 365,
  times = c(180, 365)
)

head(as_curves_df(vs_fit))
#>   mode subgroup method arm time  survival n_risk
#> 1   vs      All     km   A    0 1.0000000    119
#> 2   vs      All     km   A    2 1.0000000    119
#> 3   vs      All     km   A    4 0.9915966    119
#> 4   vs      All     km   A    9 0.9831933    118
#> 5   vs      All     km   A   12 0.9663866    117
#> 6   vs      All     km   A   16 0.9579832    115
head(as_diagnostics_df(vs_fit))
#>   mode   subgroup arm metric value
#> 1   vs        All   A      n   119
#> 2   vs        All   A events    59
#> 3   vs        All SOC      n   121
#> 4   vs        All SOC events    78
#> 5   vs All-comers   A      n    41
#> 6   vs All-comers   A events    24
```

## Run SOC-only planning

``` r
soc_fit <- project_soc_only(
  data = demo,
  subgroup = "subgroup",
  tau = 365,
  hr_scenarios = c(0.65, 0.80, 1.00),
  target_delta_rmst = 25
)

head(as_effects_df(soc_fit))
#>       mode   subgroup          contrast      metric time   estimate
#> 1 soc_only        All     HR_0.65 - SOC  delta_rmst  365 65.9071438
#> 2 soc_only        All     HR_0.80 - SOC  delta_rmst  365 34.6511694
#> 3 soc_only        All     HR_1.00 - SOC  delta_rmst  365  0.0000000
#> 4 soc_only        All required_hr - SOC required_hr  365  0.8516335
#> 5 soc_only All-comers     HR_0.65 - SOC  delta_rmst  365 49.1891369
#> 6 soc_only All-comers     HR_0.80 - SOC  delta_rmst  365 25.8445192
```

## Write results

``` r
out_dir <- tempfile("cce-demo-")
write_cce_results(soc_fit, out_dir)
list.files(out_dir)
#> [1] "curves.csv"      "diagnostics.csv" "effects.csv"     "results.json"
```
