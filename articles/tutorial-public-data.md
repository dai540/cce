# Tutorial: public survival data

## Use `survival::veteran`

[`survival::veteran`](https://rdrr.io/pkg/survival/man/veteran.html) is
a public dataset already distributed with the `survival` package, so
this tutorial needs no external download.

``` r
data(veteran, package = "survival")
#> Warning in data(veteran, package = "survival"): data set 'veteran' not found

veteran2 <- transform(
  veteran,
  arm = ifelse(trt == 1, "SOC", "A"),
  time = time,
  event = status,
  subgroup = as.character(celltype)
)

head(veteran2[, c("arm", "time", "event", "subgroup")])
#>   arm time event subgroup
#> 1 SOC   72     1 squamous
#> 2 SOC  411     1 squamous
#> 3 SOC  228     1 squamous
#> 4 SOC  126     1 squamous
#> 5 SOC  118     1 squamous
#> 6 SOC   10     1 squamous
```

## Run the observed comparison

``` r
vs_fit <- fit_cce_vs(
  data = veteran2,
  subgroup = "subgroup",
  tau = 180,
  times = c(90, 180)
)

head(as_effects_df(vs_fit))
#>   mode subgroup contrast         metric time     estimate
#> 1   vs      All  A - SOC delta_survival   90  -0.16657817
#> 2   vs      All  A - SOC delta_survival  180   0.02042615
#> 3   vs      All  A - SOC     delta_rmst  180  18.13311506
#> 4   vs    adeno  A - SOC delta_survival   90  -0.41666667
#> 5   vs    adeno  A - SOC delta_survival  180   0.06944444
#> 6   vs    adeno  A - SOC     delta_rmst  180 -10.83333333
```

## Run the SOC-only projection

``` r
soc_fit <- project_soc_only(
  data = veteran2,
  subgroup = "subgroup",
  tau = 180,
  hr_scenarios = c(0.75, 0.90, 1.00),
  target_delta_rmst = 10
)

head(as_effects_df(soc_fit))
#>       mode subgroup          contrast      metric time   estimate
#> 1 soc_only      All     HR_0.75 - SOC  delta_rmst  180 37.0727933
#> 2 soc_only      All     HR_0.90 - SOC  delta_rmst  180 12.8114755
#> 3 soc_only      All     HR_1.00 - SOC  delta_rmst  180  0.0000000
#> 4 soc_only      All required_hr - SOC required_hr  180  0.9204667
#> 5 soc_only    adeno     HR_0.75 - SOC  delta_rmst  180 13.5622187
#> 6 soc_only    adeno     HR_0.90 - SOC  delta_rmst  180  5.0077686
```
