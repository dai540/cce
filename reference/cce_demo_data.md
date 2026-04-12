# Generate demo survival data

Generates a small synthetic survival dataset for examples and tests.

## Usage

``` r
cce_demo_data(n = 200, seed = 1)
```

## Arguments

- n:

  Number of rows.

- seed:

  Random seed.

## Value

A data frame with columns `arm`, `time`, `event`, `subgroup`, `age`, and
`ps`.
