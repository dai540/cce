# Validate normalized CCE source tables without stopping the workflow

`validate_cce_tables()` inspects the normalized source tables required
by
[`build_analysis_dataset()`](https://dai540.github.io/cce/reference/build_analysis_dataset.md)
and returns a machine-readable issue table. The function is designed for
preflight checks before running a full analysis.

## Usage

``` r
validate_cce_tables(
  patient_baseline,
  treatment_episodes,
  outcomes,
  biomarkers = NULL,
  spec
)
```

## Arguments

- patient_baseline:

  Baseline patient table.

- treatment_episodes:

  Treatment-episode table.

- outcomes:

  Outcome table.

- biomarkers:

  Optional biomarker table.

- spec:

  A [`cce_spec()`](https://dai540.github.io/cce/reference/cce_spec.md)
  object.

## Value

A data frame with one row per issue.
