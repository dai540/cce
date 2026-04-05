# Build an analysis-ready dataset from normalized CCE tables

Build an analysis-ready dataset from normalized CCE tables

## Usage

``` r
build_analysis_dataset(
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

  Treatment-episode table with a single index row per patient.

- outcomes:

  Outcome table.

- biomarkers:

  Optional biomarker table.

- spec:

  A [`cce_spec()`](https://dai540.github.io/cce/reference/cce_spec.md)
  object.

## Value

A data frame with class `cce_dataset`.
