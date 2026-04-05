# Profile an analysis-ready CCE dataset

`profile_cce_dataset()` returns compact summaries that are useful for
data QA, audit trails, and result metadata.

## Usage

``` r
profile_cce_dataset(
  data,
  arm = "arm",
  time = "time",
  event = "event",
  subgroup = NULL
)
```

## Arguments

- data:

  Analysis-ready data frame.

- arm, time, event:

  Column names for treatment, follow-up time, and event.

- subgroup:

  Optional subgroup column.

## Value

A list of summary data frames with class `cce_profile`.
