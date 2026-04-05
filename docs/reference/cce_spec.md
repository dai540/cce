# Create a Counterfactual Comparator Engine specification

`cce_spec()` defines the minimal schema needed to convert normalized
input tables into an analysis-ready dataset. The returned object is
intentionally simple and can be stored as YAML for reproducible runs.

## Usage

``` r
cce_spec(
  covariates,
  subgroup_biomarker = NULL,
  endpoint = "os",
  id_col = "patient_id",
  index_date_col = "index_date",
  regimen_col = "regimen_name",
  treatment_start_col = "start_date",
  index_flag_col = "is_index_treatment",
  endpoint_col = "endpoint",
  time_col = "time",
  event_col = "event",
  follow_up_col = "last_follow_up_date",
  biomarker_name_col = "biomarker_name",
  biomarker_value_col = "biomarker_value",
  biomarker_baseline_flag_col = "is_baseline",
  arm_map = c(SOC = "SOC", A = "A"),
  missing_strategy = "complete_case",
  time_zero_tolerance_days = 0L
)
```

## Arguments

- covariates:

  Character vector of baseline covariate column names.

- subgroup_biomarker:

  Optional biomarker name used to derive a subgroup column from the
  `biomarkers` table.

- endpoint:

  Endpoint name to retain from the `outcomes` table.

- id_col:

  Patient identifier column name.

- index_date_col:

  Baseline index-date column name.

- regimen_col:

  Treatment regimen column name.

- treatment_start_col:

  Index treatment start-date column name.

- index_flag_col:

  Logical column marking the index treatment row.

- endpoint_col:

  Outcome endpoint-name column.

- time_col:

  Follow-up time column.

- event_col:

  Event indicator column.

- follow_up_col:

  Last follow-up date column.

- biomarker_name_col:

  Biomarker-name column.

- biomarker_value_col:

  Biomarker-value column.

- biomarker_baseline_flag_col:

  Logical baseline-biomarker flag column.

- arm_map:

  Named character vector mapping raw regimen labels to output labels.
  The output labels must be exactly `c("SOC", "A")`.

- missing_strategy:

  Missing-data rule. Only `"complete_case"` is implemented in v0.1.0.

- time_zero_tolerance_days:

  Allowed difference between `index_date` and the index treatment start
  date.

## Value

A `cce_spec` object.
