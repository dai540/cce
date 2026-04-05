# Package index

## Specifications and data assembly

- [`cce_spec()`](https://dai540.github.io/cce/reference/cce_spec.md) :
  Create a Counterfactual Comparator Engine specification
- [`read_cce_spec()`](https://dai540.github.io/cce/reference/read_cce_spec.md)
  : Read a CCE specification from YAML
- [`write_cce_spec()`](https://dai540.github.io/cce/reference/write_cce_spec.md)
  : Write a CCE specification to YAML
- [`validate_cce_tables()`](https://dai540.github.io/cce/reference/validate_cce_tables.md)
  : Validate normalized CCE source tables without stopping the workflow
- [`profile_cce_dataset()`](https://dai540.github.io/cce/reference/profile_cce_dataset.md)
  : Profile an analysis-ready CCE dataset
- [`build_analysis_dataset()`](https://dai540.github.io/cce/reference/build_analysis_dataset.md)
  : Build an analysis-ready dataset from normalized CCE tables
- [`cce_demo_data()`](https://dai540.github.io/cce/reference/cce_demo_data.md)
  : Generate bundled synthetic example data

## Estimation

- [`fit_cce_vs()`](https://dai540.github.io/cce/reference/fit_cce_vs.md)
  : Estimate a counterfactual VS comparison
- [`project_soc_only()`](https://dai540.github.io/cce/reference/project_soc_only.md)
  : Project SOC-only curves under proportional hazards scenarios
- [`required_hr()`](https://dai540.github.io/cce/reference/required_hr.md)
  : Reverse-solve the hazard ratio needed to hit an RMST target
- [`estimate_pos_proxy()`](https://dai540.github.io/cce/reference/estimate_pos_proxy.md)
  : Estimate an assumption-based probability-of-success proxy

## Outputs

- [`as_curves_df()`](https://dai540.github.io/cce/reference/as_curves_df.md)
  : Extract curves from a CCE result
- [`as_effects_df()`](https://dai540.github.io/cce/reference/as_effects_df.md)
  : Extract effect summaries from a CCE result
- [`as_diagnostics_df()`](https://dai540.github.io/cce/reference/as_diagnostics_df.md)
  : Extract diagnostics from a CCE result
- [`write_cce_results()`](https://dai540.github.io/cce/reference/write_cce_results.md)
  : Write CCE result files to disk
