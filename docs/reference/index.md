# Package index

## Specifications and data assembly

- [`cce_spec()`](cce_spec.md) : Create a Counterfactual Comparator
  Engine specification
- [`read_cce_spec()`](read_cce_spec.md) : Read a CCE specification from
  YAML
- [`write_cce_spec()`](write_cce_spec.md) : Write a CCE specification to
  YAML
- [`build_analysis_dataset()`](build_analysis_dataset.md) : Build an
  analysis-ready dataset from normalized CCE tables
- [`cce_demo_data()`](cce_demo_data.md) : Generate bundled synthetic
  example data

## Estimation

- [`fit_cce_vs()`](fit_cce_vs.md) : Estimate a counterfactual VS
  comparison
- [`project_soc_only()`](project_soc_only.md) : Project SOC-only curves
  under proportional hazards scenarios
- [`required_hr()`](required_hr.md) : Reverse-solve the hazard ratio
  needed to hit an RMST target
- [`estimate_pos_proxy()`](estimate_pos_proxy.md) : Estimate an
  assumption-based probability-of-success proxy

## Outputs

- [`as_curves_df()`](as_curves_df.md) : Extract curves from a CCE result
- [`as_effects_df()`](as_effects_df.md) : Extract effect summaries from
  a CCE result
- [`as_diagnostics_df()`](as_diagnostics_df.md) : Extract diagnostics
  from a CCE result
- [`write_cce_results()`](write_cce_results.md) : Write CCE result files
  to disk
