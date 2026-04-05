test_that("VS mode returns expected tables", {
  demo <- cce_demo_data(n = 140, seed = 2)
  fit <- fit_cce_vs(
    data = demo$analysis_data,
    arm = "arm",
    time = "time",
    event = "event",
    covariates = c("age", "sex", "stage_or_risk", "ps"),
    subgroup = "subgroup",
    tau = 300,
    landmark_times = c(150, 300),
    bootstrap = 4,
    seed = 12
  )

  expect_s3_class(fit, "cce_vs_result")
  expect_true(all(c("method", "subgroup", "arm", "time", "survival") %in% names(as_curves_df(fit))))
  expect_true(all(c("delta_rmst", "delta_survival") %in% names(as_effects_df(fit))))
  expect_true(nrow(as_diagnostics_df(fit)) > 0)
  expect_true(all(is.finite(as_effects_df(fit)$delta_rmst)))
})
