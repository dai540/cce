test_that("results are written to disk", {
  demo <- cce_demo_data(n = 100, seed = 4)
  fit <- fit_cce_vs(
    data = demo$analysis_data,
    arm = "arm",
    time = "time",
    event = "event",
    covariates = c("age", "sex", "stage_or_risk", "ps"),
    tau = 250,
    bootstrap = 0
  )

  out_dir <- tempfile(pattern = "cce-out-")
  write_cce_results(fit, out_dir)

  expect_true(file.exists(file.path(out_dir, "results.json")))
  expect_true(file.exists(file.path(out_dir, "curves.csv")))
  expect_true(file.exists(file.path(out_dir, "effects.csv")))
  expect_true(file.exists(file.path(out_dir, "diagnostics.csv")))
})
