test_that("demo data has expected columns", {
  x <- cce_demo_data(n = 50, seed = 1)
  expect_true(all(c("arm", "time", "event", "subgroup") %in% names(x)))
  expect_equal(length(unique(x$arm)), 2)
})

test_that("vs workflow returns tables", {
  x <- cce_demo_data(n = 80, seed = 2)
  fit <- fit_cce_vs(x, tau = 365, times = c(180, 365), subgroup = "subgroup")
  expect_s3_class(fit, "cce_fit")
  expect_true(nrow(as_curves_df(fit)) > 0)
  expect_true(nrow(as_effects_df(fit)) > 0)
  expect_true(nrow(as_diagnostics_df(fit)) > 0)
})

test_that("soc-only workflow returns tables", {
  x <- cce_demo_data(n = 80, seed = 3)
  fit <- project_soc_only(
    x,
    tau = 365,
    subgroup = "subgroup",
    hr_scenarios = c(0.7, 0.9, 1.0),
    target_delta_rmst = 20
  )
  expect_s3_class(fit, "cce_fit")
  expect_true(nrow(as_curves_df(fit)) > 0)
  expect_true(any(as_effects_df(fit)$metric == "required_hr"))
})

test_that("results writer creates expected files", {
  x <- cce_demo_data(n = 60, seed = 4)
  fit <- fit_cce_vs(x, tau = 365, times = c(180, 365))
  out_dir <- tempfile("cce-out-")
  write_cce_results(fit, out_dir)
  expect_true(all(c("curves.csv", "effects.csv", "diagnostics.csv", "results.json") %in% list.files(out_dir)))
})
