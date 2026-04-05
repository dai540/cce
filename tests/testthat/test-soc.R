test_that("SOC-only projections are monotone in HR", {
  demo <- cce_demo_data(n = 140, seed = 3)
  fit <- project_soc_only(
    data = demo$analysis_data,
    arm = "arm",
    soc_level = "SOC",
    time = "time",
    event = "event",
    subgroup = "subgroup",
    tau = 300,
    hr_scenarios = c(0.7, 0.9, 1.0),
    target_delta_rmst = 20,
    prior_mean_log_hr = log(0.85),
    prior_sd_log_hr = 0.2,
    bootstrap = 4,
    seed = 33
  )

  eff <- as_effects_df(fit)
  eff_all <- eff[eff$subgroup == "All" & eff$landmark_time == max(eff$landmark_time), ]
  eff_all <- eff_all[order(eff_all$scenario_hr), ]

  expect_s3_class(fit, "cce_soc_result")
  expect_true(all(diff(eff_all$delta_rmst) <= 0))
  expect_true(any(!is.na(eff$required_hr)))
  expect_true(any(!is.na(eff$pos_proxy)))
})
