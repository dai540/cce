test_that("validate_cce_tables returns warnings before build", {
  demo <- cce_demo_data(n = 100, seed = 21)
  demo$patient_baseline$ps[1:5] <- NA

  report <- validate_cce_tables(
    patient_baseline = demo$patient_baseline,
    treatment_episodes = demo$treatment_episodes,
    outcomes = demo$outcomes,
    biomarkers = demo$biomarkers,
    spec = demo$spec
  )

  expect_true(is.data.frame(report))
  expect_true(any(report$issue == "complete_case_exclusions"))
})

test_that("profile_cce_dataset returns summary tables", {
  demo <- cce_demo_data(n = 100, seed = 22)
  prof <- profile_cce_dataset(
    data = demo$analysis_data,
    arm = "arm",
    time = "time",
    event = "event",
    subgroup = "subgroup"
  )

  expect_s3_class(prof, "cce_profile")
  expect_true(all(c("overall", "by_arm", "by_subgroup", "by_arm_subgroup", "missingness") %in% names(prof)))
  expect_equal(prof$overall$n, nrow(demo$analysis_data))
})
