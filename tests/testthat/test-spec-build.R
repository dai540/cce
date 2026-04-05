test_that("spec roundtrip and analysis dataset work", {
  demo <- cce_demo_data(n = 120, seed = 1)
  path <- tempfile(fileext = ".yml")
  write_cce_spec(demo$spec, path)
  spec2 <- read_cce_spec(path)

  rebuilt <- build_analysis_dataset(
    patient_baseline = demo$patient_baseline,
    treatment_episodes = demo$treatment_episodes,
    outcomes = demo$outcomes,
    biomarkers = demo$biomarkers,
    spec = spec2
  )

  expect_s3_class(rebuilt, "cce_dataset")
  expect_true(all(c("arm", "time", "event", "subgroup") %in% names(rebuilt)))
  expect_equal(sort(unique(as.character(rebuilt$arm))), c("A", "SOC"))
})
