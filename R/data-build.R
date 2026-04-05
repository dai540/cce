#' Build an analysis-ready dataset from normalized CCE tables
#'
#' @param patient_baseline Baseline patient table.
#' @param treatment_episodes Treatment-episode table with a single index row
#'   per patient.
#' @param outcomes Outcome table.
#' @param biomarkers Optional biomarker table.
#' @param spec A [cce_spec()] object.
#'
#' @return A data frame with class `cce_dataset`.
#' @export
build_analysis_dataset <- function(
    patient_baseline,
    treatment_episodes,
    outcomes,
    biomarkers = NULL,
    spec) {
  assert_data_frame(patient_baseline, "patient_baseline")
  assert_data_frame(treatment_episodes, "treatment_episodes")
  assert_data_frame(outcomes, "outcomes")
  if (!is.null(biomarkers)) {
    assert_data_frame(biomarkers, "biomarkers")
  }
  if (!inherits(spec, "cce_spec")) {
    stop("`spec` must inherit from `cce_spec`.", call. = FALSE)
  }

  req_patient <- c(spec$id_col, spec$index_date_col, spec$covariates)
  req_treat <- c(spec$id_col, spec$regimen_col, spec$treatment_start_col, spec$index_flag_col)
  req_outcomes <- c(spec$id_col, spec$endpoint_col, spec$time_col, spec$event_col, spec$follow_up_col)
  assert_named_columns(patient_baseline, req_patient, "patient_baseline")
  assert_named_columns(treatment_episodes, req_treat, "treatment_episodes")
  assert_named_columns(outcomes, req_outcomes, "outcomes")
  if (!is.null(biomarkers) && !is.null(spec$subgroup_biomarker)) {
    assert_named_columns(biomarkers, c(
      spec$id_col,
      spec$biomarker_name_col,
      spec$biomarker_value_col,
      spec$biomarker_baseline_flag_col
    ), "biomarkers")
  }

  id <- spec$id_col
  patient_ids <- patient_baseline[[id]]
  if (anyDuplicated(patient_ids)) {
    stop("`patient_baseline` must have one row per patient.", call. = FALSE)
  }

  treat <- treatment_episodes[treatment_episodes[[spec$index_flag_col]] %in% TRUE, , drop = FALSE]
  if (nrow(treat) == 0L) {
    stop("No index treatments were found.", call. = FALSE)
  }
  if (anyDuplicated(treat[[id]])) {
    stop("`treatment_episodes` must contain exactly one index treatment per patient.", call. = FALSE)
  }
  treat <- treat[treat[[spec$regimen_col]] %in% names(spec$arm_map), , drop = FALSE]
  if (nrow(treat) == 0L) {
    stop("No treatment rows matched `spec$arm_map`.", call. = FALSE)
  }
  treat$arm <- unname(spec$arm_map[as.character(treat[[spec$regimen_col]])])
  treat$arm <- factor(treat$arm, levels = c("SOC", "A"))

  outcome_keep <- outcomes[outcomes[[spec$endpoint_col]] == spec$endpoint, , drop = FALSE]
  if (nrow(outcome_keep) == 0L) {
    stop("No outcomes matched `spec$endpoint`.", call. = FALSE)
  }
  if (anyDuplicated(outcome_keep[[id]])) {
    stop("`outcomes` must contain at most one row per patient for the selected endpoint.", call. = FALSE)
  }

  merged <- merge(patient_baseline, treat[, c(id, "arm", spec$treatment_start_col), drop = FALSE],
    by = id, all = FALSE
  )
  merged <- merge(merged, outcome_keep[, c(id, spec$time_col, spec$event_col, spec$follow_up_col), drop = FALSE],
    by = id, all = FALSE
  )

  index_date <- normalize_date(merged[[spec$index_date_col]])
  start_date <- normalize_date(merged[[spec$treatment_start_col]])
  if (anyNA(index_date) || anyNA(start_date)) {
    stop("`index_date` and treatment start dates must be coercible to Date.", call. = FALSE)
  }
  offset_days <- abs(as.integer(start_date - index_date))
  if (any(offset_days > spec$time_zero_tolerance_days)) {
    stop("Time-zero validation failed: index_date and treatment start are inconsistent.", call. = FALSE)
  }

  if (!is.null(spec$subgroup_biomarker)) {
    if (is.null(biomarkers)) {
      stop("`biomarkers` is required when `spec$subgroup_biomarker` is set.", call. = FALSE)
    }
    bm <- biomarkers[
      biomarkers[[spec$biomarker_name_col]] == spec$subgroup_biomarker &
        biomarkers[[spec$biomarker_baseline_flag_col]] %in% TRUE,
      ,
      drop = FALSE
    ]
    if ("assay_time" %in% colnames(bm)) {
      assay_time <- suppressWarnings(as.POSIXct(bm$assay_time, tz = "UTC"))
      assay_time[is.na(assay_time)] <- as.POSIXct("1900-01-01", tz = "UTC")
      bm <- bm[order(bm[[id]], assay_time, decreasing = TRUE), , drop = FALSE]
    }
    bm <- bm[!duplicated(bm[[id]]), c(id, spec$biomarker_value_col), drop = FALSE]
    names(bm)[names(bm) == spec$biomarker_value_col] <- "subgroup"
    merged <- merge(merged, bm, by = id, all.x = TRUE)
  } else {
    merged$subgroup <- "All"
  }

  merged[[spec$event_col]] <- as.integer(merged[[spec$event_col]])
  merged[[spec$time_col]] <- as.numeric(merged[[spec$time_col]])
  if (any(!merged[[spec$event_col]] %in% c(0L, 1L), na.rm = TRUE)) {
    stop("`event` must contain 0/1 values.", call. = FALSE)
  }
  if (any(merged[[spec$time_col]] < 0, na.rm = TRUE)) {
    stop("`time` must be non-negative.", call. = FALSE)
  }

  needed <- c("arm", spec$time_col, spec$event_col, spec$covariates, "subgroup")
  keep <- stats::complete.cases(merged[, needed, drop = FALSE])
  exclusion <- data.frame(
    reason = c("missing_required_fields"),
    n = c(sum(!keep)),
    stringsAsFactors = FALSE
  )
  if (identical(spec$missing_strategy, "complete_case")) {
    merged <- merged[keep, , drop = FALSE]
  } else {
    stop("Only `complete_case` missing-data handling is implemented.", call. = FALSE)
  }
  rownames(merged) <- NULL
  attr(merged, "spec") <- spec
  attr(merged, "exclusions") <- exclusion
  class(merged) <- c("cce_dataset", class(merged))
  merged
}

#' Generate bundled synthetic example data
#'
#' @param n Number of synthetic patients.
#' @param seed Random seed.
#'
#' @return A list with normalized source tables and an analysis-ready dataset.
#' @export
cce_demo_data <- function(n = 250L, seed = 42L) {
  assert_scalar_numeric(n, "n", lower = 50)
  assert_scalar_numeric(seed, "seed", lower = 1)
  set.seed(seed)

  patient_id <- sprintf("P%04d", seq_len(n))
  index_date <- as.Date("2022-01-01") + sample(0:180, n, replace = TRUE)
  age <- round(stats::rnorm(n, mean = 64, sd = 9))
  sex <- sample(c("F", "M"), n, replace = TRUE, prob = c(0.4, 0.6))
  histology <- sample(c("Adenocarcinoma", "Squamous"), n, replace = TRUE, prob = c(0.65, 0.35))
  stage_or_risk <- sample(c("III", "IV"), n, replace = TRUE, prob = c(0.25, 0.75))
  ps <- sample(0:2, n, replace = TRUE, prob = c(0.45, 0.4, 0.15))
  bm_status <- sample(c("Low", "High"), n, replace = TRUE, prob = c(0.55, 0.45))

  lin_ps <- -0.5 + 0.6 * (ps >= 1) + 0.4 * (stage_or_risk == "IV") - 0.5 * (bm_status == "High")
  trt_prob <- stats::plogis(lin_ps)
  regimen_name <- ifelse(stats::runif(n) < trt_prob, "A", "SOC")

  base_hazard <- 1 / 420
  lp_surv <- 0.3 * (stage_or_risk == "IV") + 0.25 * ps - 0.2 * (bm_status == "High") -
    0.35 * (regimen_name == "A") - 0.15 * (regimen_name == "A" & bm_status == "High")
  event_time <- stats::rexp(n, rate = base_hazard * exp(lp_surv))
  censor_time <- stats::runif(n, min = 180, max = 720)
  time <- pmin(event_time, censor_time)
  event <- as.integer(event_time <= censor_time)
  time <- ceiling(time)

  patient_baseline <- data.frame(
    patient_id = patient_id,
    index_date = index_date,
    age = age,
    sex = sex,
    histology = histology,
    stage_or_risk = stage_or_risk,
    ps = ps,
    stringsAsFactors = FALSE
  )

  treatment_episodes <- data.frame(
    patient_id = patient_id,
    regimen_name = regimen_name,
    start_date = index_date,
    end_date = index_date + pmin(time, 180),
    line_of_therapy = 1L,
    intent = "palliative",
    is_index_treatment = TRUE,
    stringsAsFactors = FALSE
  )

  outcomes <- data.frame(
    patient_id = patient_id,
    endpoint = "os",
    time = time,
    event = event,
    last_follow_up_date = index_date + time,
    stringsAsFactors = FALSE
  )

  biomarkers <- data.frame(
    patient_id = patient_id,
    biomarker_name = "PDL1",
    biomarker_value = bm_status,
    assay_time = index_date - sample(0:30, n, replace = TRUE),
    is_baseline = TRUE,
    stringsAsFactors = FALSE
  )

  spec <- cce_spec(
    covariates = c("age", "sex", "stage_or_risk", "ps"),
    subgroup_biomarker = "PDL1"
  )
  analysis <- build_analysis_dataset(
    patient_baseline = patient_baseline,
    treatment_episodes = treatment_episodes,
    outcomes = outcomes,
    biomarkers = biomarkers,
    spec = spec
  )

  list(
    patient_baseline = patient_baseline,
    treatment_episodes = treatment_episodes,
    outcomes = outcomes,
    biomarkers = biomarkers,
    spec = spec,
    analysis_data = analysis
  )
}
