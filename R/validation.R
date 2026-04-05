new_issue_row <- function(component, severity, issue, n = NA_real_) {
  data.frame(
    component = component,
    severity = severity,
    issue = issue,
    n = n,
    stringsAsFactors = FALSE
  )
}

summarise_event_groups <- function(event, groups, group_names) {
  split_event <- split(event, groups, drop = TRUE)
  if (length(split_event) == 0L) {
    out <- as.data.frame(setNames(vector("list", length(group_names) + 3L), c(group_names, "n", "events", "event_rate")))
    return(out[0, , drop = FALSE])
  }
  keys <- strsplit(names(split_event), split = "\r", fixed = TRUE)
  do.call(rbind, lapply(seq_along(split_event), function(i) {
    key_vals <- keys[[i]]
    row_vals <- c(
      stats::setNames(as.list(key_vals), group_names),
      list(
        n = length(split_event[[i]]),
        events = sum(split_event[[i]], na.rm = TRUE),
        event_rate = mean(split_event[[i]], na.rm = TRUE)
      )
    )
    as.data.frame(row_vals, stringsAsFactors = FALSE)
  }))
}

empty_issue_table <- function() {
  data.frame(
    component = character(0),
    severity = character(0),
    issue = character(0),
    n = numeric(0),
    stringsAsFactors = FALSE
  )
}

append_issue <- function(issues, component, severity, issue, n = NA_real_) {
  c(issues, list(new_issue_row(component, severity, issue, n)))
}

#' Validate normalized CCE source tables without stopping the workflow
#'
#' `validate_cce_tables()` inspects the normalized source tables required by
#' [build_analysis_dataset()] and returns a machine-readable issue table. The
#' function is designed for preflight checks before running a full analysis.
#'
#' @param patient_baseline Baseline patient table.
#' @param treatment_episodes Treatment-episode table.
#' @param outcomes Outcome table.
#' @param biomarkers Optional biomarker table.
#' @param spec A [cce_spec()] object.
#'
#' @return A data frame with one row per issue.
#' @export
validate_cce_tables <- function(
    patient_baseline,
    treatment_episodes,
    outcomes,
    biomarkers = NULL,
    spec) {
  if (!inherits(spec, "cce_spec")) {
    stop("`spec` must inherit from `cce_spec`.", call. = FALSE)
  }

  issues <- list()

  if (!is.data.frame(patient_baseline)) {
    issues <- append_issue(issues, "patient_baseline", "fatal", "not_a_data_frame")
  }
  if (!is.data.frame(treatment_episodes)) {
    issues <- append_issue(issues, "treatment_episodes", "fatal", "not_a_data_frame")
  }
  if (!is.data.frame(outcomes)) {
    issues <- append_issue(issues, "outcomes", "fatal", "not_a_data_frame")
  }
  if (!is.null(spec$subgroup_biomarker) && is.null(biomarkers)) {
    issues <- append_issue(issues, "biomarkers", "fatal", "required_biomarker_table_missing")
  }
  if (!is.null(biomarkers) && !is.data.frame(biomarkers)) {
    issues <- append_issue(issues, "biomarkers", "fatal", "not_a_data_frame")
  }

  if (length(issues) > 0L) {
    return(do.call(rbind, issues))
  }

  req_patient <- c(spec$id_col, spec$index_date_col, spec$covariates)
  req_treat <- c(spec$id_col, spec$regimen_col, spec$treatment_start_col, spec$index_flag_col)
  req_outcomes <- c(spec$id_col, spec$endpoint_col, spec$time_col, spec$event_col, spec$follow_up_col)
  req_biomarkers <- c(spec$id_col, spec$biomarker_name_col, spec$biomarker_value_col, spec$biomarker_baseline_flag_col)

  missing_patient <- setdiff(req_patient, names(patient_baseline))
  missing_treat <- setdiff(req_treat, names(treatment_episodes))
  missing_outcomes <- setdiff(req_outcomes, names(outcomes))
  missing_biomarkers <- if (!is.null(spec$subgroup_biomarker) && !is.null(biomarkers)) setdiff(req_biomarkers, names(biomarkers)) else character(0)

  if (length(missing_patient) > 0L) {
    issues <- append_issue(issues, "patient_baseline", "fatal", paste("missing_columns:", paste(missing_patient, collapse = ",")), length(missing_patient))
  }
  if (length(missing_treat) > 0L) {
    issues <- append_issue(issues, "treatment_episodes", "fatal", paste("missing_columns:", paste(missing_treat, collapse = ",")), length(missing_treat))
  }
  if (length(missing_outcomes) > 0L) {
    issues <- append_issue(issues, "outcomes", "fatal", paste("missing_columns:", paste(missing_outcomes, collapse = ",")), length(missing_outcomes))
  }
  if (length(missing_biomarkers) > 0L) {
    issues <- append_issue(issues, "biomarkers", "fatal", paste("missing_columns:", paste(missing_biomarkers, collapse = ",")), length(missing_biomarkers))
  }

  if (length(issues) > 0L) {
    return(do.call(rbind, issues))
  }

  id <- spec$id_col
  dup_patient <- sum(duplicated(patient_baseline[[id]]))
  if (dup_patient > 0L) {
    issues <- append_issue(issues, "patient_baseline", "fatal", "duplicate_patient_ids", dup_patient)
  }

  treat <- treatment_episodes[treatment_episodes[[spec$index_flag_col]] %in% TRUE, , drop = FALSE]
  if (nrow(treat) == 0L) {
    issues <- append_issue(issues, "treatment_episodes", "fatal", "no_index_treatments")
  } else {
    dup_treat <- sum(duplicated(treat[[id]]))
    if (dup_treat > 0L) {
      issues <- append_issue(issues, "treatment_episodes", "fatal", "multiple_index_treatments_per_patient", dup_treat)
    }
    matched_treat <- treat[treat[[spec$regimen_col]] %in% names(spec$arm_map), , drop = FALSE]
    if (nrow(matched_treat) == 0L) {
      issues <- append_issue(issues, "treatment_episodes", "fatal", "no_rows_match_arm_map")
    }
  }

  outcome_keep <- outcomes[outcomes[[spec$endpoint_col]] == spec$endpoint, , drop = FALSE]
  if (nrow(outcome_keep) == 0L) {
    issues <- append_issue(issues, "outcomes", "fatal", "no_rows_match_endpoint")
  } else {
    dup_out <- sum(duplicated(outcome_keep[[id]]))
    if (dup_out > 0L) {
      issues <- append_issue(issues, "outcomes", "fatal", "duplicate_endpoint_rows_per_patient", dup_out)
    }
    non_binary_events <- sum(!outcome_keep[[spec$event_col]] %in% c(0L, 1L), na.rm = TRUE)
    if (non_binary_events > 0L) {
      issues <- append_issue(issues, "outcomes", "fatal", "non_binary_event_values", non_binary_events)
    }
    negative_time <- sum(outcome_keep[[spec$time_col]] < 0, na.rm = TRUE)
    if (negative_time > 0L) {
      issues <- append_issue(issues, "outcomes", "fatal", "negative_follow_up_times", negative_time)
    }
  }

  if (length(issues) == 0L) {
    merged <- merge(patient_baseline, treat[, c(id, spec$treatment_start_col), drop = FALSE], by = id, all = FALSE)
    merged <- merge(merged, outcome_keep[, c(id, spec$time_col, spec$event_col), drop = FALSE], by = id, all = FALSE)
    index_date <- normalize_date(merged[[spec$index_date_col]])
    start_date <- normalize_date(merged[[spec$treatment_start_col]])
    bad_dates <- sum(is.na(index_date) | is.na(start_date))
    if (bad_dates > 0L) {
      issues <- append_issue(issues, "patient_baseline", "fatal", "invalid_index_or_start_dates", bad_dates)
    } else {
      time_zero_fail <- sum(abs(as.integer(start_date - index_date)) > spec$time_zero_tolerance_days)
      if (time_zero_fail > 0L) {
        issues <- append_issue(issues, "cohort", "fatal", "time_zero_mismatch", time_zero_fail)
      }
    }
    needed <- c("arm", spec$time_col, spec$event_col, spec$covariates, "subgroup")
    merged$arm <- unname(spec$arm_map[as.character(treat[[spec$regimen_col]])])[match(merged[[id]], treat[[id]])]
    merged$subgroup <- "All"
    if (!is.null(spec$subgroup_biomarker) && !is.null(biomarkers)) {
      bm <- biomarkers[
        biomarkers[[spec$biomarker_name_col]] == spec$subgroup_biomarker &
          biomarkers[[spec$biomarker_baseline_flag_col]] %in% TRUE,
        c(id, spec$biomarker_value_col),
        drop = FALSE
      ]
      names(bm)[2L] <- "subgroup"
      bm <- bm[!duplicated(bm[[id]]), , drop = FALSE]
      merged <- merge(merged, bm, by = id, all.x = TRUE, suffixes = c("", ".bm"))
      if ("subgroup.bm" %in% names(merged)) {
        merged$subgroup <- merged$subgroup.bm
        merged$subgroup.bm <- NULL
      }
      missing_subgroup <- sum(is.na(merged$subgroup))
      if (missing_subgroup > 0L) {
        issues <- append_issue(issues, "biomarkers", "warning", "missing_baseline_subgroup_values", missing_subgroup)
      }
    }
    complete_case_drop <- sum(!stats::complete.cases(merged[, needed, drop = FALSE]))
    if (complete_case_drop > 0L) {
      issues <- append_issue(issues, "cohort", "warning", "complete_case_exclusions", complete_case_drop)
    }
  }

  if (length(issues) == 0L) {
    empty_issue_table()
  } else {
    do.call(rbind, issues)
  }
}

#' Profile an analysis-ready CCE dataset
#'
#' `profile_cce_dataset()` returns compact summaries that are useful for data
#' QA, audit trails, and result metadata.
#'
#' @param data Analysis-ready data frame.
#' @param arm,time,event Column names for treatment, follow-up time, and event.
#' @param subgroup Optional subgroup column.
#'
#' @return A list of summary data frames with class `cce_profile`.
#' @export
profile_cce_dataset <- function(
    data,
    arm = "arm",
    time = "time",
    event = "event",
    subgroup = NULL) {
  assert_data_frame(data, "data")
  required_cols <- c(time, event)
  if (!is.null(arm)) {
    required_cols <- c(required_cols, arm)
  }
  assert_named_columns(data, required_cols, "data")
  if (!is.null(subgroup)) {
    assert_named_columns(data, subgroup, "data")
  }

  overall <- data.frame(
    n = nrow(data),
    events = sum(data[[event]], na.rm = TRUE),
    event_rate = mean(data[[event]], na.rm = TRUE),
    median_follow_up = stats::median(data[[time]], na.rm = TRUE),
    max_follow_up = max(data[[time]], na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (!is.null(arm)) {
    arm_summary <- summarise_event_groups(
      event = data[[event]],
      groups = as.character(data[[arm]]),
      group_names = "arm"
    )
  } else {
    arm_summary <- data.frame(
      arm = "All",
      n = nrow(data),
      events = sum(data[[event]], na.rm = TRUE),
      event_rate = mean(data[[event]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  missingness <- data.frame(
    column = names(data),
    missing_n = vapply(data, function(x) sum(is.na(x)), numeric(1)),
    missing_rate = vapply(data, function(x) mean(is.na(x)), numeric(1)),
    stringsAsFactors = FALSE
  )

  subgroup_summary <- empty_issue_table()
  arm_subgroup_summary <- empty_issue_table()
  if (!is.null(subgroup)) {
    subgroup_summary <- summarise_event_groups(
      event = data[[event]],
      groups = as.character(data[[subgroup]]),
      group_names = "subgroup"
    )

    if (!is.null(arm)) {
      combo_key <- paste(as.character(data[[arm]]), as.character(data[[subgroup]]), sep = "\r")
      arm_subgroup_summary <- summarise_event_groups(
        event = data[[event]],
        groups = combo_key,
        group_names = c("arm", "subgroup")
      )
    } else {
      arm_subgroup_summary <- data.frame(
        arm = "All",
        subgroup = subgroup_summary$subgroup,
        n = subgroup_summary$n,
        events = subgroup_summary$events,
        event_rate = subgroup_summary$event_rate,
        stringsAsFactors = FALSE
      )
    }
  } else {
    subgroup_summary <- data.frame(
      subgroup = "All",
      n = nrow(data),
      events = sum(data[[event]], na.rm = TRUE),
      event_rate = mean(data[[event]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    arm_subgroup_summary <- data.frame(
      arm = arm_summary$arm,
      subgroup = "All",
      n = arm_summary$n,
      events = arm_summary$events,
      event_rate = arm_summary$event_rate,
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    overall = overall,
    by_arm = arm_summary,
    by_subgroup = subgroup_summary,
    by_arm_subgroup = arm_subgroup_summary,
    missingness = missingness
  )
  class(out) <- "cce_profile"
  out
}

#' @export
print.cce_profile <- function(x, ...) {
  cat("CCE dataset profile\n")
  print(x$overall)
  invisible(x)
}
