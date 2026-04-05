#' Create a Counterfactual Comparator Engine specification
#'
#' `cce_spec()` defines the minimal schema needed to convert normalized input
#' tables into an analysis-ready dataset. The returned object is intentionally
#' simple and can be stored as YAML for reproducible runs.
#'
#' @param covariates Character vector of baseline covariate column names.
#' @param subgroup_biomarker Optional biomarker name used to derive a subgroup
#'   column from the `biomarkers` table.
#' @param endpoint Endpoint name to retain from the `outcomes` table.
#' @param id_col Patient identifier column name.
#' @param index_date_col Baseline index-date column name.
#' @param regimen_col Treatment regimen column name.
#' @param treatment_start_col Index treatment start-date column name.
#' @param index_flag_col Logical column marking the index treatment row.
#' @param endpoint_col Outcome endpoint-name column.
#' @param time_col Follow-up time column.
#' @param event_col Event indicator column.
#' @param follow_up_col Last follow-up date column.
#' @param biomarker_name_col Biomarker-name column.
#' @param biomarker_value_col Biomarker-value column.
#' @param biomarker_baseline_flag_col Logical baseline-biomarker flag column.
#' @param arm_map Named character vector mapping raw regimen labels to output
#'   labels. The output labels must be exactly `c("SOC", "A")`.
#' @param missing_strategy Missing-data rule. Only `"complete_case"` is
#'   implemented in v0.1.0.
#' @param time_zero_tolerance_days Allowed difference between `index_date` and
#'   the index treatment start date.
#'
#' @return A `cce_spec` object.
#' @export
cce_spec <- function(
    covariates,
    subgroup_biomarker = NULL,
    endpoint = "os",
    id_col = "patient_id",
    index_date_col = "index_date",
    regimen_col = "regimen_name",
    treatment_start_col = "start_date",
    index_flag_col = "is_index_treatment",
    endpoint_col = "endpoint",
    time_col = "time",
    event_col = "event",
    follow_up_col = "last_follow_up_date",
    biomarker_name_col = "biomarker_name",
    biomarker_value_col = "biomarker_value",
    biomarker_baseline_flag_col = "is_baseline",
    arm_map = c(SOC = "SOC", A = "A"),
    missing_strategy = "complete_case",
    time_zero_tolerance_days = 0L) {
  assert_character_vector(covariates, "covariates")
  if (!is.null(subgroup_biomarker)) {
    assert_scalar_character(subgroup_biomarker, "subgroup_biomarker")
  }
  assert_scalar_character(endpoint, "endpoint")
  for (nm in c(
    "id_col", "index_date_col", "regimen_col", "treatment_start_col",
    "index_flag_col", "endpoint_col", "time_col", "event_col",
    "follow_up_col", "biomarker_name_col", "biomarker_value_col",
    "biomarker_baseline_flag_col", "missing_strategy"
  )) {
    assert_scalar_character(get(nm), nm)
  }
  if (!is.character(arm_map) || is.null(names(arm_map)) || anyNA(names(arm_map)) || anyNA(arm_map)) {
    stop("`arm_map` must be a named character vector.", call. = FALSE)
  }
  if (!identical(sort(unname(arm_map)), c("A", "SOC"))) {
    stop("`arm_map` output labels must be exactly c('SOC', 'A').", call. = FALSE)
  }
  assert_scalar_numeric(time_zero_tolerance_days, "time_zero_tolerance_days", lower = 0)
  structure(list(
    covariates = covariates,
    subgroup_biomarker = subgroup_biomarker,
    endpoint = endpoint,
    id_col = id_col,
    index_date_col = index_date_col,
    regimen_col = regimen_col,
    treatment_start_col = treatment_start_col,
    index_flag_col = index_flag_col,
    endpoint_col = endpoint_col,
    time_col = time_col,
    event_col = event_col,
    follow_up_col = follow_up_col,
    biomarker_name_col = biomarker_name_col,
    biomarker_value_col = biomarker_value_col,
    biomarker_baseline_flag_col = biomarker_baseline_flag_col,
    arm_map = arm_map,
    missing_strategy = missing_strategy,
    time_zero_tolerance_days = as.integer(time_zero_tolerance_days)
  ), class = "cce_spec")
}

#' Write a CCE specification to YAML
#'
#' @param spec A `cce_spec` object.
#' @param path File path ending in `.yml` or `.yaml`.
#'
#' @return Invisibly returns `path`.
#' @export
write_cce_spec <- function(spec, path) {
  if (!inherits(spec, "cce_spec")) {
    stop("`spec` must inherit from `cce_spec`.", call. = FALSE)
  }
  assert_scalar_character(path, "path")
  payload <- unclass(spec)
  payload$arm_map <- as.list(stats::setNames(unname(spec$arm_map), names(spec$arm_map)))
  yaml::write_yaml(payload, file = path)
  invisible(path)
}

#' Read a CCE specification from YAML
#'
#' @param path YAML file created by [write_cce_spec()].
#'
#' @return A `cce_spec` object.
#' @export
read_cce_spec <- function(path) {
  assert_scalar_character(path, "path")
  spec <- yaml::read_yaml(path)
  if (is.list(spec$arm_map) && !is.null(names(spec$arm_map))) {
    spec$arm_map <- stats::setNames(unlist(spec$arm_map, use.names = FALSE), names(spec$arm_map))
  } else if (is.character(spec$arm_map) && is.null(names(spec$arm_map)) && length(spec$arm_map) == 2L) {
    spec$arm_map <- stats::setNames(spec$arm_map, spec$arm_map)
  }
  do.call(cce_spec, spec)
}
