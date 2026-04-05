#' Extract curves from a CCE result
#'
#' @param x A `cce_vs_result` or `cce_soc_result`.
#'
#' @return A data frame.
#' @export
as_curves_df <- function(x) {
  x$curves
}

#' Extract effect summaries from a CCE result
#'
#' @param x A `cce_vs_result` or `cce_soc_result`.
#'
#' @return A data frame.
#' @export
as_effects_df <- function(x) {
  x$effects
}

#' Extract diagnostics from a CCE result
#'
#' @param x A `cce_vs_result` or `cce_soc_result`.
#'
#' @return A data frame.
#' @export
as_diagnostics_df <- function(x) {
  x$diagnostics
}

#' Write CCE result files to disk
#'
#' @param x A `cce_vs_result` or `cce_soc_result`.
#' @param path Output directory.
#'
#' @return Invisibly returns `path`.
#' @export
write_cce_results <- function(x, path = "outputs") {
  assert_scalar_character(path, "path")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(as_curves_df(x), file.path(path, "curves.csv"), row.names = FALSE)
  utils::write.csv(as_effects_df(x), file.path(path, "effects.csv"), row.names = FALSE)
  utils::write.csv(as_diagnostics_df(x), file.path(path, "diagnostics.csv"), row.names = FALSE)
  payload <- list(
    run_id = x$meta$run_id,
    mode = x$meta$mode,
    generated_at = x$meta$generated_at,
    tau = x$meta$tau,
    landmark_times = x$meta$landmark_times,
    label = x$label,
    warnings = x$warnings,
    fail_flags = x$fail_flags,
    input_hash = if (!is.null(x$meta$input_hash)) x$meta$input_hash else NA_character_,
    config_hash = if (!is.null(x$meta$config_hash)) x$meta$config_hash else NA_character_,
    code_hash = if (!is.null(x$meta$code_hash)) x$meta$code_hash else NA_character_
  )
  jsonlite::write_json(payload, path = file.path(path, "results.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(path)
}
