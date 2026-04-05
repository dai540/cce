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
  payload <- c(
    list(
      label = x$label,
      warnings = x$warnings,
      fail_flags = x$fail_flags
    ),
    x$meta
  )
  json_txt <- jsonlite::toJSON(
    payload,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null",
    dataframe = "rows",
    na = "null"
  )
  writeLines(json_txt, con = file.path(path, "results.json"), useBytes = TRUE)
  invisible(path)
}
