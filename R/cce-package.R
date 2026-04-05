#' cce: Counterfactual Comparator Engine
#'
#' The package ships a compact but production-oriented toolkit for
#' counterfactual comparator analyses in oncology-style survival workflows.
#' It covers normalized input validation, analysis cohort assembly,
#' standardized survival estimation, SOC-only projection, diagnostics, and
#' machine-readable exports.
#'
#' @keywords internal
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("iptw_weight", "km_weight"))
}
