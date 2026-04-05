make_surv_formula <- function(time, event, rhs) {
  rhs_txt <- if (length(rhs) == 0L) "1" else paste(rhs, collapse = " + ")
  stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time, event, rhs_txt))
}

make_ps_formula <- function(covariates) {
  rhs_txt <- if (length(covariates) == 0L) "1" else paste(covariates, collapse = " + ")
  stats::as.formula(sprintf("arm_binary ~ %s", rhs_txt))
}

n_risk_at <- function(time, at) {
  vapply(at, function(tt) sum(time >= tt), numeric(1))
}

diagnostic_row <- function(method, subgroup, metric, value, threshold = NA_real_, status = "info") {
  data.frame(
    method = method,
    subgroup = subgroup,
    metric = metric,
    value = value,
    threshold = threshold,
    status = status,
    stringsAsFactors = FALSE
  )
}

fit_vs_subgroup <- function(
    data,
    arm,
    time,
    event,
    covariates,
    subgroup_name,
    tau,
    landmark_times,
    times,
    warning_max_weight,
    fail_max_weight,
    warning_smd,
    fail_smd,
    weight_cap) {
  arm_levels <- unique(as.character(data[[arm]]))
  if (all(c("SOC", "A") %in% arm_levels)) {
    arm_levels <- c("SOC", "A")
  } else {
    arm_levels <- sort(arm_levels)
  }
  data[[arm]] <- factor(as.character(data[[arm]]), levels = arm_levels)
  if (nlevels(data[[arm]]) != 2L || any(table(data[[arm]]) == 0L)) {
    stop(sprintf("Subgroup `%s` does not contain two treatment arms.", subgroup_name),
      call. = FALSE
    )
  }
  active_covariates <- covariates[vapply(data[covariates], function(x) {
    length(unique(x[!is.na(x)])) > 1L
  }, logical(1))]

  g_model <- survival::coxph(
    formula = make_surv_formula(time, event, c(arm, active_covariates)),
    data = data,
    ties = "efron",
    x = TRUE,
    model = TRUE
  )
  g_curves <- lapply(arm_levels, function(level) {
    surv <- standardized_survival_from_cox(g_model, data, arm = arm, arm_level = level, times = times)
    data.frame(
      mode = "vs",
      method = "gformula",
      subgroup = subgroup_name,
      arm = level,
      time = times,
      survival = surv,
      n_risk = n_risk_at(data[data[[arm]] == level, time], times),
      stringsAsFactors = FALSE
    )
  })

  data$arm_binary <- as.integer(data[[arm]] == arm_levels[2L])
  ps_model <- stats::glm(
    formula = make_ps_formula(active_covariates),
    family = stats::binomial(),
    data = data
  )
  ps <- stats::predict(ps_model, type = "response")
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  pa <- mean(data$arm_binary)
  weights <- ifelse(data$arm_binary == 1L, pa / ps, (1 - pa) / (1 - ps))
  weights <- pmin(weights, weight_cap)
  data$iptw_weight <- weights

  iptw_model <- survival::coxph(
    formula = make_surv_formula(time, event, arm),
    data = data,
    weights = iptw_weight,
    robust = TRUE,
    ties = "efron",
    x = TRUE,
    model = TRUE
  )
  iptw_curves <- lapply(arm_levels, function(level) {
    surv <- standardized_survival_from_cox(iptw_model, data, arm = arm, arm_level = level, times = times)
    data.frame(
      mode = "vs",
      method = "iptw",
      subgroup = subgroup_name,
      arm = level,
      time = times,
      survival = surv,
      n_risk = n_risk_at(data[data[[arm]] == level, time], times),
      stringsAsFactors = FALSE
    )
  })

  smd_before <- compute_smd(data, arm, active_covariates)
  smd_after <- compute_smd(data, arm, active_covariates, weights = weights)
  max_before <- if (length(smd_before) == 0L) 0 else max(abs(smd_before))
  max_after <- if (length(smd_after) == 0L) 0 else max(abs(smd_after))

  diagnostics <- do.call(rbind, list(
    diagnostic_row("iptw", subgroup_name, "max_weight", max(weights), fail_max_weight,
      if (max(weights) > fail_max_weight) "fail" else if (max(weights) > warning_max_weight) "warning" else "ok"
    ),
    diagnostic_row("iptw", subgroup_name, "ess_total", ess(weights), NA_real_, "info"),
    diagnostic_row("iptw", subgroup_name, "ess_arm0", ess(weights[data[[arm]] == arm_levels[1L]]), NA_real_, "info"),
    diagnostic_row("iptw", subgroup_name, "ess_arm1", ess(weights[data[[arm]] == arm_levels[2L]]), NA_real_, "info"),
    diagnostic_row("iptw", subgroup_name, "max_abs_smd_before", max_before, warning_smd,
      if (max_before > fail_smd) "fail" else if (max_before > warning_smd) "warning" else "ok"
    ),
    diagnostic_row("iptw", subgroup_name, "max_abs_smd_after", max_after, warning_smd,
      if (max_after > fail_smd) "fail" else if (max_after > warning_smd) "warning" else "ok"
    ),
    diagnostic_row("iptw", subgroup_name, "ps_min_arm0", min(ps[data[[arm]] == arm_levels[1L]]), NA_real_, "info"),
    diagnostic_row("iptw", subgroup_name, "ps_max_arm0", max(ps[data[[arm]] == arm_levels[1L]]), NA_real_, "info"),
    diagnostic_row("iptw", subgroup_name, "ps_min_arm1", min(ps[data[[arm]] == arm_levels[2L]]), NA_real_, "info"),
    diagnostic_row("iptw", subgroup_name, "ps_max_arm1", max(ps[data[[arm]] == arm_levels[2L]]), NA_real_, "info"),
    diagnostic_row("iptw", subgroup_name, "n_arm0", sum(data[[arm]] == arm_levels[1L]), NA_real_, "info"),
    diagnostic_row("iptw", subgroup_name, "n_arm1", sum(data[[arm]] == arm_levels[2L]), NA_real_, "info"),
    diagnostic_row("iptw", subgroup_name, "events_arm0", sum(data[[event]][data[[arm]] == arm_levels[1L]]), NA_real_, "info"),
    diagnostic_row("iptw", subgroup_name, "events_arm1", sum(data[[event]][data[[arm]] == arm_levels[2L]]), NA_real_, "info")
  ))

  curves <- do.call(rbind, c(g_curves, iptw_curves))
  effects <- effect_rows_from_curves(curves, subgroup = subgroup_name, tau = tau, landmark_times = landmark_times)
  list(curves = curves, effects = effects, diagnostics = diagnostics)
}

fit_vs_once <- function(
    data,
    arm,
    time,
    event,
    covariates,
    subgroup,
    tau,
    landmark_times,
    times,
    warning_max_weight,
    fail_max_weight,
    warning_smd,
    fail_smd,
    weight_cap) {
  group_list <- list(All = data)
  if (!is.null(subgroup)) {
    split_groups <- split(data, as.character(data[[subgroup]]))
    group_list <- c(group_list, split_groups)
  }

  pieces <- lapply(names(group_list), function(name) {
    fit_vs_subgroup(
      data = group_list[[name]],
      arm = arm,
      time = time,
      event = event,
      covariates = covariates,
      subgroup_name = name,
      tau = tau,
      landmark_times = landmark_times,
      times = times,
      warning_max_weight = warning_max_weight,
      fail_max_weight = fail_max_weight,
      warning_smd = warning_smd,
      fail_smd = fail_smd,
      weight_cap = weight_cap
    )
  })

  curves <- do.call(rbind, lapply(pieces, `[[`, "curves"))
  effects <- do.call(rbind, lapply(pieces, `[[`, "effects"))
  diagnostics <- do.call(rbind, lapply(pieces, `[[`, "diagnostics"))

  discordance <- FALSE
  all_effects <- effects[effects$subgroup == "All", , drop = FALSE]
  if (length(unique(all_effects$method)) == 2L) {
    delta_by_method <- stats::aggregate(delta_rmst ~ method, data = all_effects[all_effects$landmark_time == max(landmark_times), , drop = FALSE], FUN = mean)
    if (nrow(delta_by_method) == 2L) {
      discordance <- sign(delta_by_method$delta_rmst[1L]) != sign(delta_by_method$delta_rmst[2L])
    }
  }
  fail_flags <- list(
    max_weight_fail = any(diagnostics$metric == "max_weight" & diagnostics$status == "fail"),
    balance_fail = any(diagnostics$metric == "max_abs_smd_after" & diagnostics$status == "fail"),
    method_discordance = discordance
  )
  warnings <- unique(c(
    if (any(diagnostics$metric == "max_weight" & diagnostics$status == "warning")) "High IPTW weights detected." else NULL,
    if (any(diagnostics$metric == "max_abs_smd_after" & diagnostics$status == "warning")) "Residual covariate imbalance detected." else NULL,
    if (discordance) "g-formula and IPTW disagree on the RMST direction." else NULL
  ))
  list(curves = curves, effects = effects, diagnostics = diagnostics, fail_flags = fail_flags, warnings = warnings)
}

bootstrap_vs <- function(
    data,
    arm,
    time,
    event,
    covariates,
    subgroup,
    tau,
    landmark_times,
    times,
    bootstrap,
    seed,
    warning_max_weight,
    fail_max_weight,
    warning_smd,
    fail_smd,
    weight_cap) {
  set.seed(seed)
  curve_parts <- vector("list", bootstrap)
  effect_parts <- vector("list", bootstrap)
  for (b in seq_len(bootstrap)) {
    idx <- sample.int(nrow(data), replace = TRUE)
    res <- fit_vs_once(
      data = data[idx, , drop = FALSE],
      arm = arm,
      time = time,
      event = event,
      covariates = covariates,
      subgroup = subgroup,
      tau = tau,
      landmark_times = landmark_times,
      times = times,
      warning_max_weight = warning_max_weight,
      fail_max_weight = fail_max_weight,
      warning_smd = warning_smd,
      fail_smd = fail_smd,
      weight_cap = weight_cap
    )
    curve_parts[[b]] <- cbind(res$curves, .replicate = b)
    effect_parts[[b]] <- cbind(res$effects, .replicate = b)
  }
  list(
    curves = do.call(rbind, curve_parts),
    effects = do.call(rbind, effect_parts)
  )
}

#' Estimate a counterfactual VS comparison
#'
#' `fit_cce_vs()` fits two complementary estimators for a binary treatment
#' comparison: a Cox-model-based g-formula standardization and an inverse
#' probability of treatment weighting analysis. The function returns tidy
#' curves, effects, diagnostics, and machine-readable metadata.
#'
#' @param data Analysis-ready data frame.
#' @param arm,time,event Column names identifying treatment assignment,
#'   follow-up time, and event indicator.
#' @param covariates Character vector of baseline adjustment covariates.
#' @param subgroup Optional subgroup column name. When supplied, the result
#'   includes overall and subgroup-specific summaries.
#' @param tau RMST truncation horizon. Defaults to the 90th percentile of
#'   observed follow-up.
#' @param landmark_times Survival-difference time points.
#' @param n_grid Number of time points used to summarize curves.
#' @param bootstrap Number of bootstrap resamples used to derive intervals.
#' @param seed Random seed used for bootstrap resampling.
#' @param weight_cap Hard upper cap applied to stabilized IPTW weights.
#' @param warning_max_weight,fail_max_weight Diagnostic thresholds for the
#'   largest stabilized weight.
#' @param warning_smd,fail_smd Diagnostic thresholds for absolute SMD.
#'
#' @return An object of class `cce_vs_result`.
#' @export
fit_cce_vs <- function(
    data,
    arm = "arm",
    time = "time",
    event = "event",
    covariates,
    subgroup = NULL,
    tau = NULL,
    landmark_times = NULL,
    n_grid = 100L,
    bootstrap = 0L,
    seed = 1L,
    weight_cap = 50,
    warning_max_weight = 10,
    fail_max_weight = 50,
    warning_smd = 0.10,
    fail_smd = 0.20) {
  assert_data_frame(data, "data")
  assert_scalar_character(arm, "arm")
  assert_scalar_character(time, "time")
  assert_scalar_character(event, "event")
  assert_character_vector(covariates, "covariates")
  if (!is.null(subgroup)) {
    assert_scalar_character(subgroup, "subgroup")
  }
  assert_named_columns(data, c(arm, time, event, covariates), "data")
  if (!is.null(subgroup)) {
    assert_named_columns(data, subgroup, "data")
  }

  keep_cols <- c(arm, time, event, covariates, subgroup)
  keep_cols <- keep_cols[!is.na(keep_cols) & nzchar(keep_cols)]
  data <- data[stats::complete.cases(data[, keep_cols, drop = FALSE]), , drop = FALSE]
  if (!is.factor(data[[arm]])) {
    data[[arm]] <- factor(data[[arm]])
  }
  if (nlevels(data[[arm]]) != 2L) {
    stop("`arm` must contain exactly two levels.", call. = FALSE)
  }
  if (is.null(tau)) {
    tau <- unname(stats::quantile(data[[time]], probs = 0.9, na.rm = TRUE))
  }
  assert_scalar_numeric(tau, "tau", lower = 1)
  if (is.null(landmark_times)) {
    landmark_times <- unique(sort(c(round(tau / 2), round(tau))))
  }
  landmark_times <- sort(unique(as.numeric(landmark_times)))
  times <- default_time_grid(data[[time]], tau = tau, n_grid = n_grid)

  base_res <- fit_vs_once(
    data = data,
    arm = arm,
    time = time,
    event = event,
    covariates = covariates,
    subgroup = subgroup,
    tau = tau,
    landmark_times = landmark_times,
    times = times,
    warning_max_weight = warning_max_weight,
    fail_max_weight = fail_max_weight,
    warning_smd = warning_smd,
    fail_smd = fail_smd,
    weight_cap = weight_cap
  )

  curves <- base_res$curves
  curves$lower_ci <- NA_real_
  curves$upper_ci <- NA_real_
  effects <- base_res$effects
  effects$delta_rmst_lower_ci <- NA_real_
  effects$delta_rmst_upper_ci <- NA_real_
  effects$delta_survival_lower_ci <- NA_real_
  effects$delta_survival_upper_ci <- NA_real_

  if (bootstrap > 0L) {
    boot <- bootstrap_vs(
      data = data,
      arm = arm,
      time = time,
      event = event,
      covariates = covariates,
      subgroup = subgroup,
      tau = tau,
      landmark_times = landmark_times,
      times = times,
      bootstrap = bootstrap,
      seed = seed,
      warning_max_weight = warning_max_weight,
      fail_max_weight = fail_max_weight,
      warning_smd = warning_smd,
      fail_smd = fail_smd,
      weight_cap = weight_cap
    )
    curve_intervals <- bootstrap_quantiles(
      boot$curves,
      key_cols = c("mode", "method", "subgroup", "arm", "time"),
      value_cols = "survival"
    )
    effect_intervals <- bootstrap_quantiles(
      boot$effects,
      key_cols = c("mode", "method", "subgroup", "tau", "landmark_time"),
      value_cols = c("delta_rmst", "delta_survival")
    )
    curves <- merge_curve_intervals(curves, curve_intervals)
    effects <- merge_effect_intervals(effects, effect_intervals)
  }

  result <- list(
    curves = curves,
    effects = effects,
    diagnostics = base_res$diagnostics,
    warnings = base_res$warnings,
    fail_flags = base_res$fail_flags,
    label = label_from_fail_flags(base_res$fail_flags),
    meta = list(
      mode = "vs",
      tau = tau,
      landmark_times = landmark_times,
      bootstrap = bootstrap,
      generated_at = as.character(Sys.time()),
      run_id = run_stamp()
    )
  )
  class(result) <- "cce_vs_result"
  result
}

#' @export
summary.cce_vs_result <- function(object, ...) {
  cat("CCE VS result\n")
  cat("Label:", object$label, "\n")
  cat("Warnings:", if (length(object$warnings) == 0L) "none" else paste(object$warnings, collapse = " | "), "\n")
  print(utils::head(object$effects, n = min(6L, nrow(object$effects))))
  invisible(object)
}

#' @export
plot.cce_vs_result <- function(x, method = "gformula", subgroup = "All", ...) {
  curves <- x$curves[x$curves$method == method & x$curves$subgroup == subgroup, , drop = FALSE]
  if (nrow(curves) == 0L) {
    stop("No matching curves found for the requested method/subgroup.", call. = FALSE)
  }
  arms <- unique(curves$arm)
  cols <- c("#1b9e77", "#d95f02")
  plot(NA,
    xlim = range(curves$time), ylim = c(0, 1), xlab = "Time", ylab = "Survival",
    main = sprintf("CCE VS curves: %s / %s", method, subgroup), ...
  )
  for (i in seq_along(arms)) {
    df <- curves[curves$arm == arms[i], , drop = FALSE]
    graphics::lines(df$time, df$survival, lwd = 2, col = cols[i])
  }
  graphics::legend("topright", legend = arms, col = cols[seq_along(arms)], lwd = 2, bty = "n")
}
