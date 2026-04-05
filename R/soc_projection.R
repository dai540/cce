km_survival <- function(data, time, event, times) {
  fit <- survival::survfit(make_surv_formula(time, event, "1"), data = data)
  evaluate_step(fit$time, fit$surv, times, initial = 1)
}

#' Reverse-solve the hazard ratio needed to hit an RMST target
#'
#' @param soc_time Time grid for the SOC survival curve.
#' @param soc_survival SOC survival values at `soc_time`.
#' @param tau Truncation horizon used for RMST.
#' @param target_delta_rmst Target RMST gain relative to SOC.
#' @param lower,upper Search interval for the hazard ratio.
#'
#' @return A scalar hazard ratio or `NA_real_` when no root is bracketed.
#' @export
required_hr <- function(soc_time, soc_survival, tau, target_delta_rmst, lower = 0.05, upper = 2.5) {
  assert_scalar_numeric(tau, "tau", lower = 1)
  assert_scalar_numeric(target_delta_rmst, "target_delta_rmst")
  fn <- function(theta) {
    rmst_from_curve(soc_time, soc_survival^theta, tau) -
      rmst_from_curve(soc_time, soc_survival, tau) -
      target_delta_rmst
  }
  if (identical(target_delta_rmst, 0)) {
    return(1)
  }
  f_lower <- fn(lower)
  f_upper <- fn(upper)
  if (!is.finite(f_lower) || !is.finite(f_upper) || f_lower * f_upper > 0) {
    return(NA_real_)
  }
  stats::uniroot(fn, interval = c(lower, upper))$root
}

#' Estimate an assumption-based probability-of-success proxy
#'
#' @param soc_time Time grid for the SOC curve.
#' @param soc_survival SOC survival values at `soc_time`.
#' @param tau Truncation horizon.
#' @param target_delta_rmst Target RMST gain.
#' @param mean_log_hr,sd_log_hr Mean and standard deviation of the log-HR prior.
#' @param draws Number of Monte Carlo draws.
#' @param seed Optional seed.
#'
#' @return A scalar probability.
#' @export
estimate_pos_proxy <- function(
    soc_time,
    soc_survival,
    tau,
    target_delta_rmst,
    mean_log_hr,
    sd_log_hr,
    draws = 2000L,
    seed = NULL) {
  assert_scalar_numeric(mean_log_hr, "mean_log_hr")
  assert_scalar_numeric(sd_log_hr, "sd_log_hr", lower = 0)
  assert_scalar_numeric(draws, "draws", lower = 100)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  theta <- exp(stats::rnorm(draws, mean = mean_log_hr, sd = sd_log_hr))
  deltas <- vapply(theta, function(tt) {
    rmst_from_curve(soc_time, soc_survival^tt, tau) - rmst_from_curve(soc_time, soc_survival, tau)
  }, numeric(1))
  mean(deltas >= target_delta_rmst)
}

fit_soc_subgroup <- function(
    data,
    time,
    event,
    subgroup_name,
    tau,
    landmark_times,
    times,
    hr_scenarios,
    target_delta_rmst,
    prior_mean_log_hr,
    prior_sd_log_hr,
    prior_draws,
    seed) {
  soc_surv <- km_survival(data, time = time, event = event, times = times)
  curves <- vector("list", length(hr_scenarios) + 1L)
  curves[[1L]] <- data.frame(
    mode = "soc_only",
    method = "projection_ph",
    subgroup = subgroup_name,
    arm = "SOC",
    scenario_hr = NA_real_,
    time = times,
    survival = soc_surv,
    n_risk = n_risk_at(data[[time]], times),
    stringsAsFactors = FALSE
  )
  effect_rows <- vector("list", length(hr_scenarios))
  req_hr <- if (!is.null(target_delta_rmst)) {
    required_hr(times, soc_surv, tau = tau, target_delta_rmst = target_delta_rmst)
  } else {
    NA_real_
  }
  pos_proxy <- if (!is.null(target_delta_rmst) && !is.null(prior_mean_log_hr) && !is.null(prior_sd_log_hr)) {
    estimate_pos_proxy(
      soc_time = times,
      soc_survival = soc_surv,
      tau = tau,
      target_delta_rmst = target_delta_rmst,
      mean_log_hr = prior_mean_log_hr,
      sd_log_hr = prior_sd_log_hr,
      draws = prior_draws,
      seed = seed
    )
  } else {
    NA_real_
  }

  rmst_soc <- rmst_from_curve(times, soc_surv, tau)
  for (i in seq_along(hr_scenarios)) {
    theta <- hr_scenarios[i]
    proj <- soc_surv^theta
    curves[[i + 1L]] <- data.frame(
      mode = "soc_only",
      method = "projection_ph",
      subgroup = subgroup_name,
      arm = sprintf("SOC+X(theta=%.2f)", theta),
      scenario_hr = theta,
      time = times,
      survival = proj,
      n_risk = n_risk_at(data[[time]], times),
      stringsAsFactors = FALSE
    )
    rows <- vector("list", length(landmark_times))
    rmst_proj <- rmst_from_curve(times, proj, tau)
    for (j in seq_along(landmark_times)) {
      tt <- landmark_times[j]
      s_soc <- evaluate_step(times, soc_surv, tt, initial = 1)
      s_proj <- evaluate_step(times, proj, tt, initial = 1)
      rows[[j]] <- data.frame(
        mode = "soc_only",
        method = "projection_ph",
        subgroup = subgroup_name,
        scenario_hr = theta,
        tau = tau,
        rmst_arm0 = rmst_soc,
        rmst_arm1 = rmst_proj,
        delta_rmst = rmst_proj - rmst_soc,
        landmark_time = tt,
        survival_arm0 = s_soc,
        survival_arm1 = s_proj,
        delta_survival = s_proj - s_soc,
        required_hr = req_hr,
        pos_proxy = pos_proxy,
        stringsAsFactors = FALSE
      )
    }
    effect_rows[[i]] <- do.call(rbind, rows)
  }

  list(
    curves = do.call(rbind, curves),
    effects = do.call(rbind, effect_rows),
    diagnostics = diagnostic_row("projection_ph", subgroup_name, "n_soc", nrow(data), NA_real_, "info")
  )
}

fit_soc_once <- function(
    data,
    time,
    event,
    subgroup,
    tau,
    landmark_times,
    times,
    hr_scenarios,
    target_delta_rmst,
    prior_mean_log_hr,
    prior_sd_log_hr,
    prior_draws,
    seed) {
  group_list <- list(All = data)
  if (!is.null(subgroup)) {
    group_list <- c(group_list, split(data, as.character(data[[subgroup]])))
  }
  pieces <- lapply(names(group_list), function(name) {
    fit_soc_subgroup(
      data = group_list[[name]],
      time = time,
      event = event,
      subgroup_name = name,
      tau = tau,
      landmark_times = landmark_times,
      times = times,
      hr_scenarios = hr_scenarios,
      target_delta_rmst = target_delta_rmst,
      prior_mean_log_hr = prior_mean_log_hr,
      prior_sd_log_hr = prior_sd_log_hr,
      prior_draws = prior_draws,
      seed = seed
    )
  })
  list(
    curves = do.call(rbind, lapply(pieces, `[[`, "curves")),
    effects = do.call(rbind, lapply(pieces, `[[`, "effects")),
    diagnostics = do.call(rbind, lapply(pieces, `[[`, "diagnostics"))
  )
}

bootstrap_soc <- function(
    data,
    time,
    event,
    subgroup,
    tau,
    landmark_times,
    times,
    hr_scenarios,
    target_delta_rmst,
    prior_mean_log_hr,
    prior_sd_log_hr,
    prior_draws,
    bootstrap,
    seed) {
  set.seed(seed)
  curve_parts <- vector("list", bootstrap)
  effect_parts <- vector("list", bootstrap)
  for (b in seq_len(bootstrap)) {
    idx <- sample.int(nrow(data), replace = TRUE)
    res <- fit_soc_once(
      data = data[idx, , drop = FALSE],
      time = time,
      event = event,
      subgroup = subgroup,
      tau = tau,
      landmark_times = landmark_times,
      times = times,
      hr_scenarios = hr_scenarios,
      target_delta_rmst = target_delta_rmst,
      prior_mean_log_hr = prior_mean_log_hr,
      prior_sd_log_hr = prior_sd_log_hr,
      prior_draws = prior_draws,
      seed = seed + b
    )
    curve_parts[[b]] <- cbind(res$curves, .replicate = b)
    effect_parts[[b]] <- cbind(res$effects, .replicate = b)
  }
  list(
    curves = do.call(rbind, curve_parts),
    effects = do.call(rbind, effect_parts)
  )
}

#' Project SOC-only curves under proportional hazards scenarios
#'
#' @param data Analysis-ready data frame.
#' @param time,event Column names identifying follow-up and event status.
#' @param subgroup Optional subgroup column name.
#' @param arm Optional treatment column. When supplied, only rows matching
#'   `soc_level` are retained before projection.
#' @param soc_level Reference level used when `arm` is supplied.
#' @param tau RMST truncation horizon.
#' @param landmark_times Survival-difference time points.
#' @param hr_scenarios Numeric vector of proportional-hazard multipliers.
#' @param target_delta_rmst Optional target RMST gain used to reverse-solve
#'   `required_hr`.
#' @param prior_mean_log_hr,prior_sd_log_hr Optional normal prior parameters
#'   for the log hazard ratio. When both are supplied alongside
#'   `target_delta_rmst`, a PoS proxy is reported.
#' @param prior_draws Number of Monte Carlo draws for the PoS proxy.
#' @param bootstrap Number of bootstrap resamples used to derive intervals.
#' @param seed Random seed.
#' @param n_grid Number of time points used to summarize curves.
#'
#' @return An object of class `cce_soc_result`.
#' @export
project_soc_only <- function(
    data,
    time = "time",
    event = "event",
    subgroup = NULL,
    arm = NULL,
    soc_level = "SOC",
    tau = NULL,
    landmark_times = NULL,
    hr_scenarios = c(0.70, 0.85, 1.00),
    target_delta_rmst = NULL,
    prior_mean_log_hr = NULL,
    prior_sd_log_hr = NULL,
    prior_draws = 2000L,
    bootstrap = 0L,
    seed = 1L,
    n_grid = 100L) {
  source_attrs <- capture_dataset_attributes(data)
  assert_data_frame(data, "data")
  assert_scalar_character(time, "time")
  assert_scalar_character(event, "event")
  assert_named_columns(data, c(time, event), "data")
  if (!is.null(subgroup)) {
    assert_scalar_character(subgroup, "subgroup")
    assert_named_columns(data, subgroup, "data")
  }
  if (!is.null(arm)) {
    assert_scalar_character(arm, "arm")
    assert_named_columns(data, arm, "data")
    data <- data[as.character(data[[arm]]) == soc_level, , drop = FALSE]
  }
  data <- data[stats::complete.cases(data[, c(time, event, subgroup), drop = FALSE]), , drop = FALSE]
  if (nrow(data) == 0L) {
    stop("No SOC rows are available for projection.", call. = FALSE)
  }
  if (is.null(tau)) {
    tau <- unname(stats::quantile(data[[time]], probs = 0.9, na.rm = TRUE))
  }
  if (is.null(landmark_times)) {
    landmark_times <- unique(sort(c(round(tau / 2), round(tau))))
  }
  hr_scenarios <- sort(unique(as.numeric(hr_scenarios)))
  times <- default_time_grid(data[[time]], tau = tau, n_grid = n_grid)
  fit_profile <- profile_cce_dataset(
    data = data,
    arm = if (!is.null(arm)) arm else NULL,
    time = time,
    event = event,
    subgroup = subgroup
  )

  base_res <- fit_soc_once(
    data = data,
    time = time,
    event = event,
    subgroup = subgroup,
    tau = tau,
    landmark_times = landmark_times,
    times = times,
    hr_scenarios = hr_scenarios,
    target_delta_rmst = target_delta_rmst,
    prior_mean_log_hr = prior_mean_log_hr,
    prior_sd_log_hr = prior_sd_log_hr,
    prior_draws = prior_draws,
    seed = seed
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
    boot <- bootstrap_soc(
      data = data,
      time = time,
      event = event,
      subgroup = subgroup,
      tau = tau,
      landmark_times = landmark_times,
      times = times,
      hr_scenarios = hr_scenarios,
      target_delta_rmst = target_delta_rmst,
      prior_mean_log_hr = prior_mean_log_hr,
      prior_sd_log_hr = prior_sd_log_hr,
      prior_draws = prior_draws,
      bootstrap = bootstrap,
      seed = seed
    )
    curve_intervals <- bootstrap_quantiles(
      boot$curves,
      key_cols = c("mode", "method", "subgroup", "arm", "scenario_hr", "time"),
      value_cols = "survival"
    )
    effect_intervals <- bootstrap_quantiles(
      boot$effects,
      key_cols = c("mode", "method", "subgroup", "scenario_hr", "tau", "landmark_time"),
      value_cols = c("delta_rmst", "delta_survival")
    )
    curves <- merge_interval_columns(curves, curve_intervals,
      c("mode", "method", "subgroup", "arm", "scenario_hr", "time"), "survival"
    )
    names(curves)[names(curves) == "survival_lower_ci"] <- "lower_ci"
    names(curves)[names(curves) == "survival_upper_ci"] <- "upper_ci"
    effects <- merge_interval_columns(effects, effect_intervals,
      c("mode", "method", "subgroup", "scenario_hr", "tau", "landmark_time"), "delta_rmst"
    )
    effects <- merge_interval_columns(effects, effect_intervals,
      c("mode", "method", "subgroup", "scenario_hr", "tau", "landmark_time"), "delta_survival"
    )
  }

  result <- list(
    curves = curves,
    effects = effects,
    diagnostics = base_res$diagnostics,
    warnings = "Projection (assumption-based)",
    fail_flags = list(reference_only = TRUE),
    label = "Projection (assumption-based)",
    meta = list(
      mode = "soc_only",
      tau = tau,
      landmark_times = landmark_times,
      bootstrap = bootstrap,
      estimators = "projection_ph",
      columns = list(
        arm = arm,
        time = time,
        event = event,
        subgroup = subgroup
      ),
      soc_level = soc_level,
      subgroup_levels = if (!is.null(subgroup)) sort(unique(as.character(data[[subgroup]]))) else "All",
      hr_scenarios = hr_scenarios,
      target_delta_rmst = target_delta_rmst,
      prior = list(
        mean_log_hr = prior_mean_log_hr,
        sd_log_hr = prior_sd_log_hr,
        draws = prior_draws
      ),
      spec = if (!is.null(source_attrs$spec)) unclass(source_attrs$spec) else NULL,
      exclusions = source_attrs$exclusions,
      validation_report = source_attrs$validation_report,
      profile = profile_list_for_meta(fit_profile),
      generated_at = as.character(Sys.time()),
      run_id = run_stamp()
    )
  )
  class(result) <- "cce_soc_result"
  result
}

#' @export
summary.cce_soc_result <- function(object, ...) {
  cat("CCE SOC-only projection\n")
  cat("Label:", object$label, "\n")
  print(utils::head(object$effects, n = min(6L, nrow(object$effects))))
  invisible(object)
}

#' @export
plot.cce_soc_result <- function(x, subgroup = "All", ...) {
  curves <- x$curves[x$curves$subgroup == subgroup, , drop = FALSE]
  if (nrow(curves) == 0L) {
    stop("No matching curves found for the requested subgroup.", call. = FALSE)
  }
  arms <- unique(curves$arm)
  cols <- grDevices::hcl.colors(length(arms), "Dark 3")
  plot(NA,
    xlim = range(curves$time), ylim = c(0, 1), xlab = "Time", ylab = "Survival",
    main = sprintf("CCE SOC-only projection: %s", subgroup), ...
  )
  for (i in seq_along(arms)) {
    df <- curves[curves$arm == arms[i], , drop = FALSE]
    graphics::lines(df$time, df$survival, lwd = 2, col = cols[i])
  }
  graphics::legend("topright", legend = arms, col = cols, lwd = 2, bty = "n")
}
