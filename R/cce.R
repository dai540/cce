#' cce: Minimal Counterfactual Comparator Workflows
#'
#' `cce` provides a compact set of survival workflows for observed
#' two-arm comparisons and SOC-only projection scenarios.
#'
#' @keywords internal
"_PACKAGE"

#' Generate demo survival data
#'
#' Generates a small synthetic survival dataset for examples and tests.
#'
#' @param n Number of rows.
#' @param seed Random seed.
#'
#' @return A data frame with columns `arm`, `time`, `event`, `subgroup`,
#'   `age`, and `ps`.
#' @export
cce_demo_data <- function(n = 200, seed = 1) {
  stopifnot(length(n) == 1L, is.numeric(n), n >= 20)
  stopifnot(length(seed) == 1L, is.numeric(seed))
  set.seed(seed)

  arm <- sample(c("SOC", "A"), size = n, replace = TRUE)
  subgroup <- sample(c("All-comers", "Biomarker-high", "Biomarker-low"), size = n, replace = TRUE)
  age <- round(stats::rnorm(n, mean = 64, sd = 9))
  ps <- sample(0:2, size = n, replace = TRUE, prob = c(0.35, 0.45, 0.20))

  hazard <- ifelse(arm == "SOC", 1 / 280, 1 / 360)
  subgroup_multiplier <- ifelse(subgroup == "Biomarker-high", 0.8, ifelse(subgroup == "Biomarker-low", 1.15, 1.0))
  event_time <- stats::rexp(n, rate = hazard * subgroup_multiplier)
  censor_time <- stats::runif(n, min = 120, max = 540)

  time <- pmin(event_time, censor_time)
  event <- as.integer(event_time <= censor_time)

  data.frame(
    arm = arm,
    time = round(time),
    event = event,
    subgroup = subgroup,
    age = age,
    ps = ps,
    stringsAsFactors = FALSE
  )
}

#' Fit an observed two-arm survival comparison
#'
#' Fits a simple observed survival comparison and returns curves, effect tables,
#' diagnostics, and metadata on a shared output contract.
#'
#' @param data A data frame.
#' @param arm Arm column name.
#' @param time Follow-up time column name.
#' @param event Event indicator column name.
#' @param subgroup Optional subgroup column name.
#' @param tau Horizon for RMST calculation.
#' @param times Time points used for survival differences.
#'
#' @return A `cce_fit` object.
#' @export
fit_cce_vs <- function(data,
                       arm = "arm",
                       time = "time",
                       event = "event",
                       subgroup = NULL,
                       tau,
                       times = tau) {
  .validate_common_inputs(data, arm, time, event, subgroup, tau, times)
  arms <- .arm_levels(data[[arm]])
  soc <- "SOC"
  active <- setdiff(arms, soc)
  if (length(active) != 1L) {
    stop("`data` must contain exactly two arms and include `SOC`.", call. = FALSE)
  }

  pieces <- lapply(.group_levels(data, subgroup), function(g) {
    dat <- .subset_group(data, subgroup, g)
    .fit_vs_group(dat, arm, time, event, g, tau, times, soc, active)
  })

  .as_cce_fit(
    mode = "vs",
    label = "Observed comparison",
    curves = do.call(rbind, lapply(pieces, `[[`, "curves")),
    effects = do.call(rbind, lapply(pieces, `[[`, "effects")),
    diagnostics = do.call(rbind, lapply(pieces, `[[`, "diagnostics")),
    meta = list(
      tau = tau,
      times = sort(unique(times)),
      arms = arms
    )
  )
}

#' Project SOC-only survival scenarios
#'
#' Projects comparator curves from the observed `SOC` survival curve using a
#' proportional-hazards assumption.
#'
#' @param data A data frame.
#' @param arm Arm column name.
#' @param time Follow-up time column name.
#' @param event Event indicator column name.
#' @param subgroup Optional subgroup column name.
#' @param soc_level Label used for standard of care.
#' @param tau Horizon for RMST calculation.
#' @param hr_scenarios Numeric vector of hazard-ratio scenarios.
#' @param target_delta_rmst Optional target RMST gain used to back-solve
#'   a required hazard ratio.
#'
#' @return A `cce_fit` object.
#' @export
project_soc_only <- function(data,
                             arm = "arm",
                             time = "time",
                             event = "event",
                             subgroup = NULL,
                             soc_level = "SOC",
                             tau,
                             hr_scenarios = c(0.8, 1.0),
                             target_delta_rmst = NULL) {
  .validate_common_inputs(data, arm, time, event, subgroup, tau, times = tau)
  if (!is.numeric(hr_scenarios) || any(!is.finite(hr_scenarios)) || any(hr_scenarios <= 0)) {
    stop("`hr_scenarios` must be a positive numeric vector.", call. = FALSE)
  }
  if (!soc_level %in% data[[arm]]) {
    stop("`soc_level` was not found in `data[[arm]]`.", call. = FALSE)
  }

  pieces <- lapply(.group_levels(data, subgroup), function(g) {
    dat <- .subset_group(data, subgroup, g)
    .fit_soc_group(dat, arm, time, event, g, tau, hr_scenarios, soc_level, target_delta_rmst)
  })

  .as_cce_fit(
    mode = "soc_only",
    label = "Projection (assumption-based)",
    curves = do.call(rbind, lapply(pieces, `[[`, "curves")),
    effects = do.call(rbind, lapply(pieces, `[[`, "effects")),
    diagnostics = do.call(rbind, lapply(pieces, `[[`, "diagnostics")),
    meta = list(
      tau = tau,
      hr_scenarios = unname(hr_scenarios),
      soc_level = soc_level,
      target_delta_rmst = target_delta_rmst
    )
  )
}

#' Extract the curves table
#'
#' @param x A `cce_fit` object.
#'
#' @return A data frame.
#' @export
as_curves_df <- function(x) {
  .assert_fit(x)
  x$curves
}

#' Extract the effects table
#'
#' @param x A `cce_fit` object.
#'
#' @return A data frame.
#' @export
as_effects_df <- function(x) {
  .assert_fit(x)
  x$effects
}

#' Extract the diagnostics table
#'
#' @param x A `cce_fit` object.
#'
#' @return A data frame.
#' @export
as_diagnostics_df <- function(x) {
  .assert_fit(x)
  x$diagnostics
}

#' Write CCE results to disk
#'
#' @param x A `cce_fit` object.
#' @param dir Output directory.
#'
#' @return Invisibly returns `dir`.
#' @export
write_cce_results <- function(x, dir) {
  .assert_fit(x)
  stopifnot(length(dir) == 1L, is.character(dir), nzchar(dir))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(x$curves, file.path(dir, "curves.csv"), row.names = FALSE)
  utils::write.csv(x$effects, file.path(dir, "effects.csv"), row.names = FALSE)
  utils::write.csv(x$diagnostics, file.path(dir, "diagnostics.csv"), row.names = FALSE)

  payload <- list(
    mode = x$mode,
    label = x$label,
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    meta = x$meta
  )
  jsonlite::write_json(payload, path = file.path(dir, "results.json"), auto_unbox = TRUE, pretty = TRUE)
  invisible(dir)
}

#' @export
print.cce_fit <- function(x, ...) {
  cat("cce_fit\n")
  cat("  mode:", x$mode, "\n")
  cat("  label:", x$label, "\n")
  cat("  curves:", nrow(x$curves), "rows\n")
  cat("  effects:", nrow(x$effects), "rows\n")
  cat("  diagnostics:", nrow(x$diagnostics), "rows\n")
  invisible(x)
}

#' @export
plot.cce_fit <- function(x, subgroup = "All", ...) {
  .assert_fit(x)
  dat <- x$curves[x$curves$subgroup == subgroup, , drop = FALSE]
  if (!nrow(dat)) {
    stop("Requested subgroup was not found.", call. = FALSE)
  }
  arms <- unique(dat$arm)
  cols <- grDevices::hcl.colors(length(arms), palette = "Dark 3")

  plot(
    NA,
    xlim = range(dat$time),
    ylim = c(0, 1),
    xlab = "Time",
    ylab = "Survival",
    main = paste(x$label, "-", subgroup),
    ...
  )

  for (i in seq_along(arms)) {
    arm_i <- arms[[i]]
    d <- dat[dat$arm == arm_i, , drop = FALSE]
    lines(d$time, d$survival, type = "s", lwd = 2, col = cols[[i]])
  }
  legend("topright", legend = arms, col = cols, lwd = 2, bty = "n")
  invisible(x)
}

.validate_common_inputs <- function(data, arm, time, event, subgroup, tau, times) {
  stopifnot(is.data.frame(data))
  cols <- c(arm, time, event)
  if (!is.null(subgroup)) {
    cols <- c(cols, subgroup)
  }
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols)) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (length(tau) != 1L || !is.numeric(tau) || !is.finite(tau) || tau <= 0) {
    stop("`tau` must be one positive number.", call. = FALSE)
  }
  if (!is.numeric(times) || any(!is.finite(times)) || any(times <= 0)) {
    stop("`times` must be a positive numeric vector.", call. = FALSE)
  }
}

.arm_levels <- function(x) {
  levels <- sort(unique(as.character(x)))
  if (length(levels) != 2L) {
    stop("`arm` must contain exactly two unique values.", call. = FALSE)
  }
  levels
}

.group_levels <- function(data, subgroup) {
  if (is.null(subgroup)) {
    return("All")
  }
  c("All", sort(unique(as.character(data[[subgroup]]))))
}

.subset_group <- function(data, subgroup, group_label) {
  if (is.null(subgroup) || identical(group_label, "All")) {
    return(data)
  }
  data[as.character(data[[subgroup]]) == group_label, , drop = FALSE]
}

.fit_vs_group <- function(data, arm, time, event, group_label, tau, times, soc, active) {
  fit <- survival::survfit(stats::as.formula(
    paste0("survival::Surv(", time, ", ", event, ") ~ ", arm)
  ), data = data)

  grid <- .curve_grid(fit, tau, times)
  curves <- .curve_df_from_fit(fit, arm, grid, group_label, mode = "vs", method = "km")

  by_arm <- split(curves, curves$arm)
  surv_soc <- by_arm[[soc]]
  surv_act <- by_arm[[active]]
  rmst_soc <- .rmst_step(surv_soc$time, surv_soc$survival)
  rmst_act <- .rmst_step(surv_act$time, surv_act$survival)

  landmark_rows <- lapply(sort(unique(times)), function(ti) {
    s_soc <- surv_soc$survival[match(ti, surv_soc$time)]
    s_act <- surv_act$survival[match(ti, surv_act$time)]
    data.frame(
      mode = "vs",
      subgroup = group_label,
      contrast = paste(active, "-", soc),
      metric = "delta_survival",
      time = ti,
      estimate = s_act - s_soc,
      stringsAsFactors = FALSE
    )
  })

  effects <- rbind(
    do.call(rbind, landmark_rows),
    data.frame(
      mode = "vs",
      subgroup = group_label,
      contrast = paste(active, "-", soc),
      metric = "delta_rmst",
      time = tau,
      estimate = rmst_act - rmst_soc,
      stringsAsFactors = FALSE
    )
  )

  diagnostics <- .diagnostics_df(data, arm, event, "vs", group_label)

  list(curves = curves, effects = effects, diagnostics = diagnostics)
}

.fit_soc_group <- function(data, arm, time, event, group_label, tau, hr_scenarios, soc_level, target_delta_rmst) {
  soc_data <- data[as.character(data[[arm]]) == soc_level, , drop = FALSE]
  fit <- survival::survfit(stats::as.formula(
    paste0("survival::Surv(", time, ", ", event, ") ~ 1")
  ), data = soc_data)

  grid <- .curve_grid(fit, tau, tau)
  base_curve <- .single_curve_df(fit, grid, group_label, mode = "soc_only", method = "projection_ph", arm_label = soc_level)

  projected_curves <- lapply(hr_scenarios, function(hr) {
    data.frame(
      mode = "soc_only",
      subgroup = group_label,
      method = "projection_ph",
      arm = paste0("HR_", formatC(hr, format = "f", digits = 2)),
      time = base_curve$time,
      survival = base_curve$survival ^ hr,
      n_risk = base_curve$n_risk,
      stringsAsFactors = FALSE
    )
  })
  curves <- rbind(base_curve, do.call(rbind, projected_curves))

  rmst_soc <- .rmst_step(base_curve$time, base_curve$survival)
  effects <- do.call(rbind, lapply(hr_scenarios, function(hr) {
    scenario_label <- paste0("HR_", formatC(hr, format = "f", digits = 2))
    curve_hr <- curves[curves$arm == scenario_label, , drop = FALSE]
    rmst_hr <- .rmst_step(curve_hr$time, curve_hr$survival)
    data.frame(
      mode = "soc_only",
      subgroup = group_label,
      contrast = paste(scenario_label, "-", soc_level),
      metric = "delta_rmst",
      time = tau,
      estimate = rmst_hr - rmst_soc,
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(target_delta_rmst)) {
    required <- .required_hr(base_curve$time, base_curve$survival, target_delta_rmst)
    effects <- rbind(
      effects,
      data.frame(
        mode = "soc_only",
        subgroup = group_label,
        contrast = paste("required_hr", "-", soc_level),
        metric = "required_hr",
        time = tau,
        estimate = required,
        stringsAsFactors = FALSE
      )
    )
  }

  diagnostics <- .diagnostics_df(soc_data, arm, event, "soc_only", group_label)

  list(curves = curves, effects = effects, diagnostics = diagnostics)
}

.curve_grid <- function(fit, tau, times) {
  sort(unique(c(0, fit$time, times, tau)))
}

.curve_df_from_fit <- function(fit, arm_name, grid, subgroup_label, mode, method) {
  s <- summary(fit, times = grid, extend = TRUE)
  strata <- as.character(s$strata)
  arm_value <- sub(paste0("^", arm_name, "="), "", strata)
  data.frame(
    mode = mode,
    subgroup = subgroup_label,
    method = method,
    arm = arm_value,
    time = s$time,
    survival = s$surv,
    n_risk = s$n.risk,
    stringsAsFactors = FALSE
  )
}

.single_curve_df <- function(fit, grid, subgroup_label, mode, method, arm_label) {
  s <- summary(fit, times = grid, extend = TRUE)
  data.frame(
    mode = mode,
    subgroup = subgroup_label,
    method = method,
    arm = arm_label,
    time = s$time,
    survival = s$surv,
    n_risk = s$n.risk,
    stringsAsFactors = FALSE
  )
}

.rmst_step <- function(time, survival) {
  if (length(time) < 2L) {
    return(NA_real_)
  }
  sum(diff(time) * survival[-length(survival)])
}

.required_hr <- function(time, survival, target_delta_rmst) {
  stopifnot(length(target_delta_rmst) == 1L, is.numeric(target_delta_rmst), is.finite(target_delta_rmst))
  rmst_soc <- .rmst_step(time, survival)
  objective <- function(hr) {
    rmst_hr <- .rmst_step(time, survival ^ hr)
    (rmst_hr - rmst_soc) - target_delta_rmst
  }
  lower <- 0.05
  upper <- 2.00
  if (objective(lower) * objective(upper) > 0) {
    return(NA_real_)
  }
  stats::uniroot(objective, interval = c(lower, upper))$root
}

.diagnostics_df <- function(data, arm, event, mode, subgroup_label) {
  split_data <- split(data, as.character(data[[arm]]))
  do.call(rbind, lapply(names(split_data), function(a) {
    dat <- split_data[[a]]
    rbind(
      data.frame(mode = mode, subgroup = subgroup_label, arm = a, metric = "n", value = nrow(dat), stringsAsFactors = FALSE),
      data.frame(mode = mode, subgroup = subgroup_label, arm = a, metric = "events", value = sum(dat[[event]]), stringsAsFactors = FALSE)
    )
  }))
}

.as_cce_fit <- function(mode, label, curves, effects, diagnostics, meta) {
  structure(
    list(
      mode = mode,
      label = label,
      curves = curves,
      effects = effects,
      diagnostics = diagnostics,
      meta = meta
    ),
    class = "cce_fit"
  )
}

.assert_fit <- function(x) {
  if (!inherits(x, "cce_fit")) {
    stop("Expected a `cce_fit` object.", call. = FALSE)
  }
}
