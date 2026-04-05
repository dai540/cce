assert_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
}

assert_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("`%s` must be a single non-empty character string.", name),
      call. = FALSE
    )
  }
}

assert_character_vector <- function(x, name, min_len = 1L) {
  if (!is.character(x) || length(x) < min_len || anyNA(x) || any(!nzchar(x))) {
    stop(sprintf("`%s` must be a character vector of length >= %d.", name, min_len),
      call. = FALSE
    )
  }
}

assert_data_frame <- function(x, name) {
  if (!is.data.frame(x)) {
    stop(sprintf("`%s` must be a data.frame.", name), call. = FALSE)
  }
}

assert_named_columns <- function(data, cols, name) {
  missing_cols <- setdiff(cols, colnames(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "`%s` is missing required columns: %s.",
      name,
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }
}

assert_scalar_numeric <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower || x > upper) {
    stop(sprintf("`%s` must be a single numeric value in [%s, %s].",
      name, lower, upper
    ), call. = FALSE)
  }
}

normalize_date <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }
  as.Date(x)
}

safe_factor <- function(x, levels = NULL) {
  if (is.factor(x)) {
    if (!is.null(levels)) {
      return(factor(as.character(x), levels = levels))
    }
    return(x)
  }
  factor(x, levels = levels)
}

evaluate_step <- function(step_times, step_values, at, initial = 0) {
  if (length(step_times) == 0L || length(step_values) == 0L) {
    return(rep(initial, length(at)))
  }
  idx <- findInterval(at, step_times)
  out <- rep(initial, length(at))
  non_zero <- idx > 0L
  out[non_zero] <- step_values[idx[non_zero]]
  out
}

rmst_from_curve <- function(times, surv, tau) {
  grid <- sort(unique(c(0, times[times <= tau], tau)))
  if (length(grid) < 2L) {
    return(0)
  }
  step_vals <- evaluate_step(times, surv, grid[-length(grid)], initial = 1)
  sum(diff(grid) * step_vals)
}

weighted_mean_var <- function(x, w) {
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    return(list(mean = NA_real_, var = NA_real_))
  }
  mu <- sum(w * x) / sw
  vv <- sum(w * (x - mu)^2) / sw
  list(mean = mu, var = vv)
}

covariate_matrix <- function(data, covariates) {
  if (length(covariates) == 0L) {
    return(matrix(numeric(0), nrow = nrow(data), ncol = 0L))
  }
  mm <- stats::model.matrix(stats::reformulate(covariates), data = data)
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}

compute_smd <- function(data, arm, covariates, weights = NULL) {
  x <- covariate_matrix(data, covariates)
  if (ncol(x) == 0L) {
    return(stats::setNames(numeric(0), character(0)))
  }
  g <- as.integer(data[[arm]] == levels(data[[arm]])[2L])
  if (is.null(weights)) {
    weights <- rep(1, nrow(data))
  }
  out <- numeric(ncol(x))
  names(out) <- colnames(x)
  for (j in seq_len(ncol(x))) {
    stats0 <- weighted_mean_var(x[g == 0L, j], weights[g == 0L])
    stats1 <- weighted_mean_var(x[g == 1L, j], weights[g == 1L])
    pooled <- sqrt((stats0$var + stats1$var) / 2)
    out[j] <- if (is.na(pooled) || pooled == 0) 0 else (stats1$mean - stats0$mean) / pooled
  }
  out
}

ess <- function(weights) {
  sw <- sum(weights)
  sw2 <- sum(weights^2)
  if (!is.finite(sw) || !is.finite(sw2) || sw2 <= 0) {
    return(NA_real_)
  }
  (sw^2) / sw2
}

default_time_grid <- function(time, tau = NULL, n_grid = 100L) {
  max_time <- max(time, na.rm = TRUE)
  if (!is.null(tau)) {
    max_time <- min(max_time, tau)
  }
  unique(sort(c(0, seq(0, max_time, length.out = n_grid), max_time)))
}

standardized_survival_from_cox <- function(model, template_data, arm, arm_level, times) {
  bh <- survival::basehaz(model, centered = FALSE)
  hazard <- evaluate_step(bh$time, bh$hazard, times, initial = 0)
  newdata <- template_data
  newdata[[arm]] <- factor(rep(arm_level, nrow(template_data)), levels = levels(template_data[[arm]]))
  lp <- stats::predict(model, newdata = newdata, type = "lp")
  surv_mat <- exp(-outer(exp(lp), hazard))
  colMeans(surv_mat)
}

effect_rows_from_curves <- function(curves, subgroup, tau, landmark_times) {
  methods <- unique(curves$method)
  out <- vector("list", length(methods))
  arm_levels <- unique(curves$arm)
  if (all(c("SOC", "A") %in% arm_levels)) {
    arm_levels <- c("SOC", "A")
  }
  for (i in seq_along(methods)) {
    method_df <- curves[curves$method == methods[i] & curves$subgroup == subgroup, , drop = FALSE]
    soc_df <- method_df[method_df$arm == arm_levels[1L], , drop = FALSE]
    trt_df <- method_df[method_df$arm == arm_levels[2L], , drop = FALSE]
    rmst_soc <- rmst_from_curve(soc_df$time, soc_df$survival, tau)
    rmst_trt <- rmst_from_curve(trt_df$time, trt_df$survival, tau)
    rows <- vector("list", length(landmark_times))
    for (j in seq_along(landmark_times)) {
      tt <- landmark_times[j]
      s_soc <- evaluate_step(soc_df$time, soc_df$survival, tt, initial = 1)
      s_trt <- evaluate_step(trt_df$time, trt_df$survival, tt, initial = 1)
      rows[[j]] <- data.frame(
        mode = unique(method_df$mode),
        method = methods[i],
        subgroup = subgroup,
        tau = tau,
        rmst_arm0 = rmst_soc,
        rmst_arm1 = rmst_trt,
        delta_rmst = rmst_trt - rmst_soc,
        landmark_time = tt,
        survival_arm0 = s_soc,
        survival_arm1 = s_trt,
        delta_survival = s_trt - s_soc,
        stringsAsFactors = FALSE
      )
    }
    out[[i]] <- do.call(rbind, rows)
  }
  do.call(rbind, out)
}

bootstrap_quantiles <- function(df, key_cols, value_cols) {
  split_key <- do.call(paste, c(df[key_cols], sep = "\r"))
  parts <- split(df, split_key)
  out <- vector("list", length(parts))
  nms <- names(parts)
  for (i in seq_along(parts)) {
    part <- parts[[i]]
    base <- part[1L, key_cols, drop = FALSE]
    for (col in value_cols) {
      qs <- stats::quantile(part[[col]], probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 8)
      base[[paste0(col, "_lower_ci")]] <- qs[1L]
      base[[paste0(col, "_upper_ci")]] <- qs[2L]
    }
    out[[i]] <- base
  }
  do.call(rbind, out)
}

merge_interval_columns <- function(target, source, key_cols, value_col) {
  key_target <- do.call(paste, c(target[key_cols], sep = "\r"))
  key_source <- do.call(paste, c(source[key_cols], sep = "\r"))
  idx <- match(key_target, key_source)
  target[[paste0(value_col, "_lower_ci")]] <- source[[paste0(value_col, "_lower_ci")]][idx]
  target[[paste0(value_col, "_upper_ci")]] <- source[[paste0(value_col, "_upper_ci")]][idx]
  target
}

merge_curve_intervals <- function(curves, intervals) {
  curves <- merge_interval_columns(curves, intervals, c("mode", "method", "subgroup", "arm", "time"), "survival")
  names(curves)[names(curves) == "survival_lower_ci"] <- "lower_ci"
  names(curves)[names(curves) == "survival_upper_ci"] <- "upper_ci"
  curves
}

merge_effect_intervals <- function(effects, intervals) {
  effects <- merge_interval_columns(effects, intervals, c("mode", "method", "subgroup", "tau", "landmark_time"), "delta_rmst")
  effects <- merge_interval_columns(effects, intervals, c("mode", "method", "subgroup", "tau", "landmark_time"), "delta_survival")
  effects
}

label_from_fail_flags <- function(fail_flags) {
  if (any(isTRUE(unlist(fail_flags)))) {
    return("reference-only")
  }
  "ok"
}

run_stamp <- function() {
  format(Sys.time(), "%Y%m%d-%H%M%S")
}

weighted_km_survival <- function(data, time, event, weights, times) {
  data$km_weight <- weights
  fit <- survival::survfit(
    make_surv_formula(time, event, "1"),
    data = data,
    weights = km_weight
  )
  evaluate_step(fit$time, fit$surv, times, initial = 1)
}

named_count_table <- function(x) {
  tab <- table(x, useNA = "ifany")
  out <- as.list(as.integer(tab))
  names(out) <- names(tab)
  out
}

capture_dataset_attributes <- function(data) {
  list(
    spec = attr(data, "spec"),
    exclusions = attr(data, "exclusions"),
    validation_report = attr(data, "validation_report"),
    profile = attr(data, "profile")
  )
}

profile_list_for_meta <- function(profile) {
  if (inherits(profile, "cce_profile")) {
    return(unclass(profile))
  }
  profile
}
