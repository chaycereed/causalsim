# Null-coalescing: return x unless it is NULL, then return y
`%||%` <- function(x, y) if (is.null(x)) y else x

# ── Input validation ──────────────────────────────────────────────────────────

.validate_positive <- function(x, name, as_int = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a finite positive number.", name), call. = FALSE)
  }
  if (as_int) as.integer(x) else x
}

.validate_nonneg <- function(x, name, as_int = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
    stop(
      sprintf("`%s` must be a finite non-negative number.", name),
      call. = FALSE
    )
  }
  if (as_int) as.integer(x) else x
}

# Keep old names as aliases so existing callers continue to work
.validate_positive_int <- function(x, name) {
  .validate_positive(x, name, as_int = TRUE)
}
.validate_nonneg_int <- function(x, name) {
  .validate_nonneg(x, name, as_int = TRUE)
}

# ── Covariate spec construction ───────────────────────────────────────────────

# Merge shorthand count arguments with an explicit named covariate list.
.build_covariate_spec <- function(n_confounders, n_effect_modifiers,
                                   n_instruments, n_noise, covariates) {
  auto <- c(
    .auto_covariates(n_confounders,      "confounder",      "W"),
    .auto_covariates(n_effect_modifiers, "effect_modifier", "V"),
    .auto_covariates(n_instruments,      "instrument",      "Z"),
    .auto_covariates(n_noise,            "noise",           "X")
  )

  if (!is.list(covariates)) {
    stop("`covariates` must be a named list of covar() objects.", call. = FALSE)
  }
  if (length(covariates) > 0L) {
    if (is.null(names(covariates)) || any(names(covariates) == "")) {
      stop("Every element of `covariates` must be named.", call. = FALSE)
    }
    not_covar <- !vapply(
      covariates, inherits, logical(1L), what = "causalsim_covar"
    )
    if (any(not_covar)) {
      stop(
        "All elements of `covariates` must be covar() objects.",
        call. = FALSE
      )
    }
    overlap <- intersect(names(auto), names(covariates))
    if (length(overlap) > 0L) {
      stop(
        sprintf(
          paste0(
            "Covariate name(s) appear in both auto-generated",
            " and explicit spec: %s"
          ),
          paste(overlap, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  c(auto, covariates)
}

# Naming convention: single item → bare prefix (W, Z, V, X);
# multiple items → numbered suffix (W1, W2, ...).
.auto_covariates <- function(n, role, prefix) {
  if (n == 0L) return(list())
  nms  <- if (n == 1L) prefix else paste0(prefix, seq_len(n))
  covs <- replicate(n, covar(dist = "normal", role = role), simplify = FALSE)
  stats::setNames(covs, nms)
}

# ── Covariate data generation ─────────────────────────────────────────────────

# Generate an n-row data frame from a named covariate spec list.
.generate_covariate_data <- function(covar_spec, n) {
  if (length(covar_spec) == 0L) return(data.frame())
  cols <- lapply(covar_spec, .draw_one_covariate, n = n)
  as.data.frame(cols, check.names = FALSE)
}

# Draw n observations for a single causalsim_covar object.
.draw_one_covariate <- function(cv, n) {
  p <- cv$params
  switch(cv$dist,
    normal  = stats::rnorm(n,  mean = p$mean %||% 0,   sd   = p$sd   %||% 1),
    binary  = stats::rbinom(n, size = 1L, prob = p$prob %||% 0.5),
    uniform = stats::runif(n,  min  = p$min  %||% 0,   max  = p$max  %||% 1),
    stop(sprintf("Unknown distribution: '%s'", cv$dist), call. = FALSE)
  )
}

# ── Function normalization ────────────────────────────────────────────────────

# Preset table: maps (type, level) → linear coefficient on confounder columns.
# Propensity uses plogis(); baseline uses a linear sum.
.fn_presets <- list(
  propensity = c(low = 0.25, moderate = 0.5, high = 1.0),
  baseline   = c(low = 0.25, moderate = 0.5, high = 1.0)
)

# Resolve a preset string to a closure over the covariate spec.
# The returned function uses `...` and expects to be called via .apply_effect.
.resolve_preset <- function(x, type, covar_spec) {
  table <- .fn_presets[[type]]
  if (is.null(table)) {
    stop(
      sprintf(
        "No presets defined for `%s`. Use a scalar or function.", type
      ),
      call. = FALSE
    )
  }
  if (!x %in% names(table)) {
    stop(
      sprintf(
        "`%s` preset '%s' not recognized. Choose from: %s.",
        type, x, paste(names(table), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  coeff         <- unname(table[[x]])
  confounder_nms <- names(
    Filter(function(cv) "confounder" %in% cv$role, covar_spec)
  )

  # No confounders: preset has no structural effect
  if (length(confounder_nms) == 0L) {
    default_val <- if (type == "propensity") 0.5 else 0
    return(function(...) default_val)
  }

  # Return closure capturing coeff and confounder names
  coeff_         <- coeff
  confounder_nms_ <- confounder_nms
  if (type == "propensity") {
    function(...) {
      args <- list(...)
      terms <- lapply(confounder_nms_, function(nm) coeff_ * args[[nm]])
      lin   <- Reduce(`+`, terms)
      stats::plogis(lin)
    }
  } else {
    function(...) {
      args <- list(...)
      Reduce(`+`, lapply(confounder_nms_, function(nm) coeff_ * args[[nm]]))
    }
  }
}

# Normalize effect/propensity/baseline to a callable function.
# Accepts a scalar (constant), a preset string, or a user-supplied function.
.normalize_fn <- function(x, type, covar_spec) {
  if (is.numeric(x) && length(x) == 1L && is.finite(x)) {
    val <- x
    return(function(...) val)
  }
  if (is.character(x) && length(x) == 1L) {
    return(.resolve_preset(x, type, covar_spec))
  }
  if (is.function(x)) {
    return(x)
  }
  stop(
    sprintf("`%s` must be a scalar, preset string, or function.", type),
    call. = FALSE
  )
}

# ── Function application ──────────────────────────────────────────────────────

# Apply a DGP function to a covariate data frame, always returning a
# length-n vector.
#
# Dispatch rules:
#   no-arg function  (formals = "")       -> scalar, repeat n times
#   ...-only function (formals = "...")   -> do.call with full covariate df
#   named-arg function                    -> do.call with matching columns only
#
# The `n` argument is used when covariate_df has 0 rows (no-covariate DGPs).
.apply_effect <- function(effect_fn, covariate_df, n = nrow(covariate_df)) {
  fn_formals  <- names(formals(effect_fn))
  has_dots    <- "..." %in% fn_formals
  named_args  <- setdiff(fn_formals, "...")

  result <- if (!has_dots && length(named_args) == 0L) {
    effect_fn()
  } else if (has_dots && length(named_args) == 0L) {
    do.call(effect_fn, as.list(covariate_df))
  } else {
    do.call(effect_fn, as.list(covariate_df[named_args]))
  }

  # Recycle scalars to length n so callers always get a vector
  if (length(result) == 1L && n > 1L) rep(result, n) else result
}

# ── Argument validation ───────────────────────────────────────────────────────

# Validate that every named argument of fn exists in covariate_names.
# Skips ...-only and no-arg functions (no names to validate).
.validate_fn_args <- function(fn, covariate_names, arg_label) {
  named_args <- setdiff(names(formals(fn)), "...")
  if (length(named_args) == 0L) return(invisible(NULL))
  missing_vars <- setdiff(named_args, covariate_names)
  if (length(missing_vars) > 0L) {
    stop(
      sprintf(
        "`%s` function references covariate(s) not defined in the DGP: %s",
        arg_label, paste(missing_vars, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

# Validate that the propensity function returns values strictly in [0, 1]
# by evaluating it on a small test draw at construction time.
.validate_propensity_fn <- function(propensity_fn, covar_spec,
                                    n_check = 50L) {
  test_df <- if (length(covar_spec) > 0L) {
    .generate_covariate_data(covar_spec, n = n_check)
  } else {
    data.frame()
  }

  probs <- tryCatch(
    .apply_effect(propensity_fn, test_df, n = n_check),
    error = function(e) {
      stop(
        "Failed to evaluate `propensity` function: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (length(probs) > 0L &&
        (any(!is.finite(probs)) || any(probs < 0) || any(probs > 1))) {
    stop(
      sprintf(
        paste0("`propensity` function must return values in [0, 1]. ",
               "Got range [%.3g, %.3g]."),
        min(probs, na.rm = TRUE), max(probs, na.rm = TRUE)
      ),
      call. = FALSE
    )
  }
}

# ── Monte Carlo ATE ───────────────────────────────────────────────────────────

# Approximate the true ATE via Monte Carlo over the covariate distribution.
# For non-heterogeneous (scalar) effects, returns the exact value immediately.
.mc_ate <- function(effect_fn, covar_spec, mc_draws, heterogeneous) {
  if (!heterogeneous) return(effect_fn())
  covariate_df <- .generate_covariate_data(covar_spec, n = mc_draws)
  mean(.apply_effect(effect_fn, covariate_df, n = mc_draws))
}

# ── Estimator evaluation ──────────────────────────────────────────────────────

# Call a user estimator on one dataset and return a named numeric vector.
# Accepts a named numeric vector, a named all-numeric list, or a one-row data
# frame. Errors with a clear message if 'estimate' is not present.
.call_estimator <- function(estimator, data) {
  result <- estimator(data)
  if (is.data.frame(result)) {
    if (nrow(result) != 1L) {
      stop(
        "estimator must return a named numeric vector, named list, or ",
        "one-row data frame.", call. = FALSE
      )
    }
    result <- unlist(result[1L, ], use.names = TRUE)
  }
  if (is.list(result)) {
    all_scalar_numeric <- !is.null(names(result)) &&
      all(vapply(result, function(x) is.numeric(x) && length(x) == 1L,
                 logical(1L)))
    if (!all_scalar_numeric) {
      stop(
        "estimator returned a list, but not all elements are named numeric ",
        "scalars. Use c() or return a one-row data frame instead.",
        call. = FALSE
      )
    }
    result <- stats::setNames(as.numeric(result), names(result))
  }
  if (!is.numeric(result) || is.null(names(result))) {
    stop(
      sprintf(
        "estimator must return a named numeric vector, named list, or ",
        "one-row data frame; got %s.", class(result)[1L]
      ),
      call. = FALSE
    )
  }
  if (!"estimate" %in% names(result)) {
    stop(
      "estimator output must include a field named 'estimate'.",
      call. = FALSE
    )
  }
  result
}

# Run `reps` draws from `dgp`, apply `estimator` to each, and return a
# data frame with one row per replication.
.collect_draws <- function(dgp, estimator, reps) {
  rows <- vector("list", reps)
  for (i in seq_len(reps)) {
    data    <- causalsim_draw(dgp)
    est     <- .call_estimator(estimator, data)
    rows[[i]] <- as.list(est)
  }
  # Bind rows; missing fields (e.g. ci_lower absent in some reps) become NA
  all_nms <- unique(unlist(lapply(rows, names)))
  out <- as.data.frame(
    do.call(rbind, lapply(rows, function(r) {
      r[setdiff(all_nms, names(r))] <- NA_real_
      as.data.frame(r[all_nms])
    })),
    stringsAsFactors = FALSE
  )
  out
}

# Compute requested metrics from a draws data frame and the true ATE.
# Returns a tidy data frame: metric | value | se
.compute_metrics <- function(draws, true_ate, metrics) {
  est  <- draws[["estimate"]]
  errs <- est - true_ate

  rows <- lapply(metrics, function(m) {
    switch(m,
      bias = {
        val <- mean(errs)
        se  <- stats::sd(est) / sqrt(length(est))
        data.frame(metric = "bias", value = val, se = se,
                   stringsAsFactors = FALSE)
      },
      rmse = {
        sq  <- errs^2
        val <- sqrt(mean(sq))
        # Delta method: SE(sqrt(mean(X))) ≈ SD(X) / (2*sqrt(mean(X))*sqrt(n))
        se  <- if (val > 0) {
          stats::sd(sq) / (2 * val * sqrt(length(sq)))
        } else {
          0
        }
        data.frame(metric = "rmse", value = val, se = se,
                   stringsAsFactors = FALSE)
      },
      coverage = {
        covered <- draws[["ci_lower"]] <= true_ate &
                   true_ate            <= draws[["ci_upper"]]
        p   <- mean(covered)
        se  <- sqrt(p * (1 - p) / length(covered))
        data.frame(metric = "coverage", value = p, se = se,
                   stringsAsFactors = FALSE)
      },
      power = {
        rejected <- draws[["ci_lower"]] > 0 | draws[["ci_upper"]] < 0
        p  <- mean(rejected)
        se <- sqrt(p * (1 - p) / length(rejected))
        data.frame(metric = "power", value = p, se = se,
                   stringsAsFactors = FALSE)
      }
    )
  })

  do.call(rbind, rows)
}
