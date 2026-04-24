#' Evaluate a causal estimator against a known DGP
#'
#' Repeatedly draws datasets from a [causalsim_dgp()] object, applies a
#' user-supplied estimator to each draw, and returns a tidy summary of
#' estimator performance against the known ground truth.
#'
#' @param dgp A `causalsim_dgp` object.
#' @param estimator A function that accepts a data frame produced by
#'   [causalsim_draw()] and returns a named numeric vector or single-row
#'   data frame. Must include an `estimate` field. For `"coverage"` and
#'   `"power"` metrics, must also include `ci_lower` and `ci_upper`. See
#'   Details.
#' @param reps Positive integer. Number of simulation replications.
#'   Default `200L`.
#' @param metrics Character vector. Subset of `"bias"`, `"rmse"`,
#'   `"coverage"`, `"power"`. Default: all four.
#' @param seed Integer or `NULL`. Passed to [set.seed()] before the first
#'   replication. Default `NULL`.
#'
#' @details
#' ## Estimator convention
#'
#' The estimator receives the full data frame returned by [causalsim_draw()],
#' including ground-truth columns `.tau` and `.p`. It should return a named
#' numeric vector or one-row data frame with at minimum:
#'
#' ```r
#' my_estimator <- function(data) {
#'   fit <- lm(Y ~ A + W, data = data)
#'   est <- coef(fit)["A"]
#'   se  <- sqrt(vcov(fit)["A", "A"])
#'   c(estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
#' }
#' ```
#'
#' ## Metrics
#'
#' \describe{
#'   \item{`"bias"`}{Monte Carlo estimate of bias: `mean(estimate) - true_ate`.}
#'   \item{`"rmse"`}{Root mean squared error: `sqrt(mean((estimate - true_ate)^2))`.}
#'   \item{`"coverage"`}{Proportion of replications where the CI brackets
#'     `true_ate`. Requires `ci_lower` and `ci_upper`.}
#'   \item{`"power"`}{Proportion of replications where the CI excludes zero
#'     — i.e., the rejection rate of H0: ATE = 0. When `true_ate != 0`,
#'     this is power; when `true_ate == 0`, it is the Type I error rate.
#'     Requires `ci_lower` and `ci_upper`.}
#' }
#'
#' Monte Carlo standard errors (MCSE) are reported alongside each metric.
#' For bias: `sd(estimates) / sqrt(reps)`. For RMSE: delta-method
#' approximation. For coverage and power: Bernoulli SE.
#'
#' ## Raw draws
#'
#' Per-replication estimates are stored in `result$draws` and can be used
#' for custom analysis or plotting.
#'
#' @return An S3 object of class `causalsim_eval` with components:
#' \describe{
#'   \item{`metrics`}{Tidy data frame: `metric`, `value`, `se`.}
#'   \item{`draws`}{Data frame of per-replication estimates, one row per rep.}
#'   \item{`true_ate`}{The DGP's true ATE.}
#'   \item{`reps`}{Number of replications run.}
#' }
#'
#' @examples
#' dgp <- causalsim_dgp(n = 300, n_confounders = 1, effect = 2,
#'                      propensity = 0.5)
#'
#' # OLS estimator (unbiased under RCT)
#' ols_estimator <- function(data) {
#'   fit <- lm(Y ~ A + W, data = data)
#'   est <- coef(fit)[["A"]]
#'   se  <- sqrt(vcov(fit)["A", "A"])
#'   c(estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
#' }
#'
#' result <- causalsim_eval(dgp, ols_estimator, reps = 100L, seed = 1L)
#' result
#'
#' @export
causalsim_eval <- function(
  dgp,
  estimator,
  reps    = 200L,
  metrics = c("bias", "rmse", "coverage", "power"),
  seed    = NULL
) {
  if (!inherits(dgp, "causalsim_dgp")) {
    stop("`dgp` must be a causalsim_dgp object.", call. = FALSE)
  }
  if (!is.function(estimator)) {
    stop("`estimator` must be a function.", call. = FALSE)
  }
  reps    <- .validate_positive(reps, "reps", as_int = TRUE)
  metrics <- match.arg(
    metrics,
    choices    = c("bias", "rmse", "coverage", "power"),
    several.ok = TRUE
  )

  if (!is.null(seed)) set.seed(seed)

  draws <- .collect_draws(dgp, estimator, reps)

  needs_ci <- any(metrics %in% c("coverage", "power"))
  if (needs_ci && !all(c("ci_lower", "ci_upper") %in% names(draws))) {
    stop(
      paste0(
        "Metrics 'coverage' and 'power' require the estimator to return ",
        "'ci_lower' and 'ci_upper'."
      ),
      call. = FALSE
    )
  }

  metrics_df <- .compute_metrics(draws, dgp$true_ate, metrics)

  structure(
    list(
      metrics  = metrics_df,
      draws    = draws,
      true_ate = dgp$true_ate,
      reps     = reps
    ),
    class = "causalsim_eval"
  )
}

#' @export
print.causalsim_eval <- function(x, digits = 4L, ...) {
  cat(sprintf(
    "<causalsim_eval>  reps: %d  true ATE: %s\n\n",
    x$reps, format(round(x$true_ate, digits))
  ))
  m      <- x$metrics
  m$value <- round(m$value, digits)
  m$se    <- round(m$se,    digits + 1L)
  print(m, row.names = FALSE)
  invisible(x)
}

#' Summarise a causalsim_eval result
#'
#' Returns a structured summary of estimator performance including the metrics
#' table and the distribution of per-replication estimates (mean, SD, median,
#' and 10th / 90th percentiles).
#'
#' @param object A `causalsim_eval` object.
#' @param ... Ignored.
#'
#' @return An S3 object of class `causalsim_eval_summary` with components:
#' \describe{
#'   \item{`metrics`}{The metrics data frame from the original result.}
#'   \item{`mean_estimate`}{Mean of per-replication estimates.}
#'   \item{`sd_estimate`}{Standard deviation of per-replication estimates.}
#'   \item{`median_estimate`}{Median of per-replication estimates.}
#'   \item{`p10_estimate`}{10th percentile of per-replication estimates.}
#'   \item{`p90_estimate`}{90th percentile of per-replication estimates.}
#'   \item{`true_ate`}{True ATE from the DGP.}
#'   \item{`reps`}{Number of replications.}
#' }
#'
#' @export
summary.causalsim_eval <- function(object, ...) {
  est <- object$draws$estimate
  structure(
    list(
      metrics         = object$metrics,
      mean_estimate   = mean(est),
      sd_estimate     = stats::sd(est),
      median_estimate = stats::median(est),
      p10_estimate    = unname(stats::quantile(est, 0.10)),
      p90_estimate    = unname(stats::quantile(est, 0.90)),
      true_ate        = object$true_ate,
      reps            = object$reps
    ),
    class = "causalsim_eval_summary"
  )
}

#' @export
print.causalsim_eval_summary <- function(x, digits = 4L, ...) {
  cat(sprintf(
    "<causalsim_eval_summary>  reps: %d  true ATE: %s\n\n",
    x$reps, format(round(x$true_ate, digits))
  ))
  cat("Estimate distribution:\n")
  cat(sprintf("  mean  : %s\n",         format(round(x$mean_estimate,   digits))))
  cat(sprintf("  sd    : %s\n",         format(round(x$sd_estimate,     digits))))
  cat(sprintf("  median: %s\n",         format(round(x$median_estimate, digits))))
  cat(sprintf("  [p10, p90]: [%s, %s]\n\n",
              format(round(x$p10_estimate, digits)),
              format(round(x$p90_estimate, digits))))
  cat("Metrics:\n")
  m        <- x$metrics
  m$value  <- round(m$value, digits)
  m$se     <- round(m$se,    digits + 1L)
  print(m, row.names = FALSE)
  invisible(x)
}

#' Plot the distribution of estimates from a causalsim_eval result
#'
#' Draws a histogram of per-replication estimates with vertical lines marking
#' the true ATE (solid) and the mean estimate (dashed).
#'
#' @param x A `causalsim_eval` object.
#' @param ... Additional arguments passed to [graphics::hist()].
#'
#' @return `x`, invisibly.
#'
#' @export
plot.causalsim_eval <- function(x, ...) {
  est      <- x$draws$estimate
  true_ate <- x$true_ate
  mean_est <- mean(est)

  graphics::hist(
    est,
    main   = "Distribution of Estimates",
    xlab   = "Estimate",
    col    = "grey85",
    border = "white",
    ...
  )
  graphics::abline(v = true_ate, lwd = 2)
  graphics::abline(v = mean_est, lwd = 1, lty = 2)
  graphics::legend(
    "topright",
    legend = c(
      sprintf("true ATE = %.3f", true_ate),
      sprintf("mean est = %.3f", mean_est)
    ),
    lty = c(1L, 2L),
    lwd = c(2L, 1L),
    bty = "n"
  )
  invisible(x)
}
