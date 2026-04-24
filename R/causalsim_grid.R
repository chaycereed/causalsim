#' Evaluate an estimator across a parameter grid
#'
#' Runs [causalsim_eval()] over the Cartesian product of one or more DGP
#' parameter values, returning a tidy data frame of performance metrics for
#' every cell. Designed for studying how estimator behavior changes with
#' sample size, confounding strength, effect size, or noise level.
#'
#' @param dgp A `causalsim_dgp` object. Provides fixed parameter values for
#'   all dimensions not listed in `vary`.
#' @param estimator A function as accepted by [causalsim_eval()].
#' @param vary Named list of atomic vectors. Each name must be a valid
#'   [causalsim_dgp()] argument (except `covariates`, which cannot be
#'   varied atomically). Each vector supplies the values to try for that
#'   dimension. The full grid is the Cartesian product of all dimensions.
#' @param reps Positive integer. Replications per grid cell. Default `200L`.
#' @param metrics Character vector. Passed to [causalsim_eval()].
#' @param seed Integer or `NULL`. Passed to [set.seed()] before the first
#'   cell. The same RNG stream continues across cells, so results are
#'   jointly reproducible. Default `NULL`.
#' @param verbose Logical. If `TRUE`, prints a progress message before each
#'   cell. Default `FALSE`.
#'
#' @details
#' ## How it works
#'
#' For each cell in `expand.grid(vary)`, `causalsim_grid()` takes `dgp`'s
#' stored original parameters, overrides the cell's values, reconstructs a
#' new [causalsim_dgp()], runs [causalsim_eval()], and tags the resulting
#' metrics with the cell's parameter values.
#'
#' ## Constraints on `vary`
#'
#' Elements of `vary` must be atomic vectors (character, numeric, integer).
#' Functions cannot be varied via this interface. `covariates` (a list of
#' [covar()] objects) is excluded. For complex covariate variation, construct
#' DGPs manually and use [causalsim_eval()] directly.
#'
#' Note that varying `n_confounders` from 1 to 2 changes auto-generated
#' covariate names from `W` to `W1, W2`. If the base DGP's `effect` function
#' references `W`, the reconstructed DGP will error at construction time.
#' Design the base DGP accordingly.
#'
#' @return An S3 object of class `causalsim_grid` with components:
#' \describe{
#'   \item{`results`}{Tidy data frame: one row per (cell, metric).
#'     Columns: grid parameter names, `metric`, `value`, `se`.}
#'   \item{`grid`}{Data frame of grid cells, one row per cell.}
#'   \item{`vary`}{Character vector of varied parameter names.}
#'   \item{`reps`}{Replications per cell.}
#'   \item{`metrics`}{Metrics evaluated.}
#' }
#'
#' @examples
#' dgp <- causalsim_dgp(n = 500, n_confounders = 1, effect = 2,
#'                      propensity = 0.5)
#'
#' ols_estimator <- function(data) {
#'   fit <- lm(Y ~ A + W, data = data)
#'   est <- coef(fit)[["A"]]
#'   se  <- sqrt(vcov(fit)["A", "A"])
#'   c(estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
#' }
#'
#' grid_result <- causalsim_grid(
#'   dgp       = dgp,
#'   estimator = ols_estimator,
#'   vary      = list(n = c(250L, 500L, 1000L)),
#'   reps      = 50L,
#'   seed      = 1L
#' )
#' grid_result
#'
#' @export
causalsim_grid <- function(
  dgp,
  estimator,
  vary,
  reps    = 200L,
  metrics = c("bias", "rmse", "coverage", "power"),
  seed    = NULL,
  verbose = FALSE
) {
  if (!inherits(dgp, "causalsim_dgp")) {
    stop("`dgp` must be a causalsim_dgp object.", call. = FALSE)
  }
  if (!is.function(estimator)) {
    stop("`estimator` must be a function.", call. = FALSE)
  }
  reps <- .validate_positive(reps, "reps", as_int = TRUE)

  .validate_vary(vary, dgp)

  grid_df  <- expand.grid(vary, stringsAsFactors = FALSE)
  n_cells  <- nrow(grid_df)
  vary_nms <- names(vary)

  metrics <- match.arg(
    metrics,
    choices    = c("bias", "rmse", "coverage", "power"),
    several.ok = TRUE
  )

  if (!is.null(seed)) set.seed(seed)

  all_rows <- vector("list", n_cells)
  for (i in seq_len(n_cells)) {
    if (verbose) {
      message(sprintf("Grid cell %d / %d", i, n_cells))
    }

    cell_params <- dgp$params
    for (nm in vary_nms) {
      cell_params[[nm]] <- grid_df[[nm]][[i]]
    }

    cell_dgp    <- do.call(causalsim_dgp, cell_params)
    cell_eval   <- causalsim_eval(cell_dgp, estimator,
                                  reps = reps, metrics = metrics)
    cell_result <- cell_eval$metrics

    for (nm in vary_nms) {
      cell_result[[nm]] <- grid_df[[nm]][[i]]
    }

    all_rows[[i]] <- cell_result
  }

  results          <- do.call(rbind, all_rows)
  results          <- results[, c(vary_nms, "metric", "value", "se")]
  rownames(results) <- NULL

  structure(
    list(
      results = results,
      grid    = grid_df,
      vary    = vary_nms,
      reps    = reps,
      metrics = metrics
    ),
    class = "causalsim_grid"
  )
}

#' @export
print.causalsim_grid <- function(x, n = 10L, ...) {
  n_cells <- nrow(x$grid)
  cat(sprintf(
    "<causalsim_grid>  %d cell%s  vary: %s  reps/cell: %d\n",
    n_cells,
    if (n_cells == 1L) "" else "s",
    paste(x$vary, collapse = ", "),
    x$reps
  ))
  cat(sprintf("  metrics: %s\n\n", paste(x$metrics, collapse = ", ")))
  nr <- nrow(x$results)
  print(x$results[seq_len(min(n, nr)), ], row.names = FALSE)
  if (nr > n) {
    cat(sprintf(
      "\n  ... %d more row%s. Access full results via $results.\n",
      nr - n, if (nr - n == 1L) "" else "s"
    ))
  }
  invisible(x)
}

# ── Internal ──────────────────────────────────────────────────────────────────

.validate_vary <- function(vary, dgp) {
  if (!is.list(vary) || length(vary) == 0L) {
    stop(
      "`vary` must be a non-empty named list of atomic vectors.",
      call. = FALSE
    )
  }
  if (is.null(names(vary)) || any(names(vary) == "")) {
    stop("Every element of `vary` must be named.", call. = FALSE)
  }

  valid_nms <- setdiff(names(formals(causalsim_dgp)), "covariates")
  invalid   <- setdiff(names(vary), valid_nms)
  if (length(invalid) > 0L) {
    stop(
      paste0(
        "`vary` contains name(s) that are not a valid causalsim_dgp() argument: ",
        paste(invalid, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  not_atomic <- !vapply(vary, is.atomic, logical(1L))
  if (any(not_atomic)) {
    stop(
      paste0(
        "All elements of `vary` must be atomic vectors. Non-atomic: ",
        paste(names(vary)[not_atomic], collapse = ", ")
      ),
      call. = FALSE
    )
  }
}
