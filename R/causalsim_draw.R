#' Draw a dataset from a causal DGP
#'
#' Simulates one dataset from a [causalsim_dgp()] object using the structural
#' model:
#'
#' ```
#' W ~ covariate_spec
#' A ~ Bernoulli(propensity(W))
#' Y  = baseline(W) + effect(W) * A + N(0, sigma^2)
#' ```
#'
#' @param dgp A `causalsim_dgp` object created by [causalsim_dgp()].
#' @param seed Integer or `NULL`. If non-null, passed to [set.seed()] before
#'   generating any random values, making the draw reproducible. Default `NULL`.
#'
#' @return A data frame with `dgp$n` rows and the following columns:
#' \describe{
#'   \item{Covariate columns}{One column per covariate in `dgp$covar_spec`,
#'     named to match the spec (e.g. `W`, `W1`, `Z`).}
#'   \item{`A`}{Binary treatment indicator (0/1).}
#'   \item{`Y`}{Observed outcome.}
#'   \item{`.tau`}{Individual treatment effect (CATE) — ground truth.}
#'   \item{`.p`}{Individual propensity score — ground truth.}
#' }
#'
#' The `.tau` and `.p` columns carry individual-level ground truth and are
#' prefixed with `.` to distinguish them from observed variables.
#' `causalsim_eval()` uses them directly without re-estimation.
#'
#' @examples
#' dgp <- causalsim_dgp(n = 500, n_confounders = 1, effect = 2)
#' d   <- causalsim_draw(dgp, seed = 1)
#' head(d)
#'
#' # Reproducible draws
#' d1 <- causalsim_draw(dgp, seed = 42)
#' d2 <- causalsim_draw(dgp, seed = 42)
#' identical(d1, d2)  # TRUE
#'
#' @export
causalsim_draw <- function(dgp, seed = NULL) {
  if (!inherits(dgp, "causalsim_dgp")) {
    stop("`dgp` must be a causalsim_dgp object.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  n      <- dgp$n
  cov_df <- .generate_covariate_data(dgp$covar_spec, n)

  p   <- .apply_effect(dgp$propensity_fn, cov_df, n = n)
  a   <- stats::rbinom(n, size = 1L, prob = p)
  mu  <- .apply_effect(dgp$baseline_fn,  cov_df, n = n)
  tau <- .apply_effect(dgp$effect_fn,    cov_df, n = n)
  eps <- stats::rnorm(n, mean = 0, sd = dgp$sigma)
  y   <- mu + tau * a + eps

  if (length(dgp$covar_spec) > 0L) {
    data.frame(cov_df, A = a, Y = y, .tau = tau, .p = p,
               check.names = FALSE)
  } else {
    data.frame(A = a, Y = y, .tau = tau, .p = p)
  }
}
