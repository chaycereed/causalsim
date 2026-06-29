#' Simulate a causal dataset with known ground truth
#'
#' A convenience wrapper around [causalsim_dgp()] and [causalsim_draw()].
#' Defines a causal data generating process and returns one simulated dataset
#' in a single call.
#'
#' @inheritParams causalsim_dgp
#' @param seed Integer or `NULL`. Passed to [set.seed()] before drawing.
#'   Default `NULL`.
#'
#' @details
#' The returned data frame includes all covariate columns, the treatment
#' indicator `A`, the outcome `Y`, and two ground-truth columns:
#'
#' \describe{
#'   \item{`.tau`}{Individual-level causal effect (CATE).}
#'   \item{`.p`}{True propensity score.}
#' }
#'
#' For full control over the DGP (multiple draws, grid evaluation, estimator
#' benchmarking), use [causalsim_dgp()] and [causalsim_draw()] directly.
#'
#' @return A data frame with covariate columns, `A`, `Y`, `.tau`, and `.p`.
#'
#' @examples
#' # One confounder, constant effect of 2, moderate confounding
#' data <- causalsim(n = 500, n_confounders = 1, effect = 2, seed = 1L)
#' head(data)
#'
#' # Heterogeneous effect with an explicit covariate spec
#' data2 <- causalsim(
#'   n = 500,
#'   covariates = list(
#'     W = covar("normal", role = "confounder"),
#'     V = covar("binary", role = "effect_modifier", prob = 0.4)
#'   ),
#'   effect = function(V) 2 + 1.5 * V,
#'   propensity = function(W) plogis(0.5 * W),
#'   seed = 42L
#' )
#' head(data2)
#'
#' @export
causalsim <- function(
  n,
  effect = 1,
  propensity = "moderate",
  baseline = 0,
  sigma = 1,
  covariates = list(),
  n_confounders = 0L,
  n_effect_modifiers = 0L,
  n_instruments = 0L,
  n_noise = 0L,
  mc_draws = 10000L,
  seed = NULL
) {
  dgp <- causalsim_dgp(
    n = n,
    effect = effect,
    propensity = propensity,
    baseline = baseline,
    sigma = sigma,
    covariates = covariates,
    n_confounders = n_confounders,
    n_effect_modifiers = n_effect_modifiers,
    n_instruments = n_instruments,
    n_noise = n_noise,
    mc_draws = mc_draws
  )
  if (!is.null(seed)) set.seed(seed)
  causalsim_draw(dgp)
}
