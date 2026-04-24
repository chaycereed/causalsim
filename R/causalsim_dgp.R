#' Create a causal data generating process
#'
#' Defines a causal DGP with known ground truth. Covariates can be specified
#' via shorthand count arguments (Option B), an explicit named list of
#' [covar()] objects (Option A), or both combined.
#'
#' @param n Positive integer. Sample size for each simulated dataset.
#' @param effect Numeric scalar or function. A scalar specifies a constant
#'   (homogeneous) treatment effect — ATE = CATE everywhere. A function should
#'   accept named arguments matching covariate names defined in the DGP and
#'   return a numeric vector of individual-level causal effects (CATE).
#'   See Details.
#' @param propensity Numeric scalar, preset string, or function. Treatment
#'   assignment probability. A scalar (e.g. `0.5`) gives a constant propensity
#'   (randomized trial). Preset strings `"low"`, `"moderate"`, `"high"` generate
#'   a logistic propensity over confounders with coefficients 0.25, 0.5, and
#'   1.0 respectively. A function follows the same named-argument convention as
#'   `effect` and must return values in \[0, 1\]. Defaults to `"moderate"`.
#' @param baseline Numeric scalar, preset string, or function. Mean potential
#'   outcome under control, `E[Y(0) | W]`. Preset strings follow the same
#'   levels as `propensity` and apply a linear combination of confounders.
#'   Defaults to `0`.
#' @param sigma Positive numeric. Standard deviation of the outcome noise term.
#'   Default `1`.
#' @param covariates Named list of [covar()] objects (Option A / explicit path).
#'   Each name becomes the column name in generated data and the argument name
#'   expected by `effect`, `propensity`, and `baseline` functions. Merged with
#'   any auto-generated covariates; name collisions error.
#' @param n_confounders Non-negative integer. Standard normal confounders
#'   auto-generated as `W` (single) or `W1, W2, ...` (multiple).
#' @param n_effect_modifiers Non-negative integer. Auto-generates standard
#'   normal effect modifiers as `V` or `V1, V2, ...`.
#' @param n_instruments Non-negative integer. Auto-generates standard normal
#'   instruments as `Z` or `Z1, Z2, ...`.
#' @param n_noise Non-negative integer. Auto-generates standard normal noise
#'   covariates as `X` or `X1, X2, ...`.
#' @param mc_draws Positive integer. Monte Carlo draws for true ATE
#'   approximation. Default `10000L`. Ignored for scalar `effect`.
#'
#' @details
#' ## Structural model
#'
#' ```
#' W ~ covariate_spec
#' A ~ Bernoulli(propensity(W))
#' Y  = baseline(W) + effect(W) * A + N(0, sigma^2)
#' ```
#'
#' `A` is used for treatment throughout to avoid collision with R's built-in
#' `T` alias.
#'
#' ## Function calling convention
#'
#' The `effect`, `propensity`, and `baseline` functions are called with named
#' arguments matching covariate names — not a data frame. Write:
#'
#' ```r
#' effect = function(W) 2 + 1.5 * W
#' propensity = function(W1, W2) plogis(0.3 * W1 + 0.5 * W2)
#' ```
#'
#' Every argument name is validated against the DGP's covariate spec at
#' construction time, so mismatches surface immediately rather than at draw
#' time. The propensity function is additionally evaluated on a small test draw
#' to confirm it returns values in \[0, 1\].
#'
#' ## True ATE
#'
#' For scalar `effect` the true ATE is exact. For function `effect` it is
#' approximated via Monte Carlo over `mc_draws` draws from the covariate
#' distribution.
#'
#' @return An S3 object of class `causalsim_dgp` with components:
#' \describe{
#'   \item{`n`}{Sample size (integer)}
#'   \item{`covar_spec`}{Named list of [covar()] objects}
#'   \item{`effect_fn`}{Normalized effect function}
#'   \item{`propensity_fn`}{Normalized propensity function}
#'   \item{`baseline_fn`}{Normalized baseline function}
#'   \item{`sigma`}{Outcome noise standard deviation}
#'   \item{`true_ate`}{True ATE — exact for scalar, MC approximation otherwise}
#'   \item{`heterogeneous`}{Logical; `TRUE` if `effect` was a function}
#'   \item{`mc_draws`}{Monte Carlo draws used (integer)}
#' }
#'
#' @examples
#' # Minimal: one confounder, constant effect, moderate confounding
#' dgp <- causalsim_dgp(n = 500, n_confounders = 1, effect = 2)
#' dgp
#'
#' # Heterogeneous effect, explicit covariate spec
#' dgp2 <- causalsim_dgp(
#'   n = 500,
#'   covariates = list(
#'     W = covar("normal", role = "confounder"),
#'     V = covar("binary", role = "effect_modifier", prob = 0.4)
#'   ),
#'   effect     = function(V) 2 + 1.5 * V,
#'   propensity = function(W) plogis(0.5 * W),
#'   baseline   = function(W) 1.5 * W
#' )
#'
#' # Mixed: shorthand confounders + explicit instrument + RCT propensity
#' dgp3 <- causalsim_dgp(
#'   n = 1000,
#'   n_confounders = 2,
#'   covariates    = list(Z = covar("normal", role = "instrument")),
#'   effect        = 1,
#'   propensity    = 0.5
#' )
#'
#' @export
causalsim_dgp <- function(
  n,
  effect             = 1,
  propensity         = "moderate",
  baseline           = 0,
  sigma              = 1,
  covariates         = list(),
  n_confounders      = 0L,
  n_effect_modifiers = 0L,
  n_instruments      = 0L,
  n_noise            = 0L,
  mc_draws           = 10000L
) {
  n        <- .validate_positive(n,        "n",        as_int = TRUE)
  mc_draws <- .validate_positive(mc_draws, "mc_draws", as_int = TRUE)
  sigma    <- .validate_positive(sigma,    "sigma")
  n_confounders      <- .validate_nonneg(
    n_confounders, "n_confounders", as_int = TRUE
  )
  n_effect_modifiers <- .validate_nonneg(
    n_effect_modifiers, "n_effect_modifiers", as_int = TRUE
  )
  n_instruments <- .validate_nonneg(
    n_instruments, "n_instruments", as_int = TRUE
  )
  n_noise <- .validate_nonneg(n_noise, "n_noise", as_int = TRUE)

  covar_spec <- .build_covariate_spec(
    n_confounders, n_effect_modifiers, n_instruments, n_noise, covariates
  )

  heterogeneous  <- is.function(effect)
  effect_fn      <- .normalize_fn(effect,     "effect",     covar_spec)
  propensity_fn  <- .normalize_fn(propensity, "propensity", covar_spec)
  baseline_fn    <- .normalize_fn(baseline,   "baseline",   covar_spec)

  .validate_fn_args(effect_fn,     names(covar_spec), "effect")
  .validate_fn_args(propensity_fn, names(covar_spec), "propensity")
  .validate_fn_args(baseline_fn,   names(covar_spec), "baseline")
  .validate_propensity_fn(propensity_fn, covar_spec)

  true_ate <- .mc_ate(effect_fn, covar_spec, mc_draws, heterogeneous)

  structure(
    list(
      n             = n,
      covar_spec    = covar_spec,
      effect_fn     = effect_fn,
      propensity_fn = propensity_fn,
      baseline_fn   = baseline_fn,
      sigma         = sigma,
      true_ate      = true_ate,
      heterogeneous = heterogeneous,
      mc_draws      = mc_draws,
      params        = list(
        n                  = n,
        effect             = effect,
        propensity         = propensity,
        baseline           = baseline,
        sigma              = sigma,
        covariates         = covariates,
        n_confounders      = n_confounders,
        n_effect_modifiers = n_effect_modifiers,
        n_instruments      = n_instruments,
        n_noise            = n_noise,
        mc_draws           = mc_draws
      )
    ),
    class = "causalsim_dgp"
  )
}

#' @export
print.causalsim_dgp <- function(x, ...) {
  cat("<causalsim_dgp>\n")
  cat(sprintf("  n            : %d\n",   x$n))
  cat(sprintf("  true ATE     : %.4f\n", x$true_ate))
  cat(sprintf("  heterogeneous: %s\n",   x$heterogeneous))
  cat(sprintf("  sigma        : %.2f\n", x$sigma))
  if (length(x$covar_spec) > 0L) {
    cat("  covariates   :\n")
    for (nm in names(x$covar_spec)) {
      cv <- x$covar_spec[[nm]]
      cat(sprintf("    %-6s  %s  [%s]\n",
                  nm, cv$dist, paste(cv$role, collapse = "+")))
    }
  }
  invisible(x)
}
