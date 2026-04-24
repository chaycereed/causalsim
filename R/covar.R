#' Define a covariate for a causal DGP
#'
#' Constructs a fully-specified covariate descriptor for use in
#' [causalsim_dgp()]. Every covariate has a distribution family, optional
#' distribution parameters, and one or more causal roles that control how it
#' enters the automatically-generated treatment and outcome models.
#'
#' @param dist Character. Distribution family. One of `"normal"`, `"binary"`,
#'   `"uniform"`.
#' @param role Character vector. Causal role(s). One or more of:
#'   \describe{
#'     \item{`"confounder"`}{Enters both the propensity model and the outcome
#'       baseline — creates confounding bias in naive estimators.}
#'     \item{`"instrument"`}{Enters the propensity model only (exclusion
#'       restriction holds).}
#'     \item{`"effect_modifier"`}{Available for use in the `effect` function;
#'       excluded from the propensity model unless also a `"confounder"`.}
#'     \item{`"noise"`}{Independent of both treatment and outcome; adds
#'       variance without confounding.}
#'   }
#'   Roles are not mutually exclusive: `role = c("confounder", "effect_modifier")`
#'   specifies a variable that both confounds and moderates the effect.
#' @param ... Distribution parameters.
#'   - `"normal"`: `mean` (default 0), `sd` (default 1)
#'   - `"binary"`: `prob` (default 0.5)
#'   - `"uniform"`: `min` (default 0), `max` (default 1)
#'
#' @return An S3 object of class `causalsim_covar`.
#'
#' @examples
#' # Standard normal confounder
#' covar("normal", role = "confounder", mean = 0, sd = 1)
#'
#' # Binary effect modifier
#' covar("binary", role = "effect_modifier", prob = 0.4)
#'
#' # Variable that both confounds and moderates the effect
#' covar("normal", role = c("confounder", "effect_modifier"))
#'
#' @export
covar <- function(dist = "normal", role = "confounder", ...) {
  dist <- match.arg(dist, c("normal", "binary", "uniform"))
  role <- match.arg(
    role,
    choices  = c("confounder", "instrument", "effect_modifier", "noise"),
    several.ok = TRUE
  )
  structure(
    list(dist = dist, role = role, params = list(...)),
    class = "causalsim_covar"
  )
}

#' @export
print.causalsim_covar <- function(x, ...) {
  params_str <- if (length(x$params) > 0L) {
    paste(names(x$params), x$params, sep = "=", collapse = ", ")
  } else {
    ""
  }
  dist_str <- if (nchar(params_str) > 0L) {
    sprintf("%s(%s)", x$dist, params_str)
  } else {
    x$dist
  }
  cat(sprintf("<causalsim_covar>  dist: %-20s  role: [%s]\n",
              dist_str, paste(x$role, collapse = ", ")))
  invisible(x)
}
