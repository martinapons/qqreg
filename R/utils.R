#' @importFrom stats model.frame model.matrix model.response var quantile qnorm na.pass
#' @importFrom future plan multisession availableCores
#' @importFrom furrr future_map furrr_options
#' @importFrom quantreg rq
NULL

#' qqreg: Quantile-on-Quantiles Regression
#'
#' Implements the two-step quantile regression estimator from
#' Pons (2026) "Quantile on Quantiles". The package estimates a
#' two-dimensional quantile surface capturing heterogeneity both
#' within and between groups.
#'
#' @section Main function:
#' \describe{
#'   \item{\code{\link{qq_fit}}}{Fit a QQ regression model}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{print.qqfit}}{Print model summary}
#'   \item{\code{summary.qqfit}}{Detailed model summary}
#'   \item{\code{coef.qqfit}}{Extract coefficients}
#'   \item{\code{plot.qqfit}}{Plot coefficient surfaces}
#' }
#'
#' @docType package
#' @name qqreg-package
"_PACKAGE"


#' Check if Variable is Constant Within Groups
#'
#' @param x A vector.
#' @param group A group indicator vector.
#' @param tol Tolerance for checking variance.
#'
#' @return Logical, TRUE if x is constant within all groups.
#'
#' @keywords internal
is_between_var <- function(x, group, tol = 1e-10) {
  within_var <- tapply(x, group, stats::var, na.rm = TRUE)
  within_var[is.na(within_var)] <- 0
  all(within_var < tol)
}


#' Check if Variable Varies Within All Groups
#'
#' @param x A vector.
#' @param group A group indicator vector.
#' @param tol Tolerance for checking variance.
#'
#' @return Logical, TRUE if x varies within all groups.
#'
#' @keywords internal
is_within_var <- function(x, group, tol = 1e-10) {
  within_var <- tapply(x, group, stats::var, na.rm = TRUE)
  within_var[is.na(within_var)] <- 0
  all(within_var > tol)
}


#' Format Numbers for Display
#'
#' @param x Numeric vector.
#' @param digits Number of decimal places.
#'
#' @return Character vector.
#'
#' @keywords internal
format_num <- function(x, digits = 4) {
  formatC(x, digits = digits, format = "f")
}


#' Safe Quantile Regression
#'
#' Wrapper around \code{quantreg::rq} that handles errors gracefully.
#'
#' @param y Response vector.
#' @param X Design matrix.
#' @param tau Quantile level.
#' @param weights Optional weights.
#'
#' @return A list with coefficients and fitted values, or NAs on error.
#'
#' @keywords internal
safe_rq <- function(y, X, tau, weights = NULL) {
  tryCatch({
    if (is.null(weights)) {
      fit <- quantreg::rq(y ~ X - 1, tau = tau)
    } else {
      fit <- quantreg::rq(y ~ X - 1, tau = tau, weights = weights)
    }
    list(
      coef = fit$coefficients,
      fitted = fit$fitted.values,
      residuals = fit$residuals,
      success = TRUE
    )
  }, error = function(e) {
    k <- ncol(X)
    n <- length(y)
    list(
      coef = rep(NA_real_, k),
      fitted = rep(NA_real_, n),
      residuals = rep(NA_real_, n),
      success = FALSE
    )
  })
}
