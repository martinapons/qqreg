#' Second Stage QQ Regression
#'
#' Runs quantile regression of first-stage fitted values on the full
#' design matrix for each (u, v) pair in the quantile grid.
#'
#' @param fitted_first Matrix of first-stage fitted values (n x n_tau).
#' @param X Full design matrix including intercept and all covariates.
#' @param taus Numeric vector of quantile levels.
#' @param weights Optional observation weights.
#'
#' @return A list with:
#'   \item{coef}{Array of coefficients (k x n_tau x n_tau)}
#'   \item{fitted}{Array of fitted values (n x n_tau x n_tau)}
#'   \item{residuals}{Array of residuals (n x n_tau x n_tau)}
#'
#' @keywords internal
qq_second_stage <- function(fitted_first, X, taus, weights = NULL) {
  n_tau <- length(taus)
  n <- nrow(X)
  k <- ncol(X)

  # Initialize arrays
  coef_array <- array(NA, dim = c(k, n_tau, n_tau))
  fitted_array <- array(NA, dim = c(n, n_tau, n_tau))
  resid_array <- array(NA, dim = c(n, n_tau, n_tau))

  # Loop over u (rows of first-stage fitted values)
  for (u_idx in seq_along(taus)) {
    y_u <- fitted_first[, u_idx]

    # Loop over v
    for (v_idx in seq_along(taus)) {
      # Run quantile regression
      if (is.null(weights)) {
        qr_fit <- quantreg::rq(y_u ~ X - 1, tau = taus[v_idx])
      } else {
        qr_fit <- quantreg::rq(y_u ~ X - 1, tau = taus[v_idx], weights = weights)
      }

      coef_array[, u_idx, v_idx] <- qr_fit$coefficients
      fitted_array[, u_idx, v_idx] <- qr_fit$fitted.values
      resid_array[, u_idx, v_idx] <- qr_fit$residuals
    }
  }

  # Set dimension names
  dimnames(coef_array) <- list(
    colnames(X),
    paste0("u_", taus),
    paste0("v_", taus)
  )
  dimnames(fitted_array) <- list(
    NULL,
    paste0("u_", taus),
    paste0("v_", taus)
  )
  dimnames(resid_array) <- dimnames(fitted_array)

  list(
    coef = coef_array,
    fitted = fitted_array,
    residuals = resid_array
  )
}


#' Compute Ranks from Second Stage
#'
#' Computes within-group ranks based on second-stage residuals.
#' For each observation and each u, the rank indicates the position
#' in the within-group distribution.
#'
#' @param residuals Array of second-stage residuals (n x n_tau x n_tau).
#' @param group_vec Group indicator vector.
#' @param taus Numeric vector of quantile levels.
#'
#' @return A matrix of ranks (n x n_tau), where each column corresponds
#'   to a value of u.
#'
#' @keywords internal
compute_ranks <- function(residuals, group_vec, taus) {
  n <- dim(residuals)[1]
  n_tau <- length(taus)

  rank_mat <- matrix(0, n, n_tau)

  # For each u, count how many v have positive residuals
  for (u_idx in seq_along(taus)) {
    for (v_idx in seq_along(taus)) {
      rank_mat[, u_idx] <- rank_mat[, u_idx] + (residuals[, u_idx, v_idx] > 0)
    }
  }

  colnames(rank_mat) <- paste0("u_", taus)
  return(rank_mat)
}
