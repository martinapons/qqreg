#' First Stage QQ Regression
#'
#' Runs quantile regression within each group for each quantile in the grid.
#' Handles variable selection per group (for mixed variables) and ensures
#' monotonicity by sorting fitted values within groups.
#'
#' @param y Outcome vector.
#' @param X1 Design matrix for first stage (within-group varying variables).
#' @param group_vec Group indicator vector.
#' @param taus Numeric vector of quantile levels.
#' @param var_types List from \code{detect_var_types} with variable classifications.
#' @param parallel Logical, whether to parallelize across groups.
#' @param ncores Number of cores to use if parallel = TRUE.
#'
#' @return A matrix of fitted values (n x n_tau), with rows sorted within groups
#'   to ensure monotonicity in u.
#'
#' @keywords internal
qq_first_stage <- function(y, X1, group_vec, taus, var_types,
                           parallel = TRUE, ncores = NULL) {
  groups <- unique(group_vec)
  n <- length(y)
  n_tau <- length(taus)

  # Function to run first stage for a single group
  run_group <- function(g) {
    idx <- which(group_vec == g)
    y_g <- y[idx]
    n_g <- length(y_g)
    fitted_g <- matrix(NA, n_g, n_tau)

    # Get design matrix for this group
    if (is.null(X1) || ncol(X1) == 0) {
      X1_g <- NULL
    } else {
      X1_g <- X1[idx, , drop = FALSE]

      # Remove zero-variance columns for this group
      var_g <- apply(X1_g, 2, stats::var, na.rm = TRUE)
      var_g[is.na(var_g)] <- 0
      keep_cols <- var_g > 1e-10

      if (any(keep_cols)) {
        X1_g <- X1_g[, keep_cols, drop = FALSE]

        # Check rank and remove linearly dependent columns
        qr_X <- qr(X1_g, tol = 1e-9, LAPACK = FALSE)
        if (qr_X$rank < ncol(X1_g)) {
          X1_g <- X1_g[, qr_X$pivot[seq_len(qr_X$rank)], drop = FALSE]
        }
      } else {
        X1_g <- NULL
      }
    }

    # Run quantile regression for each tau
    for (t_idx in seq_along(taus)) {
      if (is.null(X1_g) || ncol(X1_g) == 0) {
        # No within-group covariates: compute sample quantile
        fitted_g[, t_idx] <- stats::quantile(y_g, taus[t_idx], na.rm = TRUE)
      } else {
        # Run QR with intercept
        qr_fit <- quantreg::rq(y_g ~ X1_g, tau = taus[t_idx])
        fitted_g[, t_idx] <- qr_fit$fitted.values
      }
    }

    # Sort fitted values within group to ensure monotonicity in u
    # Each row: sort across tau columns
    fitted_g <- t(apply(fitted_g, 1, sort))

    list(idx = idx, fitted = fitted_g)
  }

  # Run across all groups
  if (parallel && length(groups) > 1) {
    # Set up parallel backend
    if (is.null(ncores)) {
      ncores <- max(1, floor(0.8 * future::availableCores()))
    }
    future::plan(future::multisession, workers = ncores)

    results <- furrr::future_map(groups, run_group,
                                  .options = furrr::furrr_options(seed = TRUE))
  } else {
    results <- lapply(groups, run_group)
  }

  # Combine results into single matrix
  fitted <- matrix(NA, n, n_tau)
  for (res in results) {
    fitted[res$idx, ] <- res$fitted
  }

  colnames(fitted) <- paste0("u_", taus)
  return(fitted)
}


#' First Stage for Collapsed Data
#'
#' When no within-group varying variables exist, computes group-level quantiles
#' and returns collapsed data with one observation per group per tau.
#'
#' @param y Outcome vector.
#' @param group_vec Group indicator vector.
#' @param taus Numeric vector of quantile levels.
#'
#' @return A list with:
#'   \item{fitted}{Matrix of group-level quantiles (n_groups x n_tau)}
#'   \item{weights}{Vector of group sizes for weighting}
#'   \item{group_order}{Order of groups in the output}
#'
#' @keywords internal
qq_first_stage_collapsed <- function(y, group_vec, taus) {
  groups <- unique(group_vec)
  n_groups <- length(groups)
  n_tau <- length(taus)

  fitted <- matrix(NA, n_groups, n_tau)
  weights <- numeric(n_groups)

  for (i in seq_along(groups)) {
    g <- groups[i]
    idx <- which(group_vec == g)
    y_g <- y[idx]
    weights[i] <- length(y_g)

    for (t_idx in seq_along(taus)) {
      fitted[i, t_idx] <- stats::quantile(y_g, taus[t_idx], na.rm = TRUE)
    }
  }

  colnames(fitted) <- paste0("u_", taus)
  rownames(fitted) <- as.character(groups)

  list(
    fitted = fitted,
    weights = weights,
    group_order = groups
  )
}
