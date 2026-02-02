#' Clustered Bootstrap for QQ Regression
#'
#' Performs clustered bootstrap inference by resampling entire groups
#' with replacement. Parallelization is done over bootstrap iterations
#' (nboot) since nboot >> n_tau.
#'
#' @param fitted_first Matrix of first-stage fitted values (n x n_tau).
#' @param X Full design matrix.
#' @param group_vec Group indicator vector.
#' @param taus Numeric vector of quantile levels.
#' @param nboot Number of bootstrap replications.
#' @param parallel Logical, whether to parallelize.
#' @param ncores Number of cores to use if parallel = TRUE.
#' @param weights Optional observation weights.
#'
#' @return A list with:
#'   \item{se}{Array of bootstrap standard errors (k x n_tau x n_tau)}
#'   \item{coef_boot}{Array of bootstrap coefficient draws (k x n_tau x n_tau x nboot)}
#'
#' @keywords internal
qq_bootstrap <- function(fitted_first, X, group_vec, taus, nboot,
                          parallel = TRUE, ncores = NULL,
                          weights = NULL) {
  groups <- unique(group_vec)
  n_groups <- length(groups)
  k <- ncol(X)
  n_tau <- length(taus)

  # Pre-compute group indices for efficiency
  group_indices <- split(seq_along(group_vec), group_vec)

  # Function to run one bootstrap iteration
  run_boot <- function(b) {
    # Resample groups with replacement
    boot_groups <- sample(groups, n_groups, replace = TRUE)

    # Build bootstrap sample indices
    boot_idx <- unlist(group_indices[as.character(boot_groups)], use.names = FALSE)

    # Get bootstrap data
    fitted_boot <- fitted_first[boot_idx, , drop = FALSE]
    X_boot <- X[boot_idx, , drop = FALSE]
    weights_boot <- if (!is.null(weights)) weights[boot_idx] else NULL

    # Run second stage on bootstrap sample (sequential over tau)
    coef_b <- array(NA, dim = c(k, n_tau, n_tau))

    for (u_idx in seq_along(taus)) {
      y_u <- fitted_boot[, u_idx]

      for (v_idx in seq_along(taus)) {
        qr_fit <- tryCatch({
          if (is.null(weights_boot)) {
            quantreg::rq(y_u ~ X_boot - 1, tau = taus[v_idx])
          } else {
            quantreg::rq(y_u ~ X_boot - 1, tau = taus[v_idx], weights = weights_boot)
          }
        }, error = function(e) NULL)

        if (!is.null(qr_fit)) {
          coef_b[, u_idx, v_idx] <- qr_fit$coefficients
        }
      }
    }

    return(coef_b)
  }

  # Run bootstrap iterations
  if (parallel) {
    # Set up parallel backend
    if (is.null(ncores)) {
      ncores <- max(1, floor(0.8 * future::availableCores()))
    }
    future::plan(future::multisession, workers = ncores)

    boot_results <- furrr::future_map(
      seq_len(nboot),
      run_boot,
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
  } else {
    boot_results <- lapply(seq_len(nboot), run_boot)
  }

  # Combine into 4D array: k x n_tau x n_tau x nboot
  boot_coefs <- array(NA, dim = c(k, n_tau, n_tau, nboot))
  for (b in seq_along(boot_results)) {
    boot_coefs[, , , b] <- boot_results[[b]]
  }

  # Set dimension names
  dimnames(boot_coefs) <- list(
    colnames(X),
    paste0("u_", taus),
    paste0("v_", taus),
    paste0("boot_", seq_len(nboot))
  )

  # Compute SE from bootstrap draws using IQR / 1.34 (robust to outliers)
  se_array <- apply(boot_coefs, 1:3, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 3) return(NA)
    (stats::quantile(x, 0.75) - stats::quantile(x, 0.25)) / 1.34
  })

  dimnames(se_array) <- list(
    colnames(X),
    paste0("u_", taus),
    paste0("v_", taus)
  )

  list(
    se = se_array,
    coef_boot = boot_coefs
  )
}
