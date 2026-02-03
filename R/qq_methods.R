#' Print Method for qqfit Objects
#'
#' @param x A \code{qqfit} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.qqfit <- function(x, ...) {
  cat("\nQuantile on Quantiles Regression\n")
  cat("================================\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nObservations:", x$nobs, "\n")
  cat("Groups:", x$ngroups, "\n")
  cat("Quantile grid:", length(x$taus), "points from",
      min(x$taus), "to", max(x$taus), "\n")

  if (x$collapsed) {
    cat("\nNote: Data collapsed to group level (no within-group variation)\n")
  }

  cat("\nVariable classification:\n")
  cat("  Within-group:  ", if (length(x$var_types$within) > 0)
    paste(x$var_types$within, collapse = ", ") else "(none)", "\n")
  cat("  Between-group: ", if (length(x$var_types$between) > 0)
    paste(x$var_types$between, collapse = ", ") else "(none)", "\n")
  cat("  Mixed:         ", if (length(x$var_types$mixed) > 0)
    paste(x$var_types$mixed, collapse = ", ") else "(none)", "\n")

  cat("\nCoefficients:", nrow(x$coef), "parameters\n")

  if (!is.null(x$se)) {
    cat("Bootstrap SE: available\n")
  } else {
    cat("Bootstrap SE: not computed\n")
  }

  invisible(x)
}


#' Summary Method for qqfit Objects
#'
#' @param object A \code{qqfit} object.
#' @param which Character vector of coefficient names to summarize.
#'   If NULL, summarizes all coefficients.
#' @param taus_u Optional subset of u quantiles to display.
#' @param taus_v Optional subset of v quantiles to display.
#' @param ... Additional arguments (ignored).
#'
#' @return A list with summary information (invisibly).
#'
#' @export
summary.qqfit <- function(object, which = NULL, taus_u = NULL, taus_v = NULL, ...) {

  # Get coefficient names
  coef_names <- dimnames(object$coef)[[1]]

  if (is.null(which)) {
    which <- coef_names
  } else {
    # Validate coefficient names
    invalid <- setdiff(which, coef_names)
    if (length(invalid) > 0) {
      stop("Invalid coefficient name(s): ", paste(invalid, collapse = ", "),
           "\nAvailable: ", paste(coef_names, collapse = ", "))
    }
  }

  # Subset taus if requested
  if (is.null(taus_u)) taus_u <- object$taus
  if (is.null(taus_v)) taus_v <- object$taus

  u_idx <- which(object$taus %in% taus_u)
  v_idx <- which(object$taus %in% taus_v)

  cat("\nQuantile on Quantiles Regression Summary\n")
  cat("========================================\n\n")

  cat("Call:\n")
  print(object$call)

  cat("\nObservations:", object$nobs, "\n")
  cat("Groups:", object$ngroups, "\n")

  for (coef_name in which) {
    cat("\n---\nCoefficient:", coef_name, "\n---\n\n")

    # Get coefficient matrix for this variable
    coef_idx <- which(coef_names == coef_name)
    coef_mat <- object$coef[coef_idx, u_idx, v_idx, drop = FALSE]
    coef_mat <- matrix(coef_mat, nrow = length(u_idx), ncol = length(v_idx))
    rownames(coef_mat) <- paste0("u=", taus_u)
    colnames(coef_mat) <- paste0("v=", taus_v)

    cat("Point Estimates:\n")
    print(round(coef_mat, 4))

    if (!is.null(object$se)) {
      se_mat <- object$se[coef_idx, u_idx, v_idx, drop = FALSE]
      se_mat <- matrix(se_mat, nrow = length(u_idx), ncol = length(v_idx))
      rownames(se_mat) <- paste0("u=", taus_u)
      colnames(se_mat) <- paste0("v=", taus_v)

      cat("\nBootstrap Standard Errors:\n")
      print(round(se_mat, 4))

      # t-statistics
      t_mat <- coef_mat / se_mat
      cat("\nt-statistics:\n")
      print(round(t_mat, 2))
    }

    cat("\n")
  }

  invisible(list(
    coef = object$coef[which, u_idx, v_idx, drop = FALSE],
    se = if (!is.null(object$se)) object$se[which, u_idx, v_idx, drop = FALSE] else NULL
  ))
}


#' Extract Coefficients from qqfit Objects
#'
#' @param object A \code{qqfit} object.
#' @param which Character vector of coefficient names to extract.
#'   If NULL, extracts all coefficients.
#' @param ... Additional arguments (ignored).
#'
#' @return An array of coefficients. If \code{which} specifies a single
#'   coefficient, returns a matrix (u x v). Otherwise returns the full
#'   array or subset.
#'
#' @export
coef.qqfit <- function(object, which = NULL, ...) {
  if (is.null(which)) {
    return(object$coef)
  }

  # Validate coefficient names
  coef_names <- dimnames(object$coef)[[1]]
  invalid <- setdiff(which, coef_names)
  if (length(invalid) > 0) {
    stop("Invalid coefficient name(s): ", paste(invalid, collapse = ", "),
         "\nAvailable: ", paste(coef_names, collapse = ", "))
  }

  # Subset
  coef_idx <- match(which, coef_names)

  if (length(which) == 1) {
    # Return matrix for single coefficient
    mat <- object$coef[coef_idx, , ]
    return(mat)
  } else {
    return(object$coef[coef_idx, , , drop = FALSE])
  }
}


#' Extract Standard Errors from qqfit Objects
#'
#' @param object A \code{qqfit} object.
#' @param which Character vector of coefficient names to extract.
#'   If NULL, extracts all.
#'
#' @return An array of bootstrap standard errors, or NULL if not computed.
#'
#' @export
se.qqfit <- function(object, which = NULL) {
  if (is.null(object$se)) {
    message("No standard errors available. Run qq_fit with bootstrap = TRUE.")
    return(NULL)
  }

  if (is.null(which)) {
    return(object$se)
  }

  # Validate coefficient names
  coef_names <- dimnames(object$se)[[1]]
  invalid <- setdiff(which, coef_names)
  if (length(invalid) > 0) {
    stop("Invalid coefficient name(s): ", paste(invalid, collapse = ", "),
         "\nAvailable: ", paste(coef_names, collapse = ", "))
  }

  coef_idx <- match(which, coef_names)

  if (length(which) == 1) {
    return(object$se[coef_idx, , ])
  } else {
    return(object$se[coef_idx, , , drop = FALSE])
  }
}


#' List Available Coefficients
#'
#' @param object A \code{qqfit} object.
#'
#' @return Character vector of coefficient names.
#'
#' @export
coef_names <- function(object) {
  UseMethod("coef_names")
}

#' @export
coef_names.qqfit <- function(object) {
  dimnames(object$coef)[[1]]
}
