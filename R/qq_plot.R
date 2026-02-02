#' Plot Method for qqfit Objects
#'
#' Creates visualizations of the QQ regression coefficient surface.
#'
#' @param x A \code{qqfit} object.
#' @param which Character, coefficient name to plot. Use \code{"ALL"} to
#'   generate plots for all coefficients (returns a list). Default is the
#'   first non-intercept coefficient.
#' @param type Character, type of plot:
#'   \itemize{
#'     \item \code{"surface"}: 2D contour/tile plot of coefficient surface
#'     \item \code{"slice_u"}: Fix u, plot coefficients across v (multiple panels)
#'     \item \code{"slice_v"}: Fix v, plot coefficients across u (multiple panels)
#'   }
#' @param taus_slice Numeric vector, quantiles to use for slice plots.
#'   Default is \code{seq(0.1, 0.9, 0.1)}. Must be a subset of the model's
#'   quantile grid.
#' @param ncol Integer, number of columns for slice plot panels (default 3).
#' @param ylim Numeric vector of length 2, y-axis limits for slice plots.
#' @param ... Additional arguments passed to underlying plot functions.
#'
#' @return A ggplot object (or list of ggplot objects if \code{which = "ALL"}).
#'
#' @details
#' Requires the \pkg{ggplot2} package. For surface plots, \pkg{metR} is used
#' for filled contours if available, otherwise \code{geom_tile} is used.
#' For multi-panel slice plots, \pkg{patchwork} is required.
#'
#' @examples
#' \dontrun{
#' fit <- qq_fit(y ~ x1 + x2, data = data, group = "group")
#'
#' # Surface plot for one coefficient
#' plot(fit, which = "x1", type = "surface")
#'
#' # Slice plots fixing u
#' plot(fit, which = "x1", type = "slice_u")
#'
#' # All coefficients
#' plots <- plot(fit, which = "ALL", type = "surface")
#' }
#'
#' @export
plot.qqfit <- function(x, which = NULL, type = c("surface", "slice_u", "slice_v"),
                        taus_slice = seq(0.1, 0.9, 0.1), ncol = 3,
                        ylim = NULL, ...) {

  # Check for ggplot2

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  type <- match.arg(type)

  # Get coefficient names
  coef_names <- dimnames(x$coef)[[1]]

  # Handle which argument
  if (is.null(which)) {
    # Default to first non-intercept coefficient
    non_intercept <- setdiff(coef_names, "(Intercept)")
    if (length(non_intercept) > 0) {
      which <- non_intercept[1]
      message("Plotting coefficient: ", which)
    } else {
      which <- coef_names[1]
      message("Plotting coefficient: ", which)
    }
  }

  # Handle "ALL"
  if (identical(which, "ALL")) {
    plots <- lapply(coef_names, function(w) {
      plot_single(x, which = w, type = type, taus_slice = taus_slice,
                  ncol = ncol, ylim = ylim, ...)
    })
    names(plots) <- coef_names
    return(plots)
  }

  # Validate coefficient name
  if (!which %in% coef_names) {
    stop("Invalid coefficient name: '", which, "'.\n",
         "Available coefficients: ", paste(coef_names, collapse = ", "))
  }

  plot_single(x, which = which, type = type, taus_slice = taus_slice,
              ncol = ncol, ylim = ylim, ...)
}


#' Internal: Plot Single Coefficient
#'
#' @keywords internal
plot_single <- function(x, which, type, taus_slice, ncol, ylim, ...) {

  # Extract coefficient matrix and SE matrix for this variable
  coef_idx <- which(dimnames(x$coef)[[1]] == which)
  coefs <- x$coef[coef_idx, , ]
  se <- if (!is.null(x$se)) x$se[coef_idx, , ] else NULL
  taus <- x$taus

  switch(type,
    "surface" = plot_surface(coefs, se, taus, which, ...),
    "slice_u" = plot_slices_u(coefs, se, taus, taus_slice, which, ncol, ylim, ...),
    "slice_v" = plot_slices_v(coefs, se, taus, taus_slice, which, ncol, ylim, ...)
  )
}


#' Internal: Surface Plot
#'
#' @keywords internal
plot_surface <- function(coefs, se, taus, coef_name, ...) {

  # Create tidy data for plotting
  tau1 <- rep(taus, times = length(taus))
  tau2 <- rep(taus, each = length(taus))
  tidy_mat <- data.frame(
    u = tau1,
    v = tau2,
    effect = as.vector(coefs)
  )

  # Define color palette
  myColors <- c("#11333a", "#0faeb6", "#fda535", "#be3821")

  # Try to use metR for filled contours, fall back to geom_tile
  if (requireNamespace("metR", quietly = TRUE)) {
    p <- ggplot2::ggplot(tidy_mat, ggplot2::aes(x = .data$v, y = .data$u)) +
      metR::geom_contour_fill(ggplot2::aes(z = .data$effect),
                               binwidth = (max(coefs) - min(coefs)) / 30) +
      ggplot2::scale_fill_gradientn(name = "Effect", colours = myColors)
  } else {
    p <- ggplot2::ggplot(tidy_mat, ggplot2::aes(x = .data$v, y = .data$u)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$effect)) +
      ggplot2::scale_fill_gradientn(name = "Effect", colours = myColors)
  }

  p <- p +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 15)) +
    ggplot2::ylab(expression(u)) +
    ggplot2::xlab(expression(v)) +
    ggplot2::ggtitle(coef_name) +
    ggplot2::theme_minimal()

  return(p)
}


#' Internal: Slice Plot Fixing u
#'
#' @keywords internal
plot_slices_u <- function(coefs, se, taus, taus_slice, coef_name, ncol, ylim, ...) {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for slice plots. Please install it.")
  }

  # Check that taus_slice values are in the model's quantile grid
  valid_taus <- round(taus_slice, 3) %in% round(taus, 3)
  if (!all(valid_taus)) {
    missing <- taus_slice[!valid_taus]
    stop("Some elements of taus_slice are not in the model's quantile grid: ",
         paste(missing, collapse = ", "), ".\n",
         "Either change your taus grid in qq_fit() or specify taus_slice explicitly.")
  }

  # Determine which panels show y-axis
  show_y_idx <- seq(1, length(taus_slice), by = ncol)

  myColors <- c("#11333a", "#0faeb6", "#fda535", "#be3821")

  plots <- list()
  for (k in seq_along(taus_slice)) {
    u_val <- taus_slice[k]
    u_idx <- which(abs(taus - u_val) < 1e-6)

    if (length(u_idx) == 0) next

    beta <- coefs[u_idx, ]
    se_vals <- if (!is.null(se)) se[u_idx, ] else rep(NA, length(taus))

    # Confidence interval
    crit <- stats::qnorm(0.975)
    lb <- beta - se_vals * crit
    ub <- beta + se_vals * crit

    dd <- data.frame(v = taus, beta = beta, se = se_vals, lb = lb, ub = ub)

    title <- bquote(paste(u, " = ", .(u_val)))

    p <- ggplot2::ggplot(data = dd) +
      ggplot2::theme_minimal() +
      ggplot2::aes(x = .data$v, y = .data$beta) +
      ggplot2::geom_line(linewidth = 0.5, color = myColors[1]) +
      ggplot2::geom_point(color = myColors[1], size = 1) +
      ggplot2::labs(title = title, x = expression(v)) +
      ggplot2::coord_cartesian(ylim = ylim)

    if (!all(is.na(se_vals))) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lb, ymax = .data$ub),
                                     alpha = 0.2, fill = myColors[2])
    }

    # Only show y-axis label on first column
    if (k %in% show_y_idx) {
      p <- p + ggplot2::ylab("Point Estimate")
    } else {
      p <- p + ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
    }

    plots[[k]] <- p
  }

  # Combine with patchwork
  combined <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(title = coef_name)

  return(combined)
}


#' Internal: Slice Plot Fixing v
#'
#' @keywords internal
plot_slices_v <- function(coefs, se, taus, taus_slice, coef_name, ncol, ylim, ...) {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for slice plots. Please install it.")
  }

  # Check that taus_slice values are in the model's quantile grid
  valid_taus <- round(taus_slice, 3) %in% round(taus, 3)
  if (!all(valid_taus)) {
    missing <- taus_slice[!valid_taus]
    stop("Some elements of taus_slice are not in the model's quantile grid: ",
         paste(missing, collapse = ", "), ".\n",
         "Either change your taus grid in qq_fit() or specify taus_slice explicitly.")
  }

  show_y_idx <- seq(1, length(taus_slice), by = ncol)
  myColors <- c("#11333a", "#0faeb6", "#fda535", "#be3821")

  plots <- list()
  for (k in seq_along(taus_slice)) {
    v_val <- taus_slice[k]
    v_idx <- which(abs(taus - v_val) < 1e-6)

    if (length(v_idx) == 0) next

    beta <- coefs[, v_idx]
    se_vals <- if (!is.null(se)) se[, v_idx] else rep(NA, length(taus))

    crit <- stats::qnorm(0.975)
    lb <- beta - se_vals * crit
    ub <- beta + se_vals * crit

    dd <- data.frame(u = taus, beta = beta, se = se_vals, lb = lb, ub = ub)

    title <- bquote(paste(v, " = ", .(v_val)))

    p <- ggplot2::ggplot(data = dd) +
      ggplot2::theme_minimal() +
      ggplot2::aes(x = .data$u, y = .data$beta) +
      ggplot2::geom_line(linewidth = 0.5, color = myColors[1]) +
      ggplot2::geom_point(color = myColors[1], size = 1) +
      ggplot2::labs(title = title, x = expression(u)) +
      ggplot2::coord_cartesian(ylim = ylim)

    if (!all(is.na(se_vals))) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lb, ymax = .data$ub),
                                     alpha = 0.2, fill = myColors[2])
    }

    if (k %in% show_y_idx) {
      p <- p + ggplot2::ylab("Point Estimate")
    } else {
      p <- p + ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
    }

    plots[[k]] <- p
  }

  combined <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(title = coef_name)

  return(combined)
}
