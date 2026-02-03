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
#' @param ncol Integer, number of columns for slice plot panels. If NULL (default),
#'   automatically selects: 1 slice uses 1 column, 2/4/8 slices use 2 columns,
#'   3/6/9 or more slices use 3 columns.
#' @param shared_ylim Logical, whether all slices should share the same y-axis
#'   scale (default FALSE). If TRUE, a global y-axis range is computed across
#'   all slices and only the first column displays y-axis labels. If FALSE,
#'   each slice has its own y-axis scale and all y-axis labels are displayed.
#' @param ylim Numeric vector of length 2, y-axis limits for slice plots.
#' @param ylab Character, y-axis label. Default NULL uses \code{"Point Estimate"}
#'   for slice plots and \code{expression(u)} for surface plots.
#' @param xlab Character, x-axis label. Default NULL uses \code{expression(v)}
#'   for \code{slice_u} and surface, \code{expression(u)} for \code{slice_v}.
#' @param ggtitle Character, overall figure title. Default NULL uses the
#'   coefficient name.
#' @param panel_prefix Character, prefix for per-panel titles in slice plots.
#'   Title is built as \code{paste0(panel_prefix, value)}. Default NULL uses
#'   \code{"u = "} for \code{slice_u} and \code{"v = "} for \code{slice_v}.
#'   Set to \code{""} to suppress the prefix and show only the value.
#' @param line_color Character, colour for the line and points (default \code{"#11333a"}).
#' @param ribbon_fill Character, fill colour for the confidence ribbon (default \code{"#0faeb6"}).
#' @param point_color Character, colour for the points. Default NULL uses the
#'   same colour as \code{line_color}.
#' @param ribbon_alpha Numeric, transparency of the confidence ribbon (default 0.2).
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
                        taus_slice = seq(0.1, 0.9, 0.1), ncol = NULL,
                        shared_ylim = FALSE, ylim = NULL,
                        ylab = NULL, xlab = NULL, ggtitle = NULL, panel_prefix = NULL,
                        line_color = "#11333a", ribbon_fill = "#0faeb6",
                        point_color = NULL, ribbon_alpha = 0.2, ...) {

  # Check for ggplot2

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  type <- match.arg(type)

  # Automatic ncol selection for slice plots
  if (is.null(ncol)) {
    if (type %in% c("slice_u", "slice_v")) {
      n_slices <- length(taus_slice)
      # 1 column for 1 slice, 2 columns for 2/4/8 slices, 3 columns for 3/6/9+ slices
      if (n_slices == 1) {
        ncol <- 1
      } else if (n_slices %in% c(2, 4, 8)) {
        ncol <- 2
      } else {
        ncol <- 3
      }
    } else {
      ncol <- 3  # Default for surface plots (not used but safe)
    }
  }

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
                  ncol = ncol, shared_ylim = shared_ylim, ylim = ylim,
                  ylab = ylab, xlab = xlab, ggtitle = ggtitle, panel_prefix = panel_prefix,
                  line_color = line_color, ribbon_fill = ribbon_fill,
                  point_color = point_color, ribbon_alpha = ribbon_alpha, ...)
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
              ncol = ncol, shared_ylim = shared_ylim, ylim = ylim,
              ylab = ylab, xlab = xlab, ggtitle = ggtitle, panel_prefix = panel_prefix,
              line_color = line_color, ribbon_fill = ribbon_fill,
              point_color = point_color, ribbon_alpha = ribbon_alpha, ...)
}


#' Internal: Plot Single Coefficient
#'
#' @keywords internal
plot_single <- function(x, which, type, taus_slice, ncol, shared_ylim, ylim,
                        ylab, xlab, ggtitle, panel_prefix,
                        line_color, ribbon_fill, point_color, ribbon_alpha, ...) {

  # Extract coefficient matrix and SE matrix for this variable
  coef_idx <- which(dimnames(x$coef)[[1]] == which)
  coefs <- x$coef[coef_idx, , ]
  se <- if (!is.null(x$se)) x$se[coef_idx, , ] else NULL
  taus <- x$taus

  switch(type,
    "surface" = plot_surface(coefs, se, taus, which,
                             ggtitle = ggtitle, xlab = xlab, ylab = ylab, ...),
    "slice_u" = plot_slices_u(coefs, se, taus, taus_slice, which, ncol, shared_ylim, ylim,
                              ylab = ylab, xlab = xlab, ggtitle = ggtitle, panel_prefix = panel_prefix,
                              line_color = line_color, ribbon_fill = ribbon_fill,
                              point_color = point_color, ribbon_alpha = ribbon_alpha, ...),
    "slice_v" = plot_slices_v(coefs, se, taus, taus_slice, which, ncol, shared_ylim, ylim,
                              ylab = ylab, xlab = xlab, ggtitle = ggtitle, panel_prefix = panel_prefix,
                              line_color = line_color, ribbon_fill = ribbon_fill,
                              point_color = point_color, ribbon_alpha = ribbon_alpha, ...)
  )
}


#' Internal: Surface Plot
#'
#' @keywords internal
plot_surface <- function(coefs, se, taus, coef_name, ggtitle = NULL, xlab = NULL, ylab = NULL, ...) {

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
    message("Package 'metR' is not installed. Using geom_tile for the surface plot. ",
            "Install metR for filled contour plots: install.packages(\"metR\")")
    p <- ggplot2::ggplot(tidy_mat, ggplot2::aes(x = .data$v, y = .data$u)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$effect)) +
      ggplot2::scale_fill_gradientn(name = "Effect", colours = myColors)
  }

  # Resolve label defaults
  if (is.null(ylab)) ylab <- expression(u)
  if (is.null(xlab)) xlab <- expression(v)
  if (is.null(ggtitle)) ggtitle <- coef_name

  p <- p +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 15)) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::ggtitle(ggtitle) +
    ggplot2::theme_minimal()

  return(p)
}


#' Compute Global Y-axis Limits for Slice Plots
#'
#' @keywords internal
compute_slice_ylim <- function(coefs, se, taus, taus_slice, slice_dim = c("u", "v")) {
  slice_dim <- match.arg(slice_dim)

  all_y_values <- numeric(0)

  for (k in seq_along(taus_slice)) {
    slice_val <- taus_slice[k]
    slice_idx <- which(abs(taus - slice_val) < 1e-6)

    if (length(slice_idx) == 0) next

    # Extract coefficient slice
    if (slice_dim == "u") {
      beta <- coefs[slice_idx, ]
      se_vals <- if (!is.null(se)) se[slice_idx, ] else NULL
    } else {
      beta <- coefs[, slice_idx]
      se_vals <- if (!is.null(se)) se[, slice_idx] else NULL
    }

    all_y_values <- c(all_y_values, beta)

    # Include confidence bounds if available
    if (!is.null(se_vals)) {
      crit <- stats::qnorm(0.975)
      all_y_values <- c(all_y_values, beta - se_vals * crit, beta + se_vals * crit)
    }
  }

  # Compute range with 5% margin
  y_range <- range(all_y_values, na.rm = TRUE)
  y_margin <- 0.05 * diff(y_range)

  c(y_range[1] - y_margin, y_range[2] + y_margin)
}


#' Internal: Slice Plot Fixing u
#'
#' @keywords internal
plot_slices_u <- function(coefs, se, taus, taus_slice, coef_name, ncol, shared_ylim, ylim,
                          ylab = NULL, xlab = NULL, ggtitle = NULL, panel_prefix = NULL,
                          line_color = "#11333a", ribbon_fill = "#0faeb6",
                          point_color = NULL, ribbon_alpha = 0.2, ...) {

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

  # Resolve NULL defaults
  if (is.null(panel_prefix)) panel_prefix <- "u = "
  if (is.null(point_color))  point_color  <- line_color
  if (is.null(ylab))         ylab         <- "Point Estimate"
  if (is.null(xlab))         xlab         <- expression(v)
  if (is.null(ggtitle))      ggtitle      <- coef_name

  # Compute global y-axis limits if shared_ylim is TRUE and ylim not provided
  if (shared_ylim && is.null(ylim)) {
    ylim <- compute_slice_ylim(coefs, se, taus, taus_slice, slice_dim = "u")
  }

  # Determine which panels show y-axis (only hide if shared_ylim is TRUE)
  if (shared_ylim) {
    show_y_idx <- seq(1, length(taus_slice), by = ncol)
  } else {
    show_y_idx <- seq_along(taus_slice)  # Show all y-axes
  }

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

    title <- paste0(panel_prefix, u_val)

    p <- ggplot2::ggplot(data = dd) +
      ggplot2::theme_minimal() +
      ggplot2::aes(x = .data$v, y = .data$beta) +
      ggplot2::geom_line(linewidth = 0.5, color = line_color) +
      ggplot2::geom_point(color = point_color, size = 1) +
      ggplot2::labs(title = title, x = xlab) +
      ggplot2::coord_cartesian(ylim = ylim)

    if (!all(is.na(se_vals))) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lb, ymax = .data$ub),
                                     alpha = ribbon_alpha, fill = ribbon_fill)
    }

    # Only show y-axis label on first column
    if (k %in% show_y_idx) {
      p <- p + ggplot2::ylab(ylab)
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
    patchwork::plot_annotation(title = ggtitle)

  return(combined)
}


#' Internal: Slice Plot Fixing v
#'
#' @keywords internal
plot_slices_v <- function(coefs, se, taus, taus_slice, coef_name, ncol, shared_ylim, ylim,
                          ylab = NULL, xlab = NULL, ggtitle = NULL, panel_prefix = NULL,
                          line_color = "#11333a", ribbon_fill = "#0faeb6",
                          point_color = NULL, ribbon_alpha = 0.2, ...) {

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

  # Resolve NULL defaults
  if (is.null(panel_prefix)) panel_prefix <- "v = "
  if (is.null(point_color))  point_color  <- line_color
  if (is.null(ylab))         ylab         <- "Point Estimate"
  if (is.null(xlab))         xlab         <- expression(u)
  if (is.null(ggtitle))      ggtitle      <- coef_name

  # Compute global y-axis limits if shared_ylim is TRUE and ylim not provided
  if (shared_ylim && is.null(ylim)) {
    ylim <- compute_slice_ylim(coefs, se, taus, taus_slice, slice_dim = "v")
  }

  # Determine which panels show y-axis (only hide if shared_ylim is TRUE)
  if (shared_ylim) {
    show_y_idx <- seq(1, length(taus_slice), by = ncol)
  } else {
    show_y_idx <- seq_along(taus_slice)  # Show all y-axes
  }

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

    title <- paste0(panel_prefix, v_val)

    p <- ggplot2::ggplot(data = dd) +
      ggplot2::theme_minimal() +
      ggplot2::aes(x = .data$u, y = .data$beta) +
      ggplot2::geom_line(linewidth = 0.5, color = line_color) +
      ggplot2::geom_point(color = point_color, size = 1) +
      ggplot2::labs(title = title, x = xlab) +
      ggplot2::coord_cartesian(ylim = ylim)

    if (!all(is.na(se_vals))) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lb, ymax = .data$ub),
                                     alpha = ribbon_alpha, fill = ribbon_fill)
    }

    if (k %in% show_y_idx) {
      p <- p + ggplot2::ylab(ylab)
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
    patchwork::plot_annotation(title = ggtitle)

  return(combined)
}
