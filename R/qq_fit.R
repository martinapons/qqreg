#' Fit a Quantile on Quantiles Regression Model
#'
#' Estimates a two-dimensional quantile surface capturing heterogeneity
#' both within and between groups using the two-step quantile regression
#' estimator from Pons (2026).
#'
#' @param formula A formula specifying the model, e.g., \code{y ~ x1 + x2}.
#' @param data A data frame containing the variables in the formula.
#' @param group A character string naming the grouping variable in \code{data}.
#' @param taus Numeric vector of quantile levels, default \code{seq(0.1, 0.9, 0.1)}.
#' @param n_small Numeric, minimum group size (default NULL). Groups with fewer
#'   than \code{n_small} observations will be dropped with a message. If any group
#'   has fewer observations than the number of covariates, the function will stop
#'   with an error.
#' @param bootstrap Logical, whether to compute bootstrap standard errors (default TRUE).
#' @param nboot Number of bootstrap replications (default 1000).
#' @param parallel Logical, whether to parallelize computation (default TRUE).
#' @param parallel_first Logical, whether to parallelize the first stage across
#'   groups (default TRUE).
#' @param ncores Number of cores to use. If NULL (default), uses 80% of available cores.
#' @param weights Optional observation weights.
#'
#' @return An object of class \code{"qqfit"}, a list containing:
#'   \item{coef}{Array of coefficients (k x n_tau x n_tau)}
#'   \item{se}{Array of bootstrap standard errors (k x n_tau x n_tau), NULL if bootstrap = FALSE}
#'   \item{coef_boot}{Array of bootstrap coefficient draws (k x n_tau x n_tau x nboot), NULL if bootstrap = FALSE}
#'   \item{fitted}{Array of second-stage fitted values (n x n_tau x n_tau)}
#'   \item{fitted_first}{Matrix of first-stage fitted values (n x n_tau)}
#'   \item{rank}{Matrix of within-group ranks (n x n_tau)}
#'   \item{residuals}{Array of second-stage residuals (n x n_tau x n_tau)}
#'   \item{call}{The matched call}
#'   \item{formula}{The model formula}
#'   \item{taus}{The quantile grid}
#'   \item{nobs}{Number of observations}
#'   \item{ngroups}{Number of groups}
#'   \item{var_types}{List of variable classifications}
#'   \item{data}{The original data (sorted by group)}
#'   \item{collapsed}{Logical, whether data was collapsed to group level}
#'
#' @details
#' The model estimates:
#' \deqn{Q(v, Q(u, y_{ij} | x_{1,ij}, x_{2,j}) | x_{1,ij}, x_{2,j}) = x_1'\beta(u,v) + x_2'\gamma(u,v) + \alpha(u,v)}
#'
#' where:
#' \itemize{
#'   \item \eqn{u} is the within-group quantile index
#'   \item \eqn{v} is the between-group quantile index
#'   \item \eqn{x_1} are covariates that vary within groups
#'   \item \eqn{x_2} are covariates that vary only between groups
#' }
#'
#' Variables are automatically classified as:
#' \itemize{
#'   \item \strong{within}: varies within ALL groups (used in first stage)
#'   \item \strong{between}: constant within ALL groups (second stage only)
#'   \item \strong{mixed}: varies within SOME groups (included appropriately per group)
#' }
#'
#' If no within-group varying variables exist, the first stage computes
#' sample quantiles within each group, and the data is collapsed to the
#' group level with weights equal to group sizes.
#'
#' @examples
#' \donttest{
#' # Simulated data with heterogeneous effects
#' set.seed(123)
#' n_groups <- 50
#' n_per_group <- 50
#' n <- n_groups * n_per_group
#'
#' # Covariates
#' x1 <- 1 + rnorm(n) + rep(runif(n_groups), each = n_per_group)
#' x2 <- rep(rnorm(n_groups), each = n_per_group)
#'
#' # Shocks: eta (group-level), nu (individual-level)
#' eta <- rep(rnorm(n_groups), each = n_per_group)
#' nu <- rnorm(n)
#'
#' # DGP: y = 1 + x1 + x2 + eta*(1 - 0.1*x1 - 0.1*x2) + nu*(1 + 0.1*x1 + 0.1*x2)
#' y <- 1 + x1 + x2 + eta * (1 - 0.1*x1 - 0.1*x2) + nu * (1 + 0.1*x1 + 0.1*x2)
#'
#' data <- data.frame(
#'   group = rep(1:n_groups, each = n_per_group),
#'   x1 = x1,
#'   x2 = x2,
#'   y = y
#' )
#'
#' # Fit model
#' fit <- qq_fit(y ~ x1 + x2, data = data, group = "group",
#'               taus = seq(0.1, 0.9, 0.1), nboot = 200)
#'
#' # View results
#' summary(fit)
#' plot(fit, which = "x1", type = "surface")
#' plot(fit, which = "x1", type = "slice_u")
#' plot(fit, which = "x1", type = "slice_u", taus_slice = c(0.2, 0.8))
#'
#' # Intercept-only model (no covariates)
#' fit_intercept <- qq_fit(y ~ 1, data = data, group = "group",
#'                         taus = seq(0.1, 0.9, 0.1), nboot = 200)
#'
#' summary(fit_intercept)
#' plot(fit_intercept, which = "(Intercept)", type = "surface")
#' plot(fit_intercept, which = "(Intercept)", type = "slice_u")
#' }
#'
#' @export
qq_fit <- function(formula, data, group,
                   taus = seq(0.1, 0.9, by = 0.1),
                   n_small = NULL,
                   bootstrap = TRUE, nboot = 1000,
                   parallel = TRUE, parallel_first = TRUE,
                   ncores = NULL, weights = NULL) {

  # Capture call
  cl <- match.call()

  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }
  if (!is.character(group) || length(group) != 1) {
    stop("'group' must be a character string naming a column in data")
  }
  if (!group %in% names(data)) {
    stop("'group' variable '", group, "' not found in data")
  }
  if (!is.numeric(taus) || any(taus <= 0) || any(taus >= 1)) {
    stop("'taus' must be numeric values between 0 and 1 (exclusive)")
  }

  # Handle parallel/ncores logic
  if (parallel && is.null(ncores)) {
    ncores <- max(1, floor(0.8 * future::availableCores()))
  }
  if (parallel && !is.null(ncores) && ncores == 1) {
    parallel <- FALSE
  }

  # Sort data by group
  data <- data[order(data[[group]]), ]
  group_vec <- data[[group]]

  # Compute group sizes
  group_sizes <- table(group_vec)

  # Count covariates (need this for critical check)
  # Use a small subset to avoid memory issues
  X_temp <- stats::model.matrix(formula, data[seq_len(min(100, nrow(data))), , drop = FALSE])
  n_covariates <- ncol(X_temp)

  # CRITICAL CHECK: Groups smaller than number of covariates
  small_groups_idx <- which(group_sizes < n_covariates)
  if (length(small_groups_idx) > 0) {
    groups <- as.numeric(names(group_sizes))
    small_group_vals <- groups[small_groups_idx]

    message(sprintf("ERROR: %d group(s) have fewer observations than covariates (%d):",
                    length(small_groups_idx), n_covariates))
    message(paste("  Groups:", paste(small_group_vals, collapse = ", ")))
    message(sprintf("Either remove these groups or use n_small >= %d to drop them.", n_covariates))
    stop("Cannot estimate: groups with fewer observations than covariates")
  }

  # Drop small groups if n_small specified
  if (!is.null(n_small)) {
    if (!is.numeric(n_small) || n_small < 1) {
      stop("'n_small' must be a positive number")
    }

    drop_idx <- which(group_sizes < n_small)

    if (length(drop_idx) > 0) {
      groups <- as.numeric(names(group_sizes))
      drop_group_vals <- groups[drop_idx]

      # Keep observations NOT in dropped groups
      keep_mask <- !(group_vec %in% drop_group_vals)
      data <- data[keep_mask, ]
      group_vec <- group_vec[keep_mask]

      # Message
      message(sprintf("Dropped %d group(s) with fewer than %d observations:",
                      length(drop_idx), n_small))
      message(paste("  Groups:", paste(drop_group_vals, collapse = ", ")))
      message(sprintf("Remaining: %d groups, %d observations",
                      length(unique(group_vec)), nrow(data)))
    }
  }

  # Extract model components
  mf <- stats::model.frame(formula, data, na.action = na.pass)
  y <- stats::model.response(mf)

  # Detect variable types
  var_types <- detect_var_types(formula, data, group)

  # Build design matrices
  X_full <- stats::model.matrix(formula, data)

  # Identify columns for first stage (within + mixed)
  first_stage_vars <- c(var_types$within, var_types$mixed)
  if (var_types$has_intercept) {
    first_stage_vars <- c("(Intercept)", first_stage_vars)
  }

  # Check if we have any within-group variation
  has_within <- length(var_types$within) > 0 || length(var_types$mixed) > 0
  collapsed <- !has_within

  if (collapsed) {
    # No within-group variation: collapse to group level
    message("No within-group varying variables detected. ",
            "Collapsing to group level with weighted second stage.")

    fs_result <- qq_first_stage_collapsed(y, group_vec, taus)
    fitted_first <- fs_result$fitted
    group_weights <- fs_result$weights
    group_order <- fs_result$group_order

    # Build group-level design matrix (one row per group)
    # Take first observation from each group
    first_idx <- match(group_order, group_vec)
    X_second <- X_full[first_idx, , drop = FALSE]
    group_vec_collapsed <- group_vec[first_idx]

    # Run second stage with weights
    second <- qq_second_stage(fitted_first, X_second, taus, weights = group_weights)

    # Compute ranks (on collapsed data)
    rank <- compute_ranks(second$residuals, group_vec_collapsed, taus)

    # Store collapsed info
    nobs_original <- nrow(data)
    data_out <- data[first_idx, ]

  } else {
    # Normal case: within-group variation exists
    X_first <- X_full[, first_stage_vars, drop = FALSE]

    # Run first stage
    fitted_first <- qq_first_stage(y, X_first, group_vec, taus, var_types,
                                    parallel = parallel_first, ncores = ncores)

    # Run second stage
    second <- qq_second_stage(fitted_first, X_full, taus, weights = weights)

    # Compute ranks
    rank <- compute_ranks(second$residuals, group_vec, taus)

    nobs_original <- nrow(data)
    data_out <- data
    group_vec_collapsed <- group_vec
    X_second <- X_full
  }

  # Bootstrap (if requested)
  if (bootstrap) {
    message("Running bootstrap with ", nboot, " replications...")
    boot_result <- qq_bootstrap(
      fitted_first = fitted_first,
      X = X_second,
      group_vec = if (collapsed) group_vec_collapsed else group_vec,
      taus = taus,
      nboot = nboot,
      parallel = parallel,
      ncores = ncores,
      weights = if (collapsed) group_weights else weights
    )
    se <- boot_result$se
    coef_boot <- boot_result$coef_boot
  } else {
    se <- NULL
    coef_boot <- NULL
  }

  # Build result object
  result <- structure(
    list(
      coef = second$coef,
      se = se,
      coef_boot = coef_boot,
      fitted = second$fitted,
      fitted_first = fitted_first,
      rank = rank,
      residuals = second$residuals,
      call = cl,
      formula = formula,
      taus = taus,
      nobs = nobs_original,
      ngroups = length(unique(group_vec)),
      var_types = var_types,
      data = data_out,
      collapsed = collapsed
    ),
    class = "qqfit"
  )

  return(result)
}
