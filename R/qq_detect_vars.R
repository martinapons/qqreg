#' Detect Variable Types for QQ Regression
#'
#' Classifies variables in the model formula as within-group, between-group,
#' or mixed based on their variance within groups.
#'
#' @param formula A formula specifying the model.
#' @param data A data frame containing the variables.
#' @param group A character string naming the grouping variable.
#'
#' @return A list with three components:
#'   \item{within}{Character vector of variable names that vary within ALL groups}
#'   \item{between}{Character vector of variable names constant within ALL groups}
#'   \item{mixed}{Character vector of variable names that vary only within SOME groups}
#'
#' @details
#' Variables are classified based on their within-group variance:
#' \itemize{
#'   \item \strong{within}: variance > 0 in ALL groups
#'   \item \strong{between}: variance = 0 in ALL groups
#'   \item \strong{mixed}: variance > 0 in SOME groups but not all
#' }
#'
#' @keywords internal
detect_var_types <- function(formula, data, group) {
  # Extract model frame and design matrix
  mf <- stats::model.frame(formula, data, na.action = na.pass)
  X <- stats::model.matrix(formula, data)

  # Get group vector
  group_vec <- data[[group]]

  # Get variable names (excluding intercept)
  var_names <- colnames(X)
  has_intercept <- "(Intercept)" %in% var_names
  if (has_intercept) {
    var_names <- setdiff(var_names, "(Intercept)")
  }

  # If no covariates besides intercept, return empty lists

  if (length(var_names) == 0) {
    return(list(
      within = character(0),
      between = character(0),
      mixed = character(0),
      has_intercept = has_intercept
    ))
  }

  # Classify each variable
  var_types <- sapply(var_names, function(v) {
    # Compute within-group variance for each group
    within_var <- tapply(X[, v], group_vec, stats::var, na.rm = TRUE)

    # Handle single-observation groups (var = NA)
    within_var[is.na(within_var)] <- 0

    # Classify
    if (all(within_var > 1e-10)) {
      return("within")
    } else if (all(within_var < 1e-10)) {
      return("between")
    } else {
      return("mixed")
    }
  }, USE.NAMES = TRUE)

  # Return classified lists
  list(
    within = names(var_types)[var_types == "within"],
    between = names(var_types)[var_types == "between"],
    mixed = names(var_types)[var_types == "mixed"],
    has_intercept = has_intercept
  )
}


#' Get Group-Specific Variables
#'
#' For mixed variables, determines which variables have non-zero variance
#' within each specific group.
#'
#' @param X Design matrix.
#' @param group_vec Group indicator vector.
#' @param mixed_vars Character vector of mixed variable names.
#'
#' @return A named list where each element corresponds to a group and contains
#'   the names of mixed variables that vary within that group.
#'
#' @keywords internal
get_group_specific_vars <- function(X, group_vec, mixed_vars) {
  if (length(mixed_vars) == 0) {
    return(list())
  }

  groups <- unique(group_vec)

  group_vars <- lapply(groups, function(g) {
    idx <- which(group_vec == g)
    X_g <- X[idx, mixed_vars, drop = FALSE]

    # Find variables with non-zero variance in this group
    var_g <- apply(X_g, 2, stats::var, na.rm = TRUE)
    var_g[is.na(var_g)] <- 0
    names(var_g)[var_g > 1e-10]
  })

  names(group_vars) <- as.character(groups)
  group_vars
}
