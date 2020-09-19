#' @useDynLib "xglm", .registration = TRUE
#' @importFrom "Rcpp" "sourceCpp"
NULL

#' Fit a GLM
#'
#' @param x a design matrix
#' @param y the response
#' @param family the family
#' @param weights weights for each sample
#' @param offset offset for the model
#' @param maxit maximum number of iteration
#' @param tol tolerance for convergence
#'
#' @export
glm_fit <- function(x, y, family, weights = NULL, offset = NULL, maxit = 1000, tol = 1e-8) {
  nobs <- length(y)
  if (is.null(weights)) {
    weights <- rep(1, nobs)
  }
  if (is.null(offset)) {
    offset <- rep(0, nobs)
  }
  rcpp_glm_fit(x, y, weights, offset, family, maxit, tol)
}
