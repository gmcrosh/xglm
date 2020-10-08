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
glm_fit <- function(x, y, family, sample_weights = NULL, offset = NULL, maxit = 1000, tol = 1e-8) {
  nobs <- length(y)
  if (is.null(sample_weights)) {
    sample_weights <- rep(1, nobs)
  }
  if (is.null(offset)) {
    offset <- rep(0, nobs)
  }
  rcpp_glm_fit(x, y, sample_weights, offset, family, maxit, tol)
}


#' Fit a GLMnet
#'
#' @param x a design matrix
#' @param y the response
#' @param family the family
#' @param weights weights for each sample
#' @param offset offset for the model
#' @param maxit maximum number of iteration
#' @param tol tolerance for convergence
#' @param lambda lambda
#' @param alpha alpha
#'
#' @export
glmnet_fit <- function(x, y, family, sample_weights = NULL, offset = NULL, maxit = 1000, tol = 1e-8, lambda, alpha = 1) {
  nobs <- length(y)
  if (is.null(sample_weights)) {
    sample_weights <- rep(1, nobs)
  }
  if (is.null(offset)) {
    offset <- rep(0, nobs)
  }
  rcpp_glmnet_fit(x, y, sample_weights, offset, family, maxit, tol, lambda, alpha)
}
