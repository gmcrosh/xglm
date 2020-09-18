#' @useDynLib "xglm", .registration = TRUE
NULL

#' Fit a GLM
#'
#' @param x a design matrix
#' @param y the response
#' @param family the family
#' @param maxit maximum number of iteration
#' @param tol tolerance for convergence
#'
#' @export
glm_fit <- function(x, y, family, maxit = 1000, tol = 1e-10) {
  rcpp_glm_fit(x, y, family, maxit, tol)
}
