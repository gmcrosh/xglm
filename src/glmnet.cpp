#include "glmnet.h"

arma::colvec casl_util_soft_thresh_rcpp(arma::colvec a, double b) {
  arma::uvec ids = arma::find(abs(a) <= b);
  a.elem(ids).zeros();
  ids = arma::find(a > 0);
  a.elem(ids) = a.elem(ids) - b;
  ids = arma::find(a < 0);
  a.elem(ids) = a.elem(ids) + b;
  return(a);
}

double casl_util_soft_thresh_rcpp(double a, double b) {
  if(std::abs(a) <= b) {
    return(0);
  } else if(a < 0) {
    return(a + b);
  } else {
    return(a - b);
  }
}

arma::colvec casl_lenet_update_beta_rcpp(
    const arma::mat& x,
    const arma::colvec& y,
    double lambda,
    double alpha,
    arma::colvec& b,
    const arma::colvec& w) {
  arma::mat wx = arma::diagmat(w) * x;
  arma::mat wx2 = arma::diagmat(w) * arma::square(x);
  arma::colvec xb = x * b;
  for(int i = 0; i < b.n_elem; i++) {
    xb = xb - x.col(i) * b[i];
    b[i] = arma::accu(casl_util_soft_thresh_rcpp(
      arma::accu(wx.col(i) % (y - xb)), alpha * lambda));
    b[i] = b[i] / (arma::accu(wx2.col(i)) + lambda * (1 - alpha));
    xb = xb + (x.col(i) * b[i]);
  }
  return(b);
}

arma::colvec glmnet_fit(const arma::mat& x,
                     const arma::colvec& y,
                     const arma::colvec& sample_weights,
                     const arma::colvec& offset,
                     const Family::ExponentialFamily& family,
                     int maxit,
                     double tol,
                     double lambda,
                     double alpha) {
  const int ncols = x.n_cols;
  const int nrow = x.n_rows;
  arma::mat Q, R;
  arma::colvec s = arma::zeros<arma::colvec>(ncols);
  arma::colvec eta = x * s + offset;
  double dev, devold;
  for (int i = 0; i < maxit; i++) {
    const arma::colvec mu = family.link_inverse(eta);
    const arma::colvec mu_p = family.link_mu_eta(eta);
    const arma::colvec z = (eta - offset) + (y - mu) / mu_p;
    const arma::colvec W = ((sample_weights % arma::square(mu_p)) / family.variance(mu)) / nrow;
    s = casl_lenet_update_beta_rcpp(x, z, lambda, alpha, s, W);
    eta = x * s + offset;
    devold = dev;
    dev = family.deviance(y, family.link_inverse(eta), sample_weights);
    if (i > 0) {
      const bool is_converged = (std::abs(dev - devold) / (0.1 + std::abs(dev)) < tol);
      if (is_converged) break;
    }
  }
  return s;
}
