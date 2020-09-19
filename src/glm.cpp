#include "glm.h"

arma::mat glm_fit(const arma::mat& x, 
                  const arma::colvec& y,
                  const arma::colvec& weights,
                  const arma::colvec& offset,
                  const Family::ExponentialFamily& family,
                  int maxit, 
                  double tol) {
  const int n_cols = x.n_cols;
  arma::mat Q, R;
  arma::colvec s = arma::zeros<arma::colvec>(n_cols);
  arma::colvec s_old;
  arma::colvec eta = family.link_fun(family.initialize(y));
  arma::qr_econ(Q, R, x);
  for (int i = 0; i < maxit; i++) {
    s_old = s;
    const arma::colvec mu = family.link_inverse(eta);
    const arma::colvec mu_p = family.link_mu_eta(eta);
    const arma::colvec z = (eta - offset) + (y - mu) / mu_p;
    const arma::colvec W = (weights % arma::square(mu_p)) / family.variance(mu);
    const arma::mat C = arma::chol(Q.t() * (Q.each_col() % W));
    const arma::colvec s1 = arma::solve(arma::trimatl(C.t()), Q.t() * (W % z));
    s = arma::solve(arma::trimatu(C), s1);
    eta = Q * s + offset;

    const bool is_converged = std::sqrt(arma::accu(arma::square(s - s_old))) < tol;
    if (is_converged) break;
  }
  return arma::solve(arma::trimatu(R), Q.t() * (eta - offset));
}
