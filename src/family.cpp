#include <RcppArmadillo.h>
#include "family.h"

namespace Family {

arma::colvec Gaussian::variance(const arma::colvec& mu) const {
  return arma::ones<arma::colvec>(mu.n_elem);
}
arma::colvec Gaussian::initialize(const arma::colvec& y) const {
  return y;
}

arma::colvec Binomial::variance(const arma::colvec& mu) const {
  return mu % (1.0 - mu);
}
arma::colvec Binomial::initialize(const arma::colvec& y) const {
  return (y + 0.5)/(1);
}

arma::colvec Poisson::variance(const arma::colvec& mu) const {
  return mu;
}
arma::colvec Poisson::initialize(const arma::colvec& y) const {
  return (y + 0.1);
}

arma::colvec Gamma::variance(const arma::colvec& mu) const {
  return arma::square(mu);
}
arma::colvec Gamma::initialize(const arma::colvec& y) const {
  return y;
}

}
