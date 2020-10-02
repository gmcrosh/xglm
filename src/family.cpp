#include <RcppArmadillo.h>
#include "family.h"

arma::colvec ylogy(arma::colvec y, arma::colvec mu) {
  arma::colvec ans=arma::zeros(size(y));
  arma::uvec ids = arma::find(y != 0);
  ans.elem(ids) = (y.elem(ids) % arma::log(y.elem(ids) / mu.elem(ids)));
  return ans;
}

namespace Family {

arma::colvec Gaussian::variance(const arma::colvec& mu) const {
  return arma::ones<arma::colvec>(mu.n_elem);
}
arma::colvec Gaussian::initialize(const arma::colvec& y, const arma::colvec& weight) const {
  return y;
}
double Gaussian::deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const {
  return arma::accu(weight % arma::square(y - mu));
}


arma::colvec Binomial::variance(const arma::colvec& mu) const {
  return mu % (1.0 - mu);
}
arma::colvec Binomial::initialize(const arma::colvec& y, const arma::colvec& weight) const {
  return ((weight % y) + 0.5)/(1 + weight);
}
double Binomial::deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const {
  return arma::accu(2 * weight % (ylogy(y, mu) + ylogy(1 - y, 1 - mu)));
}

arma::colvec Poisson::variance(const arma::colvec& mu) const {
  return mu;
}
arma::colvec Poisson::initialize(const arma::colvec& y, const arma::colvec& weight) const {
  return (y + 0.1);
}
double Poisson::deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const {
  arma::colvec r = mu % weight;
  arma::uvec ids = arma::find(y > 0);
  r.elem(ids) = (weight.elem(ids) % (y.elem(ids) % arma::log(y.elem(ids) / mu.elem(ids)) - (y.elem(ids) - mu.elem(ids))));
  return arma::accu(2 * r);
}

arma::colvec Gamma::variance(const arma::colvec& mu) const {
  return arma::square(mu);
}
arma::colvec Gamma::initialize(const arma::colvec& y, const arma::colvec& weight) const {
  return y;
}
double Gamma::deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const {
  arma::uvec ids = arma::find(y == 0);
  arma::colvec logval = (y / mu);
  logval.elem(ids).ones();
  return arma::accu(-2 * weight * (arma::log(logval) - ((y - mu)/mu)));
}

arma::colvec Tweedie::variance(const arma::colvec& mu) const {
  return arma::pow(mu, varp);
}
arma::colvec Tweedie::initialize(const arma::colvec& y, const arma::colvec& weight) const {
  arma::colvec ycopy = y;
  return ycopy.replace(0, 0.1);
}
double Tweedie::deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const {
  arma::colvec y1 = y;
  y1.replace(0, 0.1);
  return arma::accu(
    2 * weight % ( -1 * ((arma::pow(y, 2 -varp) - arma::pow(mu, 2 -varp)) / (2 -varp)) +
      (y % (arma::pow(y1, 1 -varp) - arma::pow(mu, 1 -varp)) / (1 -varp))));
}


}
