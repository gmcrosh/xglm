#ifndef family_H
#define family_H

#include <RcppArmadillo.h>
#include "linkfunctions.h"

namespace Family {

class ExponentialFamily {
  // TODO: make variance a parameter as well
public:
  arma::colvec link_fun(const arma::colvec& eta) const {
    return link_function->link(eta);
  }
  arma::colvec link_inverse(const arma::colvec& eta) const {
    return link_function->inverse(eta);
  }
  arma::colvec link_mu_eta(const arma::colvec& eta) const {
    return link_function->mu_eta(eta);
  }
  virtual arma::colvec variance(const arma::colvec& mu) const = 0;
  virtual arma::colvec initialize(const arma::colvec& y, const arma::colvec& weight) const = 0;
  virtual double deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const = 0;
  ExponentialFamily(std::unique_ptr<Link::LinkFunction>& link) : link_function(std::move(link)) {
  }
  virtual ~ExponentialFamily() {}
private:
  std::unique_ptr<Link::LinkFunction> link_function;
};

class Gaussian : public ExponentialFamily {
public:
  arma::colvec variance(const arma::colvec& mu) const;
  arma::colvec initialize(const arma::colvec& y, const arma::colvec& weight) const;
  double deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const;
  Gaussian(std::unique_ptr<Link::LinkFunction>& link) : ExponentialFamily(link) {}
};

class Binomial : public ExponentialFamily {
public:
  arma::colvec variance(const arma::colvec& mu) const;
  arma::colvec initialize(const arma::colvec& y, const arma::colvec& weight) const;
  double deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const;
  Binomial(std::unique_ptr<Link::LinkFunction>& link) : ExponentialFamily(link) {}
};

class Poisson : public ExponentialFamily {
public:
  arma::colvec variance(const arma::colvec& mu) const;
  arma::colvec initialize(const arma::colvec& y, const arma::colvec& weight) const;
  double deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const;
  Poisson(std::unique_ptr<Link::LinkFunction>& link) : ExponentialFamily(link) {}
};

class Gamma : public ExponentialFamily {
public:
  arma::colvec variance(const arma::colvec& mu) const;
  arma::colvec initialize(const arma::colvec& y, const arma::colvec& weight) const;
  double deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const;
  Gamma(std::unique_ptr<Link::LinkFunction>& link) : ExponentialFamily(link) {}
};

class Tweedie : public ExponentialFamily {
public:
  double varp;
  arma::colvec variance(const arma::colvec& mu) const;
  arma::colvec initialize(const arma::colvec& y, const arma::colvec& weight) const;
  double deviance(const arma::colvec& y, const arma::colvec& mu, const arma::colvec& weight) const;
  Tweedie(double v, std::unique_ptr<Link::LinkFunction>& link) : ExponentialFamily(link) {
    varp = v;
  }
};

}
#endif
