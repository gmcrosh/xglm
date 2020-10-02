#include <RcppArmadillo.h>
#include "family.h"
#include "linkfunctions.h"
#include "glm.h"

// [[Rcpp::depends(RcppArmadillo)]]
std::unique_ptr<Link::LinkFunction> link_from_string(const std::string& link_name) {
  if (link_name == "logit") {
    std::unique_ptr<Link::LinkFunction> ptr(new Link::Logit());
    return ptr;
  } else if (link_name == "identity") {
    std::unique_ptr<Link::LinkFunction> ptr(new Link::Identity());
    return ptr;
  } else if (link_name == "log") {
    std::unique_ptr<Link::LinkFunction> ptr(new Link::Log());
    return ptr;
  } else if (link_name == "inverse") {
    std::unique_ptr<Link::LinkFunction> ptr(new Link::Inverse());
    return ptr;
  } else if (link_name == "probit") {
    std::unique_ptr<Link::LinkFunction> ptr(new Link::Probit());
    return ptr;
  } else if (link_name == "sqrt") {
    std::unique_ptr<Link::LinkFunction> ptr(new Link::Sqrt());
    return ptr;
  }
  Rcpp::stop("Link function not available.");
}

std::unique_ptr<Link::LinkFunction> link_tweedie(const double linkp) {
  if (linkp == 0) {
    std::unique_ptr<Link::LinkFunction> ptr(new Link::Log());
    return ptr;
  } else if (linkp == 1) {
    std::unique_ptr<Link::LinkFunction> ptr(new Link::Identity());
    return ptr;
  } else if (linkp == 2) {
    std::unique_ptr<Link::LinkFunction> ptr(new Link::Inverse());
    return ptr;
  }
  Rcpp::stop("Link function not available.");
}


template<class FamilyClass>
inline Rcpp::XPtr<FamilyClass> make_family(std::string link) {
  auto l = link_from_string(link);
  auto* family = new FamilyClass(l);
  Rcpp::XPtr<FamilyClass> pointer(family, true);
  return pointer;
}

// [[Rcpp::export]]
Rcpp::XPtr<Family::Gaussian> rcpp_make_gaussian(std::string link) {
  return make_family<Family::Gaussian>(link);
}

// [[Rcpp::export]]
Rcpp::XPtr<Family::Binomial> rcpp_make_binomial(std::string link) {
  return make_family<Family::Binomial>(link);
}

// [[Rcpp::export]]
Rcpp::XPtr<Family::Poisson> rcpp_make_poisson(std::string link) {
  return make_family<Family::Poisson>(link);
}

// [[Rcpp::export]]
Rcpp::XPtr<Family::Gamma> rcpp_make_gamma(std::string link) {
  return make_family<Family::Gamma>(link);
}

// [[Rcpp::export]]
Rcpp::XPtr<Family::Tweedie> rcpp_make_tweedie(double varp, double linkp) {
  auto l = link_tweedie(linkp);
  auto* family = new Family::Tweedie(varp, l);
  Rcpp::XPtr<Family::Tweedie> pointer(family, true);
  return pointer;
}

// [[Rcpp::export]]
arma::mat rcpp_glm_fit(const arma::mat& x, const arma::colvec& y,
                  const arma::colvec& sample_weights, const arma::colvec& offset,
                  Rcpp::XPtr<Family::ExponentialFamily> family,
                  int maxit, double tol) {
  return glm_fit(x, y, sample_weights, offset, *family, maxit, tol);
}
