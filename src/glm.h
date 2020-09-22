#ifndef glm_H
#define glm_H

#include "family.h"
#include "linkfunctions.h"

arma::mat glm_fit(const arma::mat& x,
                  const arma::colvec& y,
                  const arma::colvec& sample_weights,
                  const arma::colvec& offset,
                  const Family::ExponentialFamily& family,
                  int maxit,
                  double tol);

#endif
