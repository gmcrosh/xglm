#ifndef glmnet_H
#define glmnet_H

#include "family.h"
#include "linkfunctions.h"

arma::colvec glmnet_fit(const arma::mat& x,
                     const arma::colvec& y,
                     const arma::colvec& sample_weights,
                     const arma::colvec& offset,
                     const Family::ExponentialFamily& family,
                     int maxit,
                     double tol,
                     double lambda,
                     double alpha);

#endif
