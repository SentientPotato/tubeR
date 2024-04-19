#include <RcppArmadillo.h>
#include "matops.h"

// [[Rcpp::export(.sqdist)]]
arma::mat sqdist_export(const arma::mat& X, const arma::mat& Z) {
    return sqdist(X, Z);
}
