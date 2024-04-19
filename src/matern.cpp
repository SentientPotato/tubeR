#include <RcppArmadillo.h>
#include "matops.h"

#define SQRT_5 2.236067977499789805051477742381393909454345703125

// Cov(X, X)
// [[Rcpp::export(.matern_2_same)]]
arma::mat cov_matern_2(const arma::mat& X) {
    arma::uword n = X.n_rows;
    arma::mat   K = arma::eye(n, n);
    arma::mat   S = sqdist(X, X);
    arma::mat   D = arma::sqrt(S);
    for ( arma::uword j = 0; j < (n-1); ++j ) {
        for ( arma::uword i = j+1; i < n; ++i ) {
            double elem = (1.0 + SQRT_5 * D(i, j) + (5.0/3.0) * S(i, j));
            elem *= std::exp(-SQRT_5 * D(i, j));
            K(i, j) = elem;
            K(j, i) = elem;
        }
    }
    return K;
}

// Cov(X, Z)
// (can't export this one because overloading doesn't work with Rcpp's C API)
arma::mat cov_matern_2(const arma::mat& X, const arma::mat& Z) {
    arma::uword n = X.n_rows;
    arma::uword m = Z.n_rows;
    arma::mat   K(n, m);
    arma::mat   S = sqdist(X, Z);
    arma::mat   D = arma::sqrt(S);
    for ( arma::uword j = 0; j < n; ++j ) {
        for ( arma::uword i = 0; i < n; ++i ) {
            double elem = (1.0 + SQRT_5 * D(i, j) + (5.0/3.0) * S(i, j));
            elem *= std::exp(-SQRT_5 * D(i, j));
            K(i, j) = elem;
        }
    }
    return K;
}

// So we export this one with a slightly different name instead
// [[Rcpp::export(.matern_2_cross)]]
arma::mat cov_matern_2_cross(const arma::mat& X, const arma::mat& Z) {
    arma::uword n = X.n_rows;
    arma::uword m = Z.n_rows;
    arma::mat   K(n, m);
    arma::mat   S = sqdist(X, Z);
    arma::mat   D = arma::sqrt(S);
    for ( arma::uword j = 0; j < n; ++j ) {
        for ( arma::uword i = 0; i < n; ++i ) {
            double elem = (1.0 + SQRT_5 * D(i, j) + (5.0/3.0) * S(i, j));
            elem *= std::exp(-SQRT_5 * D(i, j));
            K(i, j) = elem;
        }
    }
    return K;
}
