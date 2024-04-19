#ifndef TUBER_MATOPS_H
#define TUBER_MATOPS_H

#include <RcppArmadillo.h>

inline arma::mat sqdist(const arma::mat& X, const arma::mat& Z) {
    arma::uword n = X.n_rows;
    arma::uword m = Z.n_rows;
    arma::mat res(n, m);
    arma::rowvec diff(X.n_cols);
    for ( arma::uword j = 0; j < m; ++j ) {
        for ( arma::uword i = 0; i < n; ++i ) {
            diff = X.row(i) - Z.row(j);
            res(i, j) = arma::dot(diff, diff);
        }
    }
    return res;
}

#endif
