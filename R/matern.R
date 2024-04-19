#' Matern Covariance Function
#'
#' The Matern kernel is characterized by a length scale \eqn{\ell}{l} and a
#' positive parameter \eqn{\nu}{nu} governing the function's smoothness.
#' If \eqn{\nu = d + 1/2}{nu = d + 1/2}, the corresponding GP is \eqn{d} times
#' differentiable and has a particularly nice form; currently only \eqn{d = 2}
#' is supported. The covariance between \eqn{\mathbf{x}_i}{x[i]} and
#' \eqn{\mathbf{x}_j}{x[j]} when \eqn{\nu = 5/2}{nu = 5/2} is given as
#' \deqn{
#' \left(1 + \sqrt{5}D_{ij} + (5/3)S_{ij}\right)\exp\left(-\sqrt{5}D_{ij}\right)
#' }
#' where \eqn{D_{ij}} is the distance between \eqn{\mathbf{x}_i}{x[i]} and
#' \eqn{\mathbf{x}_j}{x[j]} and \eqn{S_{ij}} is the squared distance.
#' Currently at this time \eqn{\ell}{l} is also constrained to be 1.
#'
#' @param X A numeric matrix; columns are variables and rows are observations
#' @param Z Another such matrix, optional
#'
#' @return A numeric matrix where the i, j entry gives the covariance between
#'     observation i of matrix X and observation j of matrix Z if Z is provided
#'     or observations i and j of matrix X if not.
#'
#' @export
cov_matern = function(X, Z) {
    if ( missing(Z) ) return(.matern_2_same(X))
    return(.matern_2_cross(X, Z))
}
