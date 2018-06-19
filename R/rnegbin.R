#' The negative binomial distribution
#'
#' rnegbin.alphabeta simulates data from a negative binomial distribution and dnegbin and dnegbin.alphabeta compute the density.
#'
#' The functions with \code{alphabeta} use the parameterization with alpha and beta, while \code{dnegbin} uses the parameterization with r and p.
#'
#'
#' The functions \code{rgamma} and \code{rstats} from package stats are used for the simulation of negative binomial data.
#'
#' dnegbin and dnegbin.alphabeta can handle non-integer counts (without warning), because the binomial coefficient is implemented with a gamma function.
#'
#' @param n scalar, positive. Number of samples.
#' @param k numeric. Vector of observed counts.
#' @param alpha scalar, positive.
#' @param beta scalar, positive.
#' @param p scalar, probability.
#' @param r scalar, positive.
#'
#' @author Johanna Bertl
#'
#' @examples
#'
#' # see rnegbinmix
#'
#' @export
#'
dnegbin = function(k, r, p){
  if( p < 0 | p > 1 | r <= 0) {
    stop("Make sure that 0 < p < 1 and r > 0.")
  } else {
    (gamma(k+r)/(gamma(k+1)*gamma(r)))*((1 - p)**r)*(p**k)
  }
}

#' @rdname dnegbin
#' @export
dnegbin.alphabeta = function(k, alpha, beta){
  dnegbin(k, r = alpha, p = 1/(beta + 1))
}

#' @rdname dnegbin
#' @export
rnegbin.alphabeta = function(n, alpha, beta){
  # simulate lamda from a Gamma distribution
  lambda = rgamma(n, shape = alpha, rate = beta)
  # simulate negative binomial from a Poisson(lambda) distribution
  rpois(n, lambda)
}
