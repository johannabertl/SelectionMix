#' Simulating data from a negative binomial mixture distribution
#'
#' rnegbinmix simulates data from a negative binomial mixture distribution.
#'
#' The negative binomial mixture is defined as follows:
#'
#' The function \code{rnegbin} from package MASS is used for the simulation of negative binomial data.
#'
#'
#' @param n
#' @param alpha, scalar, positive
#' @param beta positive
#' @param c numeric, positive, length k
#' @param p numeric, probability vector, length k
#'
#'
#' @author Johanna Bertl

rnegbinmix = function(n, alpha, beta, c, p){
  if(!all(p>0)){
    stop("All entries of p must be positive.")
  }
  if(length(p)!=length(c)){
    stop("p and c must have the same length.")
  }
  if(sum(p) !=1 ) {
    p = p/sum(p)
    warning("p was scaled such that sum(p)=1.")
  }

  # simulate pi (for the mixing)
  pi = sample(x = c, size=n, replace=T, prob = p)
  # simulate from the negative binomial mixture distribution
  rnegbin(n, mu=1/(beta/pi + 1), theta = alpha)
}


dnegbin = function(k, r, p){
  choose(k + r - 1, k)*((1 - p)**r)*(p**k)
}
