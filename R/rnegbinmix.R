#' Simulating data from a negative binomial mixture distribution
#'
#' rnegbinmix simulates data from a negative binomial mixture distribution with k components.
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
#'
#' @examples
#'
#' # Test the negative binomial density function and how well it fits to simulated values
#'
#' k = 0:20
#' density = dnegbin.alphabeta(k=0:20, alpha=10, beta=5)
#'
#' simulated = rnegbin.alphabeta(1000, alpha=10, beta=5)
#' plot(table(simulated)/1000)
#'
#' lines(k, density, col="purple")
#'
#'
#' # Same for a mixture of two negative binomial densities
#'
#' k = 0:50
#' c1 = 0.5; c2 = 10
#' p1 = 0.2; p2 = 0.8
#' alpha = 10
#' beta = 5
#'
#' density1 = dnegbin.alphabeta(k = k, alpha = alpha, beta = beta/c1)
#' plot(k, density1, t="b")
#' density2 = dnegbin.alphabeta(k = k, alpha = alpha, beta = beta/c2)
#' plot(k, density2, t="b")
#' density.mixture = p1 * density1 + p2 * density2
#' plot(k, density.mixture, t="b")
#' 
#' density.mixture2 = dnegbinmix(k = k, alpha = alpha, beta = beta, cvec = c(c1, c2), p = c(p1, p2))
#' lines(k, density.mixture2, col="red")
#'
#' simulated.mixture = rnegbinmix(1000, alpha = alpha, beta = beta, c = c(c1, c2), p = c(p1, p2))
#' plot(table(simulated.mixture)/1000)
#' lines(k, density.mixture, col="purple")
#'
#'
#' # Mixture with 3 components
#'
#' k = 0:500
#' c1 = 0.01; c2 = 1; c3=100
#' p1 = 0.33; p2 = 0.33; p3 = 0.34
#' alpha = 10
#' beta = 5
#'
#' density1 = dnegbin.alphabeta(k = k, alpha = alpha, beta = beta/c1)
#' plot(k, density1, t="l")
#' density2 = dnegbin.alphabeta(k = k, alpha = alpha, beta = beta/c2)
#' plot(k, density2, t="l")
#' density3 = dnegbin.alphabeta(k = k, alpha = alpha, beta = beta/c3)
#' plot(k, density3, t="l")
#' density.mixture = p1 * density1 + p2 * density2 + p3 * density3
#' plot(k, density.mixture, t="l")
#'
#' simulated.mixture = rnegbinmix(1000, alpha = alpha, beta = beta, c = c(c1, c2, c3), p = c(p1, p2, p3))
#' plot(table(simulated.mixture)/1000)
#' lines(k, density.mixture, col="purple")

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
  rnegbin.alphabeta(n, alpha, beta/pi)
}


dnegbinmix = function(k, alpha, beta, cvec, p){
  if(!all(p>0)){
    stop("All entries of p must be positive.")
  }
  if(length(p)!=length(cvec)){
    stop("p and c must have the same length.")
  }
  if(sum(p) !=1 ) {
    p = p/sum(p)
    warning("p was scaled such that sum(p)=1.")
  }
  
  densmat = matrix(NA, ncol=length(cvec), nrow=length(k))
  for(i in 1:length(cvec)){
    densmat[,i] = dnegbin.alphabeta(k, alpha, beta/cvec[i])
  }
  pmat = matrix(p, ncol=length(cvec), nrow=length(k), byrow=T)
  rowSums(densmat*pmat)
  
}


dnegbin = function(k, r, p){
  choose(k + r - 1, k)*((1 - p)**r)*(p**k)
}

dnegbin.alphabeta = function(k, alpha, beta){
  dnegbin(k, r = alpha, p = 1/(beta + 1))
}

rnegbin.alphabeta = function(n, alpha, beta){
  # simulate lamda from a Gamma distribution
  lambda = rgamma(n, shape = alpha, rate = beta)
  # simulate negative binomial from a Poisson(lambda) distribution
  rpois(n, lambda)
}
