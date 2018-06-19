#' The function Q for the EM algorithm to estimate the parameters p1, ..., pk-1 in the negative binomial mixture model
#'
#' This function is used by the function EM_proportions. The parameters alpha and beta are assumed to be known. The function negQ_proportions = -Q_proportions can be used for minimization.
#'
#' @param theta c(p1, p2, ..., pk-1)
#' @param theta.prime c(p1.prime, p2.prime, ..., p(k-1).prime)
#' @param x.syn vector of synonymous mutation count
#' @param x.non vector of non-synonymous mutation count
#' @param alpha
#' @param beta
#' @param cvec vector of k positive values c(c1, ..., ck)
#'
#' @author Johanna Bertl
#'
#' @export

Q_proportions = function(theta, theta.prime, x.non, alpha, beta, cvec){

  ### Preparing parameters ###

  # Model parameters:
  p = theta
  p = c(p, 1-sum(p))

  # Model parameters of the previous step (theta prime)
  p.prime = theta.prime
  p.prime = c(p.prime, 1-sum(p.prime))

  # further parameters

  n = length(x.non) # number of genes
  k = length(p) # number of mixture components


  #### Computation of the expected log-likelihood ####

  # matrix qmat consisting of columns q1, q2, ..., qk

  qmat.temp = matrix(NA, ncol=k, nrow=n)

  for(i in 1:k){
    qmat.temp[,i] = ( (1 - (1/(beta/cvec[i] + 1)))^alpha ) *
      ( (beta/cvec[i] + 1)^(-x.non) ) * p.prime[i]
  }

  qmat = qmat.temp/rowSums(qmat.temp)


  # expected log likelihood of x.non

  summand = numeric(k)
  for(i in 1:k){
    vec = dnegbin(x.non, alpha, 1/(beta/cvec[i] + 1))*p[i]
    # if the entries of vec are too small, the logarithm can't be computed (in the next step), so they are replaced by the smallest possible value.
    if(!all(vec>=10^(-323))){
      num = sum(vec<10^(-323))
      warning(paste0("In component ", i, ", the log-likelihood of the non-synonymous mutations had to be approximated, because the likelihood contains ", num, " entries that are too small to compute the log."))
      vec = ifelse(vec<10^(-323), 10^(-323), vec)
    }
    summand[i] = sum(log(vec)*qmat[,i])
  }

  ll.x.non = sum(summand)


  # Q

  Q = ll.x.non

  return(Q)
}


negQ_proportions = function(theta, theta.prime, x.non, alpha, beta, cvec){
  -Q_proportions(theta, theta.prime, x.non, alpha, beta, cvec)
}
