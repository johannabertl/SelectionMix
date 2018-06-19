#' The function Q for the EM algorithm for the negative binomial mixture model
#'
#' The function negQ = -Q can be used for minimization.
#'
#' @param theta c(alpha, beta, p1, p2, ..., pk-1)
#' @param theta.prime c(alpha.prime, beta.prime, p1.prime, p2.prime, ..., p(k-1).prime)
#' @param x.syn vector of synonymous mutation count
#' @param x.non vector of non-synonymous mutation count
#' @param cvec vector of k positive values c(c1, ..., ck)
#'
#' @author Johanna Bertl
#'
#' @examples
#'
#' ### Example with 2 components ###
#'
#' # Simulate dataset of synonymous and non-synonymous mutations
#'
#' c1 = 0.5; c2 = 10
#' p1 = 0.2; p2 = 0.8
#' alpha = 10
#' beta = 5
#'
#' mutations = rnegbinmix_syn_nonsyn(n=3000, alpha=alpha, beta=beta, c=c(c1, c2), p = c(p1, p2))
#'
#' # Q
#'
#' Qtest = Q(theta = c(5, 1, 0.5), theta.prime = c(5, 1, 0.5), x.syn = mutations$Syn, x.non = mutations$Non, c = c(c1, c2))
#'
#'
#' ### Example with 3 components ###
#'
#' c1 = 0.01; c2 = 1; c3 = 100
#' p1 = 0.33; p2 = 0.33; p3 = 0.34
#' alpha = 10
#' beta = 5
#'
#' mutations = rnegbinmix_syn_nonsyn(n=3000, alpha=alpha, beta=beta, c=c(c1, c2, c3), p = c(p1, p2, p3))
#'
#' # Q
#'
#' Qtest = Q(theta = c(5, 1, 0.3, 0.3), theta.prime = c(5, 1, 0.3, 0.3), x.syn = mutations$Syn, x.non = mutations$Non, c = c(c1, c2, c3))
#'
#' @export

Q = function(theta, theta.prime, x.syn, x.non, cvec){

  ### Preparing parameters ###

  # Model parameters:

  alpha = theta[1]
  beta = theta[2]
  p = theta[3:length(theta)]
  p = c(p, 1-sum(p))

  # Model parameters of the previous step (theta prime)

  alpha.prime = theta.prime[1]
  beta.prime = theta.prime[2]
  p.prime = theta.prime[3:length(theta.prime)]
  p.prime = c(p.prime, 1-sum(p.prime))

  # further parameters

  n = length(x.syn) # number of genes
  k = length(p) # number of mixture components


  #### Computation of the expected log-likelihood ####

  # matrix qmat consisting of columns q1, q2, ..., qk

  qmat.temp = matrix(NA, ncol=k, nrow=n)

  for(i in 1:k){
    qmat.temp[,i] = ( (1 - (1/(beta.prime/cvec[i] + 1)))^alpha.prime ) *
      ( (beta.prime/cvec[i] + 1)^(-x.non) ) * p.prime[i]
  }

  qmat = qmat.temp/rowSums(qmat.temp)


  # expected log likelihood of x.syn

  ll.x.syn = sum(log(dnegbin(x.syn, alpha, 1/(beta + 1))))


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

  Q = ll.x.syn + ll.x.non

  return(Q)
}


negQ = function(theta, theta.prime, x.syn, x.non, cvec){
  -Q(theta, theta.prime, x.syn, x.non, cvec)
}
