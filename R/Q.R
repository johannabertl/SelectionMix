#' The function Q for the negative binomial mixture model
#'
#' Currently only for mixtures of 2 negbinoms.
#'
#' The function negQ is -Q (for minimization).
#'
#' The function Q_fixed_prop is for the function EM_fixed_prop, where the mixture proportions are assumed fixed. Here, theta = c(alpha, beta).
#'
#' @param theta c(alpha, beta, p1)
#' @param theta.prime c(alpha.prime, beta.prime, p1.prime)
#' @param x.syn vector of synonymous mutation count
#' @param x.non vector of non-synonymous mutation count
#' @param c1 positive scalar
#' @param c2 positive scalar
#'
#' @author Johanna Bertl

Q = function(theta, theta.prime, x.syn, x.non, c1, c2){

  # Parameters:

  alpha = theta[1]
  beta = theta[2]
  p1 = theta[3]
  p2 = 1-p1

  alpha.prime = theta.prime[1]
  beta.prime = theta.prime[2]
  p1.prime = theta.prime[3]
  p2.prime = 1-p1.prime

  n = length(x.syn)


  # q1 and q2 (vectors)

  q1.el1 = ( (1 - (1/(beta.prime/c1 + 1)))^alpha.prime ) *
    ( (beta.prime/c1 + 1)^(-x.non) ) * p1.prime
  q1.el2 = ( (1 - (1/(beta.prime/c2 + 1)))^alpha.prime ) *
    ( (beta.prime/c2 + 1)^(-x.non) ) * p2.prime

  q1 = q1.el1/(q1.el1 + q1.el2)
  q2 = 1 - q1


  # expected log likelihood of x.syn

  ll.x.syn = sum(log(dnegbin(x.syn, alpha, 1/(beta + 1))))


  # expected log likelihood of x.non

  summand1 = log(dnegbin(x.non, alpha, 1/(beta/c1 + 1))*p1)*q1
  summand2 = log(dnegbin(x.non, alpha, 1/(beta/c2 + 1))*p2)*q2

  ll.x.non = sum(summand1 + summand2)


  # Q

  Q = ll.x.syn + ll.x.non

  return(Q)
}


negQ = function(theta, theta.prime, x.syn, x.non, c1, c2) {
  -Q(theta, theta.prime, x.syn, x.non, c1, c2)
}