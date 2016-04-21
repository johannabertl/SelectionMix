#' Simulating synonymous and non-synonymous mutation counts
#'
#' rnegbinmix_syn_nonsyn simulates counts of synonymous and non-synonymous mutations, using a negative binomial and a negative binomial mixture, respectively, with the same parameters alpha and beta.
#'
#' The function \link{MASS::rnegbin} from package MASS is used for the simulation of negative binomial data. The function \link{rnegbinmix} is used for the negative binomial mixture.
#'
#' Currently, only a mixture of two components is implemented for the non-synonymous mutations.
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

rnegbinmix_syn_nonsyn = function(n, alpha, beta, c, p){
  Syn = rnegbin.alphabeta(n, alpha = alpha, beta = beta)
  Non = rnegbinmix(n, alpha = alpha, beta=beta, c = c(c1, c2), p = c(p1, p2))

  data.frame(Syn = Syn, Non = Non)
}


