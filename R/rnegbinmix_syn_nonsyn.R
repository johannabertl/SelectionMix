#' Simulating synonymous and non-synonymous mutation counts
#'
#' rnegbinmix_syn_nonsyn simulates counts of synonymous and non-synonymous mutations, using a negative binomial and a negative binomial mixture, respectively, with the same parameters alpha and beta.
#'
#' The function \link{MASS::rnegbin} from package MASS is used for the simulation of negative binomial data. The function \link{rnegbinmix} is used for the negative binomial mixture.
#'
#' @param n number of genes
#' @param alpha scalar, positive
#' @param beta scalar, positive
#' @param c numeric, positive, length k
#' @param p numeric, probability vector, length k
#'
#'
#' @author Johanna Bertl
#'
#' @examples
#'
#' # Mixture with 3 components
#'
#' c1 = 0.01; c2 = 1; c3=100
#' p1 = 0.33; p2 = 0.33; p3 = 0.34
#' alpha = 10
#' beta = 5
#'
#' simulated.mixture = rnegbinmix_syn_nonsyn(1000, alpha = alpha, beta = beta, c = c(c1, c2, c3), p = c(p1, p2, p3))
#'
#' @export

rnegbinmix_syn_nonsyn = function(n, alpha, beta, c, p){
  Syn = rnegbin.alphabeta(n, alpha = alpha, beta = beta)
  Non = rnegbinmix(n, alpha = alpha, beta=beta, c = c, p = p)

  data.frame(Syn = Syn, Non = Non)
}


