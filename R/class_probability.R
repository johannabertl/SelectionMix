#' Posterior class probabilities for a gene
#'
#' The probability that a gene with a given number of non-synonymous mutations is member of each selection class is estimated given a set of (estimated) parameters.
#'
#' @param x.non Numeric. Number of non-synonymous mutations. Can be a scalar (i. e. one gene) or a vector (a set of genes).
#' @param theta Numeric. Model parameters, c(alpha, beta, p1, ..., pk) such that p1 + ... + p_k = 1
#' @param cvec Numeric. c(c_1, ..., c_k)
#'
#' @author Johanna Bertl
#'
#' @examples
#'
#' # Usage on 3 genes
#' # using parameters alpha = 1, beta = 1, p1=p2=0.5, c1 = 1, c2 = 10
#' class_probabilities(c(2, 7, 100), theta = c(1, 1, 0.5, 0.5), c(1, 10))
#'
#' # All genes in the yeast data:
#' data("yeast")
#'
#' cvec = c(0.1, 1, 10, 100)
#' p = c(0.19, 0.32, 0.48, 0.007)
#' alpha = 0.9
#' beta = 45
#'
#' # class probabilities for all genes:
#' prob_mat = class_probabilities(yeast$Non, theta = c(alpha, beta, p), cvec = cvec)
#'
#' # class with the highest probability:
#' prob_max = apply(prob_mat, 1, which.max)
#'
#' @export
class_probabilities = function(x.non, theta, cvec){

  if(length(x.non)>1){
    class_probabilities_vector(x.non, theta, cvec)
  } else {
    class_probabilities_scalar(x.non, theta, cvec)
  }

}

class_probabilities_scalar = function(x.non, theta, cvec){

  alpha = theta[1]
  beta = theta[2]
  p = theta[3:length(theta)]

  if(!all(p>0)){
    stop("All entries of p must be positive.")
  }
  if(length(p)!=length(cvec)){
    stop("p and c must have the same length.")
  }
  if(sum(p) !=1 ) {
    p = p/sum(p)
    warning("p has been scaled such that sum(p)=1.")
  }

  k = length(p)

  dnegbin.alphabeta(x.non, alpha = alpha, beta = beta/cvec) * p / dnegbinmix(x.non, alpha = alpha, beta = beta, cvec = cvec, p = p)

}

class_probabilities_vector = function (x.non, theta, cvec) {

  t(sapply(X = x.non, FUN = class_probabilities_scalar, theta = theta, cvec = cvec, simplify=T))
}
