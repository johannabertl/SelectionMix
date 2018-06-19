#' The (incomplete data) loglikelihood function for the parameters p1, ..., pk-1 of the negative binomial mixture model
#'
#' The parameters alpha and beta are known.
#'
#' @param theta c(p1, p2, ..., pk-1)
#' @param x.non vector of non-synonymous mutation count
#' @param alpha
#' @param beta
#' @param cvec vector of k positive values c(c1, ..., ck)
#'
#' @author Johanna Bertl
#'
#' @examples
#'
#' c1 = 0.5; c2 = 10
#' p1 = 0.2; p2 = 0.8
#' alpha = 10
#' beta = 5
#'
#' mutations = rnegbinmix(n=3000, alpha=alpha, beta=beta, c=c(c1, c2), p = c(p1, p2))
#'
#' # likelihood surface
#'
#' p1vec = seq(0.11, 0.3, by=0.01)
#'
#' llp1 = numeric(20)
#'
#' for(i in 1:20){
#'  llp1[i] = loglikelihood_proportions(theta=p1vec[i], x.non = mutations, alpha = alpha, beta = beta, cvec=c(c1, c2))
#' }
#'
#' plot(p1vec, llp1, t="b")
#' abline(v=p1, col="red")
#'
#' @export

loglikelihood_proportions = function(theta, x.non, alpha, beta, cvec){

  p = c(theta, 1-sum(theta))

  k = length(p)
  n = length(x.non)

  ll.x.non.mat = matrix(NA, ncol=k, nrow=n)
  for(i in 1:k){
    ll.x.non.mat[,i] = p[i]*dnegbin(x.non, alpha, 1/(beta/cvec[i] + 1))
  }

  sum(log( rowSums(ll.x.non.mat) ))

}
