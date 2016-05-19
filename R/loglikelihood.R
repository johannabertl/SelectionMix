#' The (incomplete data) loglikelihood function for the negative binomial mixture model
#'
#'
#' @param theta c(alpha, beta, p1, p2, ..., pk-1)
#' @param x.syn vector of synonymous mutation count
#' @param x.non vector of non-synonymous mutation count
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
#' mutations = rnegbinmix_syn_nonsyn(n=3000, alpha=alpha, beta=beta, c=c(c1, c2), p = c(p1, p2))
#' # log likelihood at the true value
#' loglikelihood(theta=c(alpha,beta,p1), x.syn = mutations$Syn, x.non = mutations$Non, cvec=c(c1, c2))
#'
#' # likelihood surface
#'
#' alphavec = 1:20
#' betavec = seq(0.5, 10, by=0.5)
#' p1vec = seq(0.11, 0.3, by=0.01)
#'
#' llalpha = numeric(20)
#' llbeta = numeric(20)
#' llp1 = numeric(20)
#'
#' for(i in 1:20){
#'  llalpha[i] = loglikelihood(theta=c(alphavec[i],beta,p1), x.syn = mutations$Syn, x.non = mutations$Non, cvec=c(c1, c2))
#'  llbeta[i] = loglikelihood(theta=c(alpha,betavec[i],p1), x.syn = mutations$Syn, x.non = mutations$Non, cvec=c(c1, c2))
#'  llp1[i] = loglikelihood(theta=c(alpha,beta,p1vec[i]), x.syn = mutations$Syn, x.non = mutations$Non, cvec=c(c1, c2))
#' }
#'
#' plot(alphavec, llalpha, t="b")
#' abline(v=alpha, col="red")
#' plot(betavec, llbeta, t="b")
#' abline(v=beta, col="red")
#' plot(p1vec, llp1, t="b")
#' abline(v=p1, col="red")

loglikelihood = function(theta, x.syn, x.non, cvec){

  alpha = theta[1]
  beta = theta[2]
  p = theta[3:length(theta)]
  p = c(p, 1-sum(p))

  k = length(p)
  n = length(x.syn)

  ll.x.syn = sum(log(dnegbin(x.syn, alpha, 1/(beta + 1))))

  ll.x.non.mat = matrix(NA, ncol=k, nrow=n)
  for(i in 1:k){
    ll.x.non.mat[,i] = p[i]*dnegbin(x.non, alpha, 1/(beta/cvec[i] + 1))
  }

  ll.x.non = sum(log( rowSums(ll.x.non.mat) ))

  ll.x.syn + ll.x.non
}
