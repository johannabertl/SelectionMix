#' Wrapper function with random starting values for the EM algorithm
#'
#' The function EM is called multiple times, each time with a random starting value
#'
#' The starting values for the parameters alpha and beta are drawn from a uniform distribution on two user-defined intervals. The starting values for p1, ..., pk-1 are drawn from a multinomial distribution. The variance of the multinomial distribution is determined by the parameter \code{precision}: a multinomial sample of the size \code{precision} is drawn with the function \code{rmultinom}, afterwards the multinomial sample is scaled to sum up to one. The larger \code{precision}, the smaller the variance.
#'
#' @param alpha.interval
#' @param beta.interval
#' @param pvec vector of length k-1
#' @param precision positive integer to determine the variance of the multinomial draws (see details).
#' @param nstart Integer. Number of random starting values
#' @param parallel Logical. If true, the function mclapply is used for parallelization. Only makes sense when multiple cores are available.
#'
#' The other parameters are the same as in \code{\link{rnegbinmix::EM}}.
#'
#' @return
#'
#' A list consisting of \code{nstart} objects, each one a list just like the output of EM (trajectory of theta, Q and convergence code), see \code{\link{rnegbinmix::EM}}.
#'
#' @author Johanna Bertl
#'
#' @examples
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
#' # EM algorithm
#'
#' # Constraints: ui %*% theta - ci >= 0
#' # alpha > 0, beta > 0, 0 <= p1 <=1
#'
#' ui = rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(0, 0, -1))
#' ci = c(0, 0, 0, -1)
#'
#' # Starting values
#' alpha.interval = c(1, 20)
#' beta.interval = c(1, 10)
#' pvec = p1
#' precision=20
#' nstart = 5
#'
#' # Run EM algorithm
#' EMres = EM_random(alpha.interval = alpha.interval, beta.interval = beta.interval, pvec = pvec, precision = precision, nstart=nstart, x.syn = mutations$Syn, x.non = mutations$Non, cvec=c(c1, c2), iter=25, ui = ui, ci=ci, parallel=T)
#'
#' # Plot results
#' par(mfrow=c(3,2))
#'  for(i in 1:5) matplot(EMres[[i]]$theta, t="l")
#' par(mfrow=c(1,1))

EM_random = function(alpha.interval, beta.interval, pvec, precision, nstart, x.syn, x.non, cvec, iter, epsilon=NULL, ui, ci, parallel=F) {

  # prepare matrix for starting values
  k = length(pvec) + 3
  start = matrix(NA, nrow=nstart, ncol=k)

  # make sure pvec is a valid probability vector
  if(sum(pvec)>=1 | !(all(pvec>0))){
    stop("pvec is not a valid probability vector.")
  }
  #!# further checks

  # simulate starting values
  # alpha:
  start[,1] = runif(nstart, min = alpha.interval[1], max = alpha.interval[2])
  # beta:
  start[,2] = runif(nstart, min = beta.interval[1], max = beta.interval[2])
  # p:
  start[,3:k] = t(rmultinom(nstart, size=precision, prob = c(pvec, 1-sum(pvec)))/precision)


  if(parallel){
    startlist = vector("list", nstart)
    for(i in 1:nstart) startlist[[i]] = start[i,1:(k-1)]
    mclapply(startlist, EM, x.syn, x.non, cvec, iter, epsilon, ui, ci)
  } else {
    apply(start, 1, EM, x.syn, x.non, cvec, iter, epsilon, ui, ci)
  }
}
