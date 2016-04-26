#' Wrapper function with random starting values for the EM algorithm
#'
#' The function EM is called multiple times, each time with a random starting value drawn from a uniform distribution on a user-defined hypercube.
#'
#' @param theta.null 3x2 matrix, with minimal value in column one, maximal value in column 2, alpha, beta, p1 in lines 1-3.
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
#' theta.null = matrix(c(1, 1, 0.05, 20, 10, 0.5), ncol=2, nrow=3, byrow=F)
#' nstart = 5
#'
#' # Run EM algorithm
#' EMres = EM_random(theta.null = theta.null, nstart=nstart, x.syn = mutations$Syn, x.non = mutations$Non, c1 = c1, c2 = c2, iter=25, ui = ui, ci=ci, parallel=T)
#'
#' # Plot results
#' par(mfrow=c(3,2))
#'  for(i in 1:5) matplot(EMres[[i]]$theta, t="l")
#' par(mfrow=c(1,1))

EM_random = function(theta.null, nstart, x.syn, x.non, c1, c2, iter, ui, ci, parallel=F) {

  start = matrix(NA, nrow=nstart, ncol=3)
  for(i in 1:3) start[,i] = runif(nstart, min = theta.null[i,1], max = theta.null[i, 2])

  if(parallel){
    startlist = vector("list", nstart)
    for(i in 1:nstart) startlist[[i]] = start[i,]
    mclapply(startlist, EM, x.syn, x.non, c1, c2, iter, ui, ci)
  } else {
    apply(start, 1, EM, x.syn, x.non, c1, c2, iter, ui, ci)
  }
}
