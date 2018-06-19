#' EM algorithm for the negative binomial mixture model with smart starting values
#'
#' EM algorithm to estimate the parameters alpha, beta, p1, ..., pk-1 of a negative binomial mixture model with k components.
#'
#' The parameters alpha and beta are estimated on the synonymous mutations by ML. Then, the parameters p1, ..., pk-1 are estimated by an EM algorithm on the non-synonymous mutations using the estimated alpha and beta as known parameters (EM1). Finally, alpha, beta and p1, ..., pk-1 are used as starting values in another EM algorithm, where all of them are re-estimated (EM2).
#'
#' For the estimation of alpha and beta on the synonymous mutations, the function \code{MASS::fitdistr} is used.
#'
#' In the M step of both EM algorithms, constrained optimization is conducted with the function \code{constrOptim}. \code{ui} and \code{ci} are defined accordingly (see \link{constrOptim}): \code{ui \%*\% theta - ci >= 0}. It is important that ui and ci are specified such that p is a probability vector, see examples. It is not allowed to specify constraints that affect (alpha and beta) and (p1, ..., pk-1) jointly.
#'
#' As a stopping criterion, the difference in log-likelihood can be used. To make sure that there is always at least a short trajectory for visual checks, the criterion is only used after at least 5 iterations.
#'
#' The function EM_smartstart_apply can be used to run the algorithm on a list of datasets (e. g. using mclapply for parallel execution).
#'
#' @param theta.null Starting values for the parameters p1, ..., pk-1 in EM1.
#' @param x.syn Vector of synonymous mutation counts per gene.
#' @param x.non Vector of non-synonymous mutation counts per gene.
#' @param cvec Parameters c1, ..., ck.
#' @param iter1 Maximum number of iterations for EM1.
#' @param iter2 Maximum number of iterations for EM2.
#' @param epsilon1 Stopping criterion for EM1. If NULL, the maximum number of iterations is used.
#' @param epsilon2 Stopping criterion for EM2. If NULL, the maximum number of iterations is used.
#' @param ui Constraints on the parameters. See details.
#' @param ci Constraints on the parameters. See details.
#' @param full.output logical. If T, the results of the ML estimation and EM1 are output. Otherwise, only the output of EM2 is output (in the same format as by the function \code{\link{EM}}).
#'
#' @seealso
#'
#' MASS::fitdistr
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
#' # EM algorithm
#'
#' # Constraints: ui %*% theta - ci >= 0
#' # alpha > 0, beta > 0, 0 <= p1 <=1
#'
#' ui = rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(0, 0, -1))
#' ci = c(0, 0, 0, -1)
#'
#' # EM algorithm
#' EMtest = EM_smartstart(theta.null = 0.5, x.syn = mutations$Syn, x.non = mutations$Non, cvec=c(c1, c2), iter1=10, iter2=10, epsilon2 = 10^(-8), ui = ui, ci = ci)
#' matplot(EMtest$theta, t="b")
#' plot(EMtest$loglikelihood, t="b")
#'
#'
#'
#' ### Example with 4 components ###
#'
#' c1 = 0.1; c2 = 1; c3 = 10; c4 = 100;
#' cvec = c(c1, c2, c3, c4)
#' p1 = 0.2; p2 = 0.6; p3 = 0.1; p4 = 0.1;
#' pvec = c(p1, p2, p3, p4)
#' alpha = 1.72396
#' beta = 89.019
#'
#' mutations = rnegbinmix_syn_nonsyn(n=3000, alpha=alpha, beta=beta, c=cvec, p = pvec)
#'
#'
#' # EM algorithm
#'
#' # Constraints: ui %*% theta - ci >= 0
#' # alpha > 0, beta > 0, 0 <= p1 <=1
#'
#' ui = rbind(diag(1, 5),
#'      cbind(rep(0, 3), rep(0, 3), diag(-1, 3)),
#'      c(0, 0, rep(-1, 3)))
#'
#' ci = c(rep(0, 5), rep(-1, 4))
#'
#' # EM algorithm
#' EMtest = EM_smartstart(theta.null = c(0.25, 0.25, 0.25), x.syn = mutations$Syn, x.non = mutations$Non, cvec=cvec, iter1=10, iter2=50, epsilon2 = 10^(-8), ui = ui, ci=ci, full.output=T)
#' matplot(EMtest$EM_final$theta, t="b")
#' plot(EMtest$EM_final$loglikelihood, t="b")
#'
#'
#' @export

EM_smartstart = function(theta.null, x.syn, x.non, cvec, iter1, iter2, epsilon1=NULL, epsilon2=NULL, ui, ci, full.output=F){

  ### estimating alpha and beta from the synonymous mutations: ###

  Syn.negbin.fit = fitdistr(x.syn, "negative binomial")
  alpha.start = Syn.negbin.fit$estimate["size"]
  beta.start = Syn.negbin.fit$estimate["size"]/Syn.negbin.fit$estimate["mu"]

  ### Obtaining the constraints for ui and ci for p1, ..., pk-1 only

  # check if there are joint constraints on (alpha, beta) and (p1, ..., pk-1)
  constr1 = rowSums(abs(ui[,1:2]))>0 # constraints on (alpha, beta)
  constr2 = rowSums(as.matrix(abs(ui[,3:(ncol(ui))])))>0 # constraints on (p1, ..., pk-1)
  if(sum(constr1 & constr2) > 0) stop("No joint constraints on (alpha, beta) and (p1, ..., pk-1) allowed.")

  # remove the first 2 columns (they only concern alpha and beta)
  ui1 = as.matrix(ui[,-c(1,2)])
  # remove the empty lines
  empty = which(rowSums(abs(ui1))==0)
  ui1 = as.matrix(ui1[-empty,])
  ci1 = ci[-empty]

  ### estimating p1, ..., pk-1 from the nonsynonsymous mutations ###

  EM_start = EM_proportions(theta.null = theta.null, x.non = x.non, alpha=alpha.start, beta = beta.start, cvec=cvec, iter=iter1, epsilon = epsilon1, ui = ui1, ci=ci1)
  #!# avoid one-dimensional optimization by Nelder-Mead (warning says "Brent" is more reliable)
  p.start = as.matrix(EM_start$theta)[length(EM_start$Q),]

  ### starting values ###
  start = c(alpha.start, beta.start, p.start)

  ### EM algorithm to estimate the full set of parameters ###
  EM_final = EM(theta.null = start, x.syn=x.syn, x.non=x.non, cvec=cvec, iter=iter2, epsilon=epsilon2, ui=ui, ci=ci)

  ### output ###
  if(full.output){
    return(list(EM_start = EM_start, start = start, EM_final = EM_final))
  } else {
    return(EM_final)
  }

}


EM_smartstart_apply = function(x, theta.null, cvec, iter1, iter2, epsilon1, epsilon2, ui, ci, full.output){
  EM_smartstart(theta.null=theta.null, x.syn=x$Syn, x.non=x$Non, cvec=cvec, iter1=iter1, iter2=iter2, epsilon1=epsilon1, epsilon2=epsilon2, ui=ui, ci=ci, full.output=full.output)
}
