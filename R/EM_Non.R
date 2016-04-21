#' EM algorithm for the negative binomial mixture model on the non-synonymous mutations only
#'
#' EM algorithm to estimate the parameters of the negative binomial mixture model. Currently only for mixtures of 2 negbinoms as defined in \link{negbinmix::rnegbinmix}.
#'
#' The parameters alpha, beta and p1 are estimated. p2 = 1-p2. c1 and c2 are fixed.
#'
#' In the M step, contrained optimization is conducted with the constrOptim function. ui and ci are defined accordingly (see \link{constrOptim}): ui %*% theta - ci >= 0.
#'
#' @param theta.null
#' @param x.non
#' @param c1
#' @param c2
#' @param iter
#' @param ui
#' @param ci
#'
#' @author Johanna Bertl
#'
#' @return
#' \describe{
#'  \item{theta}{Matrix with the trajectories of the parameter values}
#'  \item{Q}{Vector with the trajectory of the Q values}
#'  \item{convergence}{Vector with the convergence code for each iteration of the M-part of the algorithm, see \link{stats::optim} (0 indicates successful completion).}
#' }
#'
#' @examples
#'
#' # simulating nonsynonymous mutations
#'
#' ngenes = 3408
#' c1 = 0.5; c2 = 10
#' p1 = 0.2; p2 = 0.8
#' alpha = 10
#' beta = 5
#'
#' mutations = rnegbinmix(n=ngenes, alpha, beta, c = c(c1, c2), p = c(p1, p2))
#' plot(table(mutations))
#'
#' # EM algorithm
#'
#' # Constraints: ui %*% theta - ci >= 0
#' ui = rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(0, 0, -1))
#' ci = c(0, 0, 0, -1)
#'
#' EM_output = EM_Non(c(5, 1, 0.5), x.non = mutations, c1 = c1, c2 = c2, iter=100, ui=ui, ci=ci)
#'
#' # convergence in each iteration?
#' sum(EM_output$convergence, na.rm=T)
#' # yes.
#'
#' par(mfrow=c(4,1))
#' plot(EM_output$theta[,1], t="b", main="alpha")
#' plot(EM_output$theta[,2], t="b", main="beta")
#' plot(EM_output$theta[,3], t="b", main="p1")
#' plot(EM_output$Q, t="b", main="Q")
#' par(mfrow=c(1,1))

EM_Non = function(theta.null, x.non, c1, c2, iter, ui, ci){

  theta.mat = rbind(theta.null, matrix(NA, ncol=3, nrow=iter))
  Qvec = numeric(iter+1)
  Qvec[1] = QNon(theta.null, theta.null, x.non, c1, c2)
  convergence = numeric(iter+1)
  convergence[1] = NA

  for(i in 1:iter){
    optres = constrOptim(theta = theta.mat[i,], f = negQNon, grad=NULL, ui=ui, ci=ci, theta.prime = theta.mat[i,], x.non=x.non, c1=c1, c2=c2)

    theta.mat[i+1,] = optres$par
    Qvec[i+1] = -optres$value
    convergence[i+1] = optres$convergence
  }

  return(list(theta = theta.mat, Q = Qvec, convergence = convergence))

}
