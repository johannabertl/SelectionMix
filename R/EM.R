#' EM algorithm for the negative binomial mixture model
#'
#' EM algorithm to estimate the parameters of the negative binomial mixture model. Currently only for mixtures of 2 negbinoms as defined in \link{negbinmix::rnegbinmix}.
#'
#' The parameters alpha, beta and p1 are estimated. p2 = 1-p2. c1 and c2 are fixed.
#'
#' In the M step, contrained optimization is conducted with the constrOptim function. ui and ci are defined accordingly (see \link{constrOptim}): ui %*% theta - ci >= 0.
#'
#' @param theta.null Starting value, vector consisting of alpha, beta, p1.
#' @param x.syn Vector of synonymous mutation counts per gene.
#' @param x.non Vector of non-synonymous mutation counts per gene.
#' @param c1 Parameter c1.
#' @param c2 Parameter c2.
#' @param iter Number of iterations
#' @param ui Linear constraints for theta. See details and examples.
#' @param ci Linear constraints for theta. See details and examples.
#'
#' @author Johanna Bertl
#'
#' @return
#' \describe{
#'  \item{theta}{Matrix with the trajectories of the parameter values}
#'  \item{Q}{Vector with the trajectory of the Q values}
#'  \item{convergence}{Vector with the convergence code for each iteration of the M-step of the algorithm, see \code{\link{stats::optim}} (0 indicates successful completion).}
#' }
#'
#' @seealso
#'
#' \code{\link{rnegbinmix::EM_random}} for a wrapper function with randomly drawn starting values.
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
#' # EM algorithm
#' EMres = EM(theta.null = c(5, 1, 0.5), x.syn = mutations$Syn, x.non = mutations$Non, c1 = c1, c2 = c2, iter=25, ui = ui, ci=ci)
#' matplot(EMres$theta)


EM = function(theta.null, x.syn, x.non, c1, c2, iter, ui, ci){

  theta.mat = rbind(theta.null, matrix(NA, ncol=3, nrow=iter))
  Qvec = numeric(iter+1)
  Qvec[1] = Q(theta.null, theta.null, x.syn, x.non, c1, c2)
  convergence = numeric(iter+1)
  convergence[1] = NA

  for(i in 1:iter){
    optres = constrOptim(theta = theta.mat[i,], f = negQ, grad=NULL, ui=ui, ci=ci, theta.prime = theta.mat[i,], x.syn=x.syn, x.non=x.non, c1=c1, c2=c2)

    theta.mat[i+1,] = optres$par
    Qvec[i+1] = -optres$value
    convergence[i+1] = optres$convergence
  }

  return(list(theta = theta.mat, Q = Qvec, convergence = convergence))

}
