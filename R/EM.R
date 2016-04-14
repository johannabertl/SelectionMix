#' EM algorithm for the negative binomial mixture model
#'
#' EM algorithm to estimate the parameters of the negative binomial mixture model. Currently only for mixtures of 2 negbinoms as defined in \link{negbinmix::rnegbinmix}.
#'
#' In the M step, contrained optimization is conducted with the constrOptim function. ui and ci are defined accordingly (see \link{constrOptim}): ui %*% theta - ci >= 0.
#'
#' @param theta.null
#' @param x.syn
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

EM = function(theta.null, x.syn, x.non, c1, c2, iter, ui, ci){

  theta.mat = rbind(theta.null, matrix(NA, ncol=3, nrow=iter))
  Qvec = numeric(iter+1)
  Qvec[1] = Q(theta.null, theta.null, x.syn, x.non, c1, c2)
  convergence = numeric(iter+1)
  convergence[1] = NA

  for(i in 1:iter){
    optres = constrOptim(theta = theta.null, f = negQ, grad=NULL, ui=ui, ci=ci, theta.prime = theta.mat[i,], x.syn=x.syn, x.non=x.non, c1=c1, c2=c2)

    theta.mat[i+1,] = optres$par
    Qvec[i+1] = optres$value
    convergence[i+1] = optres$convergence
  }

  return(list(theta = theta.mat, Q = Qvec, convergence = convergence))

}
