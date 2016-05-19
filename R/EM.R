#' EM algorithm for the negative binomial mixture model
#'
#' EM algorithm to estimate the parameters of a negative binomial mixture model with k components.
#'
#' The parameters alpha, beta and p1, ..., pk-1 are estimated. pk = 1 - p1 - ... -pk-1. c1, ..., ck are fixed.
#'
#' In the M step, constrained optimization is conducted with the function \code{constrOptim}. \code{ui} and \code{ci} are defined accordingly (see \link{constrOptim}): \code{ui \%*\% theta - ci >= 0}. It is important that ui and ci are specified such that p is a probability vector, see examples.
#'
#' As a stopping criterion, the difference in log-likelihood can be used. To make sure that there is always at least a short trajectory for visual checks, the criterion is only used after at least 5 iterations.
#'
#' @param theta.null Starting value, vector consisting of alpha, beta, p1, ..., pk-1.
#' @param x.syn Vector of synonymous mutation counts per gene.
#' @param x.non Vector of non-synonymous mutation counts per gene.
#' @param cvec Parameter vector c1, ..., ck.
#' @param iter Maximum number of iterations.
#' @param epsilon Log-likelihood difference stopping criterion. If NULL, the algorithm is run for \code{iter} iterations. Otherwise, it is run until the increase in the log-likelihood is smaller than \code{epsilon}, but at most for \code{iter} iterations.
#' @param ui Linear constraints for theta. See details and examples.
#' @param ci Linear constraints for theta. See details and examples.

#'
#' @author Johanna Bertl
#'
#' @return
#' \describe{
#'  \item{theta}{Matrix with the trajectories of the parameter values.}
#'  \item{Q}{Vector with the trajectory of the Q values.}
#'  \item{loglikelihood}{Vector with the trajectory of the log likelihood values.}
#'  \item{convergence}{Vector with the convergence code for each iteration of the M-step of the algorithm, see \code{\link{stats::optim}} (0 indicates successful completion).}
#' }
#'
#' @seealso
#'
#' \code{\link{rnegbinmix::EM_random}} for a wrapper function with randomly drawn starting values.
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
#' # function Q
#' Qtest = Q(theta = c(5, 1, 0.5), theta.prime = c(5, 1, 0.5), x.syn = mutations$Syn, x.non = mutations$Non, cvec = c(c1, c2))
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
#' EMres = EM(theta.null = c(5, 1, 0.5), x.syn = mutations$Syn, x.non = mutations$Non, cvec=c(c1, c2), iter=25, epsilon = 10^(-8), ui = ui, ci=ci)
#' matplot(EMres$theta, t="b")
#' plot(EMres$loglikelihood, t="b")
#'
#'
#' ### Example with 3 components ####
#'
#' c1 = 0.01; c2 = 1; c3 = 100
#' p1 = 0.33; p2 = 0.33; p3 = 0.34
#' alpha = 10
#' beta = 5
#'
#' mutations = rnegbinmix_syn_nonsyn(n=3000, alpha=alpha, beta=beta, c=c(c1, c2, c3), p = c(p1, p2, p3))
#'
#' # Q
#'
#' Qtest = Q(theta = c(5, 1, 0.33, 0.33), theta.prime = c(5, 1, 0.33, 0.33), x.syn = mutations$Syn, x.non = mutations$Non, c = c(c1, c2, c3))
#'
#' # EM algorithm
#'
#' # Constraints: ui %*% theta - ci >= 0
#' # alpha > 0, beta > 0, 0 <= p1, p2 <=1, p1 + p2 <= 1
#'
#' ui = rbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1),
#'  c(0, 0, -1, 0), c(0, 0, 0, -1), c(0, 0, -1, -1))
#' ci = c(0, 0, 0, 0, -1, -1, -1)
#'
#' EMres = EM(theta.null = c(5, 1, 0.33, 0.33), x.syn = mutations$Syn, x.non = mutations$Non, cvec=c(c1, c2, c3), iter=25, epsilon = 10^(-8), ui = ui, ci=ci)
#' matplot(EMres$theta, t="b")
#' plot(EMres$loglikelihood, t="b")


EM = function(theta.null, x.syn, x.non, cvec, iter, epsilon=NULL, ui, ci){

  # prepare the data structures

  theta.mat = rbind(theta.null, matrix(NA, ncol=length(theta.null), nrow=iter))
  Qvec = numeric(iter+1)
  Qvec[1] = Q(theta.null, theta.null, x.syn, x.non, cvec)
  convergence = numeric(iter+1)
  convergence[1] = NA
  loglik = numeric(iter+1)
  loglik[1] = loglikelihood(theta.null, x.syn, x.non, cvec)


  for(i in 1:iter){
    # E- and M-step
    optres = constrOptim(theta = theta.mat[i,], f = negQ, grad=NULL, ui=ui, ci=ci, theta.prime = theta.mat[i,], x.syn=x.syn, x.non=x.non, cvec = cvec)

    # save the results
    theta.mat[i+1,] = optres$par
    Qvec[i+1] = -optres$value
    convergence[i+1] = optres$convergence
    loglik[i+1] = loglikelihood(theta.mat[i+1,], x.syn, x.non, cvec)

    # stop if the difference in likelihood estimates is smaller than epsilon (but not earlier than after 5 iterations!)
    stopped.early=F
    if(!is.null(epsilon)){
      if(loglik[i+1] - loglik[i] < epsilon) {
        print(paste0("The convergence criterion is met in iteration ",i, "."))
        stopped.early=T
        break
      }
    }
  }

  # if the run has been stopped early, truncate the output matrices and vector
  if(stopped.early){
    theta.mat = theta.mat[1:(i+1),]
    Qvec = Qvec[1:(i+1)]
    convergence = convergence[1:(i+1)]
    loglik = loglik[1:(i+1)]
  }

  return(list(theta = theta.mat, Q = Qvec, convergence = convergence, loglikelihood = loglik))

}
