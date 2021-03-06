#' EM algorithm for the negative binomial mixture model
#'
#' EM algorithm to estimate the parameters p1, ..., pk-1 of a negative binomial mixture model with k components.
#'
#' The parameters alpha, beta, c1, ..., ck are fixed. Alpha and beta can for example be obtained from the synonymous mutations.
#'
#' In the M step, constrained optimization is conducted with the function \code{constrOptim}. \code{ui} and \code{ci} are defined accordingly (see \link{constrOptim}): \code{ui \%*\% theta - ci >= 0}. It is important that ui and ci are specified such that p is a probability vector, see examples.
#'
#' As a stopping criterion, the difference in log-likelihood can be used. To make sure that there is always at least a short trajectory for visual checks, the criterion is only used after at least 5 iterations.
#'
#' @param theta.null Starting value, vector consisting of p1, ..., pk-1.
#' @param x.non Vector of non-synonymous mutation counts per gene.
#' @param alpha
#' @param beta
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
#' \code{\link{rnegbinmix::EM}} for estimation of all parameters.
#'
#' @examples
#'
#' ### Example with 2 components ###
#'
#' # Simulate dataset non-synonymous mutations
#'
#' c1 = 0.5; c2 = 10
#' p1 = 0.2; p2 = 0.8
#' alpha = 10
#' beta = 5
#'
#' mutations = rnegbinmix(n=3000, alpha=alpha, beta=beta, c=c(c1, c2), p = c(p1, p2))
#'
#' # EM algorithm
#'
#' # Constraints: ui %*% theta - ci >= 0
#' # 0 <= p1 <=1
#'
#' ui = matrix(c(1, -1), ncol=1)
#' ci = c(0, -1)
#'
#' # EM algorithm
#' EMres = EM_proportions(theta.null = 0.5, x.non = mutations, alpha = alpha, beta = beta, cvec=c(c1, c2), iter=25, epsilon = 10^(-8), ui = ui, ci=ci)
#' plot(EMres$theta, t="b")
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
#' # Simulate mutations data (synonymous and non-synonymous)
#'
#' mutations = rnegbinmix_syn_nonsyn(n=3000, alpha=alpha, beta=beta, c=c(c1, c2, c3), p = c(p1, p2, p3))
#'
#' # Q
#'
#' Qall = Q(theta = c(5, 1, 0.33, 0.33), theta.prime = c(5, 1, 0.33, 0.33), x.syn = mutations$Syn, x.non = mutations$Non, c = c(c1, c2, c3))
#' Qnon = Q_proportions(theta = c(0.33, 0.33), theta.prime = c(0.33, 0.33), x.non = mutations$Non, alpha=10, beta=5, c = c(c1, c2, c3))
#'
#' # EM algorithm for the proportions
#'
#' #' # Constraints: ui %*% theta - ci >= 0
#' # 0 <= p1, p2 <=1, p1 + p2 <= 1
#'
#' ui = rbind(c(1, 0), c(0, 1),
#'  c(-1, 0), c( 0, -1), c(-1, -1))
#' ci = c(0, 0, -1, -1, -1)
#'
#' EM_prop = EM_proportions(theta.null = c(0.33, 0.33), x.non = mutations$Non, alpha=alpha, beta = beta, cvec=c(c1, c2, c3), iter=25, epsilon = 10^(-8), ui = ui, ci=ci)
#' matplot(EM_prop$theta, t="b")
#' plot(EM_prop$loglikelihood, t="b")
#'
#' # EM algorithm for all parameters
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
#'
#'
#' @export


EM_proportions = function(theta.null, x.non, alpha, beta, cvec, iter, epsilon=NULL, ui, ci){

  # prepare the data structures

  theta.mat = rbind(theta.null, matrix(NA, ncol=length(theta.null), nrow=iter))
  Qvec = numeric(iter+1)
  Qvec[1] = Q_proportions(theta.null, theta.null, x.non, alpha, beta, cvec)
  convergence = numeric(iter+1)
  convergence[1] = NA
  loglik = numeric(iter+1)
  loglik[1] = loglikelihood_proportions(theta.null, x.non, alpha, beta, cvec)


  for(i in 1:iter){
    # E- and M-step
    optres = constrOptim(theta = theta.mat[i,], f = negQ_proportions, grad=NULL, ui=ui, ci=ci, theta.prime = theta.mat[i,], x.non=x.non, alpha=alpha, beta=beta, cvec=cvec)

    # save the results
    theta.mat[i+1,] = optres$par
    Qvec[i+1] = -optres$value
    convergence[i+1] = optres$convergence
    loglik[i+1] = loglikelihood_proportions(theta.mat[i+1,], x.non, alpha, beta, cvec)

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
