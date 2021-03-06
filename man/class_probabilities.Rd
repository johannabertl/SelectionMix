% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_probability.R
\name{class_probabilities}
\alias{class_probabilities}
\title{Posterior class probabilities for a gene}
\usage{
class_probabilities(x.non, theta, cvec)
}
\arguments{
\item{x.non}{Numeric. Number of non-synonymous mutations. Can be a scalar (i. e. one gene) or a vector (a set of genes).}

\item{theta}{Numeric. Model parameters, c(alpha, beta, p1, ..., pk) such that p1 + ... + p_k = 1}

\item{cvec}{Numeric. c(c_1, ..., c_k)}
}
\description{
The probability that a gene with a given number of non-synonymous mutations is member of each selection class is estimated given a set of (estimated) parameters.
}
\examples{

# Usage on 3 genes
# using parameters alpha = 1, beta = 1, p1=p2=0.5, c1 = 1, c2 = 10
class_probabilities(c(2, 7, 100), theta = c(1, 1, 0.5, 0.5), c(1, 10))

# All genes in the yeast data:
data("yeast")

cvec = c(0.1, 1, 10, 100)
p = c(0.19, 0.32, 0.48, 0.007)
alpha = 0.9
beta = 45

# class probabilities for all genes:
prob_mat = class_probabilities(yeast$Non, theta = c(alpha, beta, p), cvec = cvec)

# class with the highest probability:
prob_max = apply(prob_mat, 1, which.max)

}
\author{
Johanna Bertl
}
