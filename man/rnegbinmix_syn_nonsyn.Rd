% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/rnegbinmix_syn_nonsyn.R
\name{rnegbinmix_syn_nonsyn}
\alias{rnegbinmix_syn_nonsyn}
\title{Simulating synonymous and non-synonymous mutation counts}
\usage{
rnegbinmix_syn_nonsyn(n, alpha, beta, c, p)
}
\arguments{
\item{n}{}

\item{alpha,}{scalar, positive}

\item{beta}{positive}

\item{c}{numeric, positive, length k}

\item{p}{numeric, probability vector, length k}
}
\description{
rnegbinmix_syn_nonsyn simulates counts of synonymous and non-synonymous mutations, using a negative binomial and a negative binomial mixture, respectively, with the same parameters alpha and beta.
}
\details{
The function \link{MASS::rnegbin} from package MASS is used for the simulation of negative binomial data. The function \link{rnegbinmix} is used for the negative binomial mixture.

Currently, only a mixture of two components is implemented for the non-synonymous mutations.
}
\author{
Johanna Bertl
}
