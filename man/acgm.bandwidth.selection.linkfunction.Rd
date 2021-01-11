\name{acgm.bandwidth.selection.linkfunction}
\alias{acgm.bandwidth.selection.linkfunction}

\title{Bandwidth selection of the link function under the null hypothesis}

\description{This function computes the optimal bandwidth of the link function under the null hypothesis}

\usage{acgm.bandwidth.selection.linkfunction(data.X.H0, data.Y.H0, linkfunction.H0,
		kernel.function = kernel.function.epan, verbose = TRUE)}

\arguments{
  \item{data.X.H0}{a numeric data vector used to obtain the nonparametric estimator of the regression function under the null hypothesis.}

  \item{data.Y.H0}{a numeric data vector used to obtain the nonparametric estimator of the regression function under the null hypothesis.}

  \item{linkfunction.H0}{regression function under the null hypothesis}

  \item{kernel.function}{kernel function used to obtain the nonparametric estimator of the regression function. Default option is "kernel.function.epan" which corresponds to the Epanechnikov kernel function.}

  \item{verbose}{If \code{TRUE}, the R function plots the link function (regression function) under the null hypothesis and the local linear link function estimation on a same graph. Default option is \code{TRUE}.}
}

\references{
J. T. Alcala, J. A. Cristobal, and W. Gonzalez Manteiga. Goodness-of-fit test for linear
models based on local polynomials. \emph{Statistics & Probability Letters}, 42(1), 39:46, 1999.

R. Azais, S. Ferrigno and M-J Martinez. cvmgof: An R package for Cramer-von Mises goodness-of-fit tests in regression models. Submitted. January 2021.\href{https://hal.archives-ouvertes.fr/hal-03101612}{hal-03101612}



}

\author{Romain Azais, Sandie Ferrigno and Marie-Jose Martinez}

\examples{
# Uncomment the following code block
#
# set.seed(1)
#
# # Data simulation
# n = 25 # Dataset size
# data.X = runif(n,min=0,max=5) # X
# data.Y = 0.2*data.X^2-data.X+2+rnorm(n,mean=0,sd=0.3) # Y
#
# ########################################################################
#
# # Bandwidth selection under H0
#
# # We generate a dataset under H0 to estimate the optimal bandwidth under H0
#
# linkfunction.H0 = function(x){0.2*x^2-x+2}
#
# data.X.H0 = runif(n,min=0,max=5)
# data.Y.H0 = linkfunction.H0(data.X.H0)+rnorm(n,mean=0,sd=0.3)
#
# h.opt.acgm = acgm.bandwidth.selection.linkfunction(data.X.H0, data.Y.H0,linkfunction.H0)
}
