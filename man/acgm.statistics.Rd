\name{acgm.statistics}
\alias{acgm.statistics}

\title{Local test statistic for the regression function}

\description{This function computes the local test statistic for the regression function.}

\usage{acgm.statistics(data.X, data.Y, linkfunction.H0,
		bandwidth = "optimal",
		kernel.function = kernel.function.epan,
		integration.step = 0.01,
		verbose = TRUE)}

\arguments{
  \item{data.X}{a numeric data vector used to obtain the nonparametric estimator of the regression function.}

  \item{data.Y}{a numeric data vector used to obtain the nonparametric estimator of the regression function.}

  \item{linkfunction.H0}{the regression function under the null hypothesis.}

  \item{bandwidth}{bandwidth used to obtain the nonparametric estimator of the regression function. If \code{bandwidth} = "optimal", the optimal bandwidth of the regression function under the null hypothesis is computed. Default option is "optimal".}

  \item{kernel.function}{kernel function used to obtain the nonparametric estimator of theregression function. Default option is "kernel.function.epan".}

  \item{integration.step}{a numeric value specifying integration step. Default is \code{integration.step} = 0.01.}

  \item{verbose}{If \code{TRUE}, the R function displays the optimal bandwidth value obtained under the null hypothesis. Default option is \code{TRUE}.}
}


\references{
J. T. Alcala, J. A. Cristobal, and W. Gonzalez Manteiga. Goodness-of-fit test for linear
models based on local polynomials. \emph{Statistics & Probability Letters}, 42(1), 39:46, 1999.

R. Azais, S. Ferrigno and M-J Martinez. cvmgof: An R package for Cramér-von Mises goodness-of-fit tests in regression
models. 2018. Preprint in progress.

}

\author{Romain Azais, Sandie Ferrigno and Marie-Jose Martinez}

\examples{
set.seed(1)

# Data simulation
n = 25 # Dataset size
data.X = runif(n,min=0,max=5) # X
data.Y = 0.2*data.X^2-data.X+2+rnorm(n,mean=0,sd=0.3) # Y

########################################################################

# Bandwidth selection under H0

# We want to test if the link function is f(x)=0.2*x^2-x+2
# The answer is yes (see the definition of data.Y above)
# We generate a dataset under H0 to estimate the optimal bandwidth under H0

linkfunction.H0 = function(x){0.2*x^2-x+2}

data.X.H0 = runif(n,min=0,max=5)
data.Y.H0 = linkfunction.H0(data.X.H0)+rnorm(n,mean=0,sd=0.3)

h.opt.acgm = acgm.bandwidth.selection.linkfunction(data.X.H0 , data.Y.H0,linkfunction.H0)

########################################################################

# Test statistics under H0

acgm.statistics(data.X,data.Y,linkfunction.H0,h.opt.acgm)

}


