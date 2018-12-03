\name{df.statistics}
\alias{df.statistics}

\title{Global test statistic for the conditional distribution function}

\description{This function computes the global test statistic for the conditional distribution function.}

\usage{df.statistics(data.X, data.Y, cdf.H0, bandwidth,
    kernel.function = kernel.function.epan, integration.step = 0.01)}

\arguments{
  \item{data.X}{a numeric data vector used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{data.Y}{a numeric data vector used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{cdf.H0}{the conditional distribution function under the null hypothesis.}

  \item{bandwidth}{bandwidth used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{kernel.function}{kernel function used to obtain the nonparametric estimator of the conditional distribution function. Default option is "kernel.function.epan".}

  \item{integration.step}{a numeric value specifying integration step. Default is \code{integration.step} = 0.01.}
}

\details{An inappropriate bandwidth choice can produce "NaN" values in cumulative distribution function estimates.}


\references{
G. R. Ducharme and S. Ferrigno. An omnibus test of goodness-of-fit for conditional distributions
with applications to regression models. \emph{Journal of Statistical Planning and
Inference}, 142, 2748:2761, 2012.

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

h.opt.df = df.bandwidth.selection.linkfunction(data.X.H0 , data.Y.H0,linkfunction.H0)

########################################################################

# Test statistics under H0

cond_cdf.H0 = function(x,y)
{
  out=matrix(0,nrow=length(x),ncol=length(y))
  for (i in 1:length(x)){
    x0=x[i]
    out[i,]=pnorm(y-linkfunction.H0(x0),0,0.3)
  }
  out
}
# cond_cdf.H0 is the conditional CDF associated with linkfunction.H0
# with additive Gaussian noise (standard deviation=0.3)

df.statistics(data.X,data.Y,cond_cdf.H0,h.opt.df)

}




