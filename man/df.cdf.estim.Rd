\name{df.cdf.estim}
\alias{df.cdf.estim}

\title{Local linear estimation of the conditional distribution function}

\description{This function computes the local linear estimation of the conditional distribution function.}

\usage{df.cdf.estim(x, y, data.X, data.Y, bandwidth, kernel.function = kernel.function.epan)}

\arguments{
  \item{x}{a numeric vector.}

  \item{y}{a numeric vector.}

  \item{data.X}{a numeric data vector used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{data.Y}{a numeric data vector used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{bandwidth}{bandwidth used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{kernel.function}{kernel function used to obtain the nonparametric estimator of the conditional distribution function. Default option is "kernel.function.epan" which corresponds to the Epanechnikov kernel function.}
}

\details{Inappropriate bandwidth, x or y choices can produce "NaN" values in cumulative distribution function estimates.}

\references{

G. R. Ducharme and S. Ferrigno. An omnibus test of goodness-of-fit for conditional distributions
with applications to regression models. \emph{Journal of Statistical Planning and
Inference}, 142, 2748:2761, 2012.

R. Azais, S. Ferrigno and M-J Martinez. cvmgof: An R package for Cramer-von Mises goodness-of-fit tests in regression models. Submitted. January 2021.\href{https://hal.archives-ouvertes.fr/hal-03101612}{hal-03101612}





}

\author{Romain Azais, Sandie Ferrigno and Marie-Jose Martinez}

\examples{
set.seed(1)

require(lattice) # Only for plotting conditional CDF

# Data simulation
n = 25 # Dataset size
data.X = runif(n,min=0,max=5) # X
data.Y = 0.2*data.X^2-data.X+2+rnorm(n,mean=0,sd=0.3) # Y

########################################################################

# Estimation of the link function

bandwidth = 0.75 # Here, the bandwidth is arbitrarily fixed

# Estimation of the conditional CDF
xgrid = seq(0.5,4.5,by=0.1)
ygrid = seq(-1,3,by=0.1)
cdf_df = df.cdf.estim(xgrid,ygrid,data.X,data.Y,bandwidth)

wireframe(cdf_df, drape=TRUE,
          col.regions=rainbow(100),zlab='CDF(y|x)',xlab='x',ylab='y',zlim=c(0,1.01))
}





