\name{df.linkfunction.estim}
\alias{df.linkfunction.estim}

\title{Local linear estimation of the regression function}

\description{This function computes the local linear estimation of the regression function using the local linear estimation of the conditional distribution function.}

\usage{df.linkfunction.estim(x, data.X, data.Y, bandwidth,
		kernel.function = kernel.function.epan)}

\arguments{
  \item{x}{a numeric vector.}

  \item{data.X}{a numeric data vector used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{data.Y}{a numeric data vector used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{bandwidth}{bandwidth used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{kernel.function}{kernel function used to obtain the nonparametric estimator of the conditional distribution function. Default option is "kernel.function.epan" which corresponds to the Epanechnikov kernel function.}
}

\details{Inappropriate bandwidth or x choices can produce "NaN" values in link function estimates.}

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

# Estimation of the link function

bandwidth = 0.75 # Here, the bandwidth is arbitrarily fixed

xgrid = seq(0,5,by=0.1)
ygrid_df = df.linkfunction.estim(xgrid,data.X,data.Y,bandwidth)

plot(xgrid,ygrid_df,type='l',lty=1,lwd=2,xlab='X',ylab='Y',ylim=c(0.25,2.5))
lines(xgrid,0.2*xgrid^2-xgrid+2,lwd=0.5,col='gray')

}


