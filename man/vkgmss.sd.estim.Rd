\name{vkgmss.sd.estim}
\alias{vkgmss.sd.estim}

\title{Kernel estimation of the standard deviation function}

\description{This function computes the kernel (Nadaraya-Watson) estimation of the standard deviation function.}

\usage{
vkgmss.sd.estim(x, data.X, data.Y, bandwidth,
		kernel.function = kernel.function.epan)
}

\arguments{
  \item{x}{a numeric vector.}

  \item{data.X}{a numeric data vector used to obtain the kernel estimator of the standard deviation function.}

  \item{data.Y}{a numeric data vector used to obtain the kernel estimator of the standard deviation function.}

  \item{bandwidth}{bandwidth used to obtain the kernel estimator of the standard deviation function.}

  \item{kernel.function}{kernel function used to obtain the kernel estimator of the standard deviation function. Default option is "kernel.function.epan" which corresponds to the Epanechnikov kernel function.}
}

\details{Inappropriate bandwidth or x choices can produce "NaN" values in function estimates.}

\references{

I. Van Keilegom, W. Gonzalez Manteiga, and C. Sanchez Sellero. Goodness-of-fit tests in
parametric regression based on the estimation of the error distribution. \emph{Test}, 17, 401:415,
2008.

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

# Estimation of residuals standard deviation

bandwidth = 0.75 # Here, the bandwidth is arbitrarily fixed

xgrid = seq(0,5,by=0.1)
sd = vkgmss.sd.estim(xgrid,data.X,data.Y,bandwidth)

plot(xgrid,sd , type='l',xlab='X',ylab='SD(X)')
abline(h=0.3)

}
