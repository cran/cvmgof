\name{vkgmss.linkfunction.estim}
\alias{vkgmss.linkfunction.estim}


\title{Kernel estimation of the regression function}

\description{This function computes the kernel (Nadaraya-Watson) estimation of the regression function.}

\usage{
vkgmss.linkfunction.estim(x, data.X, data.Y, bandwidth,
		kernel.function = kernel.function.epan)
}

\arguments{
  \item{x}{a numeric vector.}

  \item{data.X}{a numeric data vector used to obtain the kernel estimator of the regression function.}

  \item{data.Y}{a numeric data vector used to obtain the kernel estimator of the regression function.}

  \item{bandwidth}{bandwidth used to obtain the kernel estimator of the regression function.}

  \item{kernel.function}{kernel function used to obtain the kernel estimator of the regression function. Default option is "kernel.function.epan" which corresponds to the Epanechnikov kernel function.}
}

\details{Inappropriate bandwidth or x choices can produce "NaN" values in link function estimates.}

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

# Estimation of the link function

bandwidth = 0.75 # Here, the bandwidth is arbitrarily fixed

xgrid = seq(0,5,by=0.1)
ygrid_vkgmss = vkgmss.linkfunction.estim(xgrid,data.X,data.Y,bandwidth)

plot(xgrid,ygrid_vkgmss,type='l',lty=1,lwd=2,xlab='X',ylab='Y',ylim=c(0.25,2.5))
lines(xgrid,0.2*xgrid^2-xgrid+2,lwd=0.5,col='gray')

}

