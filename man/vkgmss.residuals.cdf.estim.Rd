\name{vkgmss.residuals.cdf.estim}
\alias{vkgmss.residuals.cdf.estim}


\title{Kernel estimation of the error distribution}

\description{This function computes the kernel (Nadaraya-Watson) estimation of the error distribution.}

\usage{vkgmss.residuals.cdf.estim(u, data.X, data.Y, bandwidth,
		kernel.function = kernel.function.epan)}

\arguments{
  \item{u}{a numeric vector.}

  \item{data.X}{a numeric data vector used to obtain the nonparametric estimator of the error distribution.}

  \item{data.Y}{a numeric data vector used to obtain the nonparametric estimator of the error distribution.}

  \item{bandwidth}{bandwidth used to obtain the nonparametric estimator of the error distribution.}

  \item{kernel.function}{kernel function used to obtain the nonparametric estimator of the error distribution. Default option is "kernel.function.epan" which corresponds to the Epanechnikov kernel function.}
}

\details{Inappropriate bandwidth or u choices can produce "NaN" values in error distribution estimates.}

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

# Estimation of residuals cdf

bandwidth = 0.75 # Here, the bandwidth is arbitrarily fixed

egrid = seq(-5,5,by=0.1)
res.cdf = vkgmss.residuals.cdf.estim(egrid,data.X,data.Y,0.5)

plot(egrid,res.cdf , type='l',xlab='e',ylab='CDF(e)')

}

