\name{vkgmss.test.bootstrap}
\alias{vkgmss.test.bootstrap}


\title{Local test for the regression function}

\description{A local test for the regression function.}

\usage{
vkgmss.test.bootstrap(data.X, data.Y, linkfunction.H0, risk,
    bandwidth = "optimal", kernel.function = kernel.function.epan,
    bootstrap = c(50, "Mammen"), verbose = TRUE)
}


\arguments{
  \item{data.X}{a numeric data vector used to obtain the nonparametric estimator of the error distribution.}

  \item{data.Y}{a numeric data vector used to obtain the nonparametric estimator of the error distribution.}

  \item{linkfunction.H0}{the regression function under the null hypothesis.}

  \item{risk}{a numeric value specifying the risk of rejecting the null hypothesis. The value (1-\code{risk}) corresponds to the confidence level of the statistical test.}

  \item{bandwidth}{the bandwidth used to obtain the nonparametric estimator of the error distribution. If \code{bandwidth}="optimal", the optimal bandwidth of the regression function under the null hypothesis is computed. Default option is "optimal".}

  \item{kernel.function}{the kernel function used to obtain the nonparametric estimator of the error distribution. Default option is "kernel.function.epan".}

  \item{bootstrap}{a numeric vector of length 2. The first value specifies the number of bootstrap datasets (default is "50"). The second value specifies the distribution used for the wild bootstrap resampling.The default is "Mammen" and the other options are "Rademacher" or "Gaussian".}

  \item{verbose}{If \code{TRUE}, the R function displays the optimal bandwidth value obtained under the null hypothesis. Default option is \code{TRUE}. }
}

\details{
From \code{data.X} and \code{data.Y} datasets, wild bootstrap datasets ("50" by default) are built. From each bootstrap dataset, a bootstrap test statistic is computed. The test statistic under the null hypothesis is compared to the distribution of the bootstrap statistics. The test is rejected if the test statistic under the null hypothesis is greater than the (1-\code{risk})-quantile of the empirical distribution of the bootstrap statistics.

An inappropriate bandwidth choice can produce "NaN" values in test statistics.
}

\value{

\code{vkgmss.test.bootstrap} returns a list containing the following components:

\item{decision}{the statistical decision made on whether to reject the null hypothesis or not.}

\item{bandwidth}{the bandwidth used to build the statistics test.}

\item{pvalue}{the p-value of the test statistics.}

\item{test_statistics}{the test statistics value.}

}


\references{
I. Van Keilegom, W. Gonzalez Manteiga, and C. Sanchez Sellero. Goodness-of-fit tests in
parametric regression based on the estimation of the error distribution. \emph{Test}, 17, 401:415,
2008.

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
# # Test (bootstrap) under H0
#
# # We want to test if the link function is f(x)=0.2*x^2-x+2
# # The answer is yes (see the definition of data.Y above)
# # We generate a dataset under H0 to estimate the optimal bandwidth under H0
#
# linkfunction.H0 = function(x){0.2*x^2-x+2}
#
# test_vkgmss.H0 = vkgmss.test.bootstrap(data.X,data.Y,linkfunction.H0,
#                                        0.05,bandwidth='optimal',
#                                        bootstrap=c(50,'Mammen'))
#
#
# ########################################################################
#
# # Test (bootstrap) under H1
#
# # We want to test if the link function is f(x)=0.5*cos(x)+1
# # The answer is no (see the definition of data.Y above)
#
# linkfunction.H1=function(x){0.8*cos(x)+1}
#
# test_vkgmss.H1 = vkgmss.test.bootstrap(data.X,data.Y,linkfunction.H1,
#                                        0.05,bandwidth='optimal',
#                                        bootstrap=c(50,'Mammen'))
}

