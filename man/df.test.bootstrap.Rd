\name{df.test.bootstrap}
\alias{df.test.bootstrap}

\title{Global test for the conditional distribution function}

\description{A global test for the conditional distribution function.}

\usage{
df.test.bootstrap(data.X, data.Y, cdf.H0, risk, bandwidth,
    kernel.function = kernel.function.epan,
    bootstrap = c(50, "Mammen"),
    integration.step = 0.01)
}

\arguments{
  \item{data.X}{a numeric data vector used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{data.Y}{a numeric data vector used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{cdf.H0}{the conditional distribution function under the null hypothesis.}

  \item{risk}{a numeric value specifying the risk of rejecting the null hypothesis. The value (1-\code{risk}) corresponds to the confidence level of the statistical test.}

  \item{bandwidth}{the bandwidth used to obtain the nonparametric estimator of the conditional distribution function.}

  \item{kernel.function}{the kernel function used to obtain the nonparametric estimator of the conditional distribution function. Default option is "kernel.function.epan".}

  \item{bootstrap}{a numeric vector of length 2. The first value specifies the number of bootstrap datasets (default is "50"). The second value specifies the distribution used for the wild bootstrap resampling (default is "Mammen").}

  \item{integration.step}{a numeric value specifying integration step. Default is \code{integration.step} = 0.01.}
}

\details{
From \code{data.X} and \code{data.Y} datasets, wild bootstrap datasets ("50" by default) are built. From each bootstrap dataset, a bootstrap test statistic is computed. The test statistic under the null hypothesis is compared to the distribution of the bootstrap statistics. The test is rejected if the test statistic under the null hypothesis is greater than the (1-\code{risk})-quantile of the empirical distribution of the bootstrap statistics.

An inappropriate bandwidth choice can produce "NaN" values in test statistics.
}


\value{

\code{df.test.bootstrap} returns a list containing the following components:

\item{decision}{the statistical decision made on whether to reject the null hypothesis or not.}

\item{bandwidth}{the bandwidth used to build the statistics test.}

\item{pvalue}{the p-value of the test statistics.}

\item{test_statistics}{the test statistics value.}

}


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

# Test (bootstrap) under H0

# Remainder:
# Ducharme and Ferrigno test is on the conditional CDF and not on the link function
# Thus we need to define the conditional CDF associated
# with the link function under H0 to evaluate this test

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

# Test (bootstrap) under H0

test_df.H0 = df.test.bootstrap(data.X,data.Y,cond_cdf.H0,
  0.05,h.opt.df,bootstrap=c(50,'Mammen'),integration.step = 0.01)

########################################################################

# Test (bootstrap) under H1

# We want to test if the link function is f(x)=0.5*cos(x)+1
# The answer is no (see the definition of data.Y above)

linkfunction.H1=function(x){0.8*cos(x)+1}

data.X.H1 = data.X.H0
data.Y.H1 = linkfunction.H1(data.X.H1)+rnorm(n,mean=0,sd=0.3)
h.opt.df = df.bandwidth.selection.linkfunction(data.X.H1 , data.Y.H1,linkfunction.H1)

cond_cdf.H1=function(x,y)
{
  out=matrix(0,nrow=length(x),ncol=length(y))
  for (i in 1:length(x)){
    x0=x[i]
    out[i,]=pnorm(y-linkfunction.H1(x0),0,0.3)
  }
  out
}

test_df.H1 = df.test.bootstrap(data.X,data.Y,cond_cdf.H1,
  0.05,h.opt.df,bootstrap=c(50,'Mammen'),integration.step = 0.01)
}
