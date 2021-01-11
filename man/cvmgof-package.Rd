\name{cvmgof}
\alias{cvmgof-package}
\alias{cvmgof}

\docType{package}

\title{cvmgof}

\description{It implements three goodness-of-fit tests to test the validity of the regression function in the regression model. Two of them (Alcala \emph{et al.} '99 and Van Keilegom \emph{et al.} '12) are ``directional'' in that they detect departures from mainly the regression function assumption of the model or ``global'' (Ducharme and Ferrigno '12) with the conditional distribution function. The establishment of such statistical tests requires nonparametric estimators and the use of wild bootstrap methods for the simulations.}

\details{
\tabular{ll}{
Package: \tab cvmgof\cr
Type: \tab Package\cr
Version: \tab 1.0.3\cr
Date: \tab 2021-01-11\cr
License: \tab Cecill\cr
}
}

\author{Romain Azais, Sandie Ferrigno and Marie-Jose Martinez}

\references{

J. T. Alcala, J. A. Cristobal, and W. Gonzalez Manteiga. Goodness-of-fit test for linear
models based on local polynomials. \emph{Statistics & Probability Letters}, 42(1), 39:46, 1999.

G. R. Ducharme and S. Ferrigno. An omnibus test of goodness-of-fit for conditional distributions
with applications to regression models. \emph{Journal of Statistical Planning and
Inference}, 142, 2748:2761, 2012.

I. Van Keilegom, W. Gonzalez Manteiga, and C. Sanchez Sellero. Goodness-of-fit tests in
parametric regression based on the estimation of the error distribution. \emph{Test}, 17, 401:415,
2008.

R. Azais, S. Ferrigno and M-J Martinez. cvmgof: An R package for Cramer-von Mises goodness-of-fit tests in regression models. Submitted. January 2021.\href{https://hal.archives-ouvertes.fr/hal-03101612}{hal-03101612}
}

\seealso{
  \code{\link{acgm.test.bootstrap}}

  \code{\link{df.test.bootstrap}}

  \code{\link{vkgmss.test.bootstrap}}
}

\examples{
# require(lattice) # Only for plotting conditional CDF

########################################################################
# Simulation

set.seed(1)
# The following example tests are computed from only 25 data points
# The seed is fixed to avoid NA in estimates

# Data simulation
n = 25 # Dataset size
data.X = runif(n,min=0,max=5) # X
data.Y = 0.2*data.X^2-data.X+2+rnorm(n,mean=0,sd=0.3) # Y
# plot(data.X,data.Y,xlab='X',ylab='Y',pch='+')

########################################################################
# Estimation of the link function (uncomment the following code block)

# bandwidth = 0.75 # Here, the bandwidth is arbitrarily fixed
#
# xgrid = seq(0,5,by=0.1)
# ygrid_df = df.linkfunction.estim(xgrid,data.X,data.Y,bandwidth)
# ygrid_acgm = acgm.linkfunction.estim(xgrid,data.X,data.Y,bandwidth)
# ygrid_vkgmss = vkgmss.linkfunction.estim(xgrid,data.X,data.Y,bandwidth)
#
# plot(xgrid,ygrid_df,type='l',col='blue',lty=1,lwd=2,xlab='X',ylab='Y',ylim=c(0.25,2.5))
# lines(xgrid,ygrid_acgm,type='l',col='red',lty=2,lwd=2)
# lines(xgrid,ygrid_vkgmss,type='l',col='dark green',lty=3,lwd=2)
# lines(xgrid,0.2*xgrid^2-xgrid+2,lwd=0.5,col='gray')
# # Ducharme and Ferrigno: blue
# # Alcala et al.: red
# # Van Keilegom et al.: dark green
# # true link function: gray
#
# # Estimation of the conditional CDF (only Ducharme and Ferrigno estimator)
# xgrid = seq(0.5,4.5,by=0.1)
# ygrid = seq(-1,3,by=0.1)
# cdf_df = df.cdf.estim(xgrid,ygrid,data.X,data.Y,bandwidth)
#
# wireframe(cdf_df, drape=TRUE,
#           col.regions=rainbow(100),
#           zlab='CDF(y|x)',
#           xlab='x',ylab='y',zlim=c(0,1.01))
#
# # Estimation of residuals cdf (only Van Keilegom et al. estimator)
#
# egrid = seq(-5,5,by=0.1)
# res.cdf_vkgmss = vkgmss.residuals.cdf.estim(egrid,data.X,data.Y,0.5)
#
# plot(egrid,res.cdf_vkgmss,type='l',xlab='e',ylab='CDF(e)')
#
# # Estimation of residuals standard deviation (only Van Keilegom et al. estimator)
#
# sd_vkgmss = vkgmss.sd.estim(xgrid,data.X,data.Y,bandwidth)
#
# plot(xgrid,sd_vkgmss, type='l',xlab='X',ylab='SD(X)')
# abline(h=0.3)

########################################################################
# Bandwidth selection under H0 (uncomment the following code block)

# # We want to test if the link function is f(x)=0.2*x^2-x+2
# # The answer is yes (see the definition of data.Y above)
# # We generate a dataset under H0 to estimate the optimal bandwidth under H0
#
# linkfunction.H0 = function(x){0.2*x^2-x+2}
#
# data.X.H0 = runif(n,min=0,max=5)
# data.Y.H0 = linkfunction.H0(data.X.H0)+rnorm(n,mean=0,sd=0.3)
#
# h.opt.df = df.bandwidth.selection.linkfunction(data.X.H0, data.Y.H0,linkfunction.H0)
# h.opt.acgm = acgm.bandwidth.selection.linkfunction(data.X.H0, data.Y.H0,linkfunction.H0)
# h.opt.vkgmss = vkgmss.bandwidth.selection.linkfunction(data.X.H0, data.Y.H0,linkfunction.H0)
# # Ducharme and Ferrigno: 1.184604
# # Alcala et al.: 0.7716453
# # Van Keilegom et al.: 0.6780543

########################################################################
# Test statistics under H0 (uncomment the following code block)

# # Remainder:
# # Ducharme and Ferrigno test is on the conditional CDF and not on the link function
# # Thus we need to define the conditional CDF associated
# # with the link function under H0 to evaluate this test
# # Alcala et al. and Van Keilegom et al. tests are on the link function
#
# # Optimal bandwidths estimated at the previous step
# h.opt.df = 1.184604
# h.opt.acgm = 0.7716453
# h.opt.vkgmss = 0.6780543
#
# cond_cdf.H0 = function(x,y)
# {
#   out=matrix(0,nrow=length(x),ncol=length(y))
#   for (i in 1:length(x)){
#     x0=x[i]
#     out[i,]=pnorm(y-linkfunction.H0(x0),0,0.3)
#   }
#   out
# }
# # cond_cdf.H0 is the conditional CDF associated with linkfunction.H0
# # with additive Gaussian noise (standard deviation=0.3)
#
# df.statistics(data.X,data.Y,cond_cdf.H0,h.opt.df)
# acgm.statistics(data.X,data.Y,linkfunction.H0,h.opt.acgm)
# vkgmss.statistics(data.X,data.Y,linkfunction.H0,h.opt.vkgmss)

########################################################################
# Test (bootstrap) under H0

h.opt.df = 1.184604 # Optimal bandwidth estimated above

linkfunction.H0 = function(x){0.2*x^2-x+2}
cond_cdf.H0 = function(x,y)
{
  out=matrix(0,nrow=length(x),ncol=length(y))
  for (i in 1:length(x)){
    x0=x[i]
    out[i,]=pnorm(y-linkfunction.H0(x0),0,0.3)
  }
  out
}

test_df.H0 = df.test.bootstrap(data.X,data.Y,cond_cdf.H0,
                               0.05,h.opt.df,bootstrap=c(20,'Mammen'),
                               integration.step = 0.1)
test_acgm.H0 = acgm.test.bootstrap(data.X,data.Y,linkfunction.H0,
                                   0.05,bandwidth='optimal',bootstrap=c(20,'Mammen'),
                                   verbose=FALSE,integration.step = 0.1)
test_vkgmss.H0 = vkgmss.test.bootstrap(data.X,data.Y,linkfunction.H0,
                                       0.05,bandwidth='optimal',bootstrap=c(20,'Mammen'),
                                       verbose=FALSE)

# test_acgm$decision is a string: 'accept H0' or 'reject H0'
# test_acgm$bandwidth is a float (optimal bandwidth under H0
# (only for Alcala and Van Keilegom tests) if bandwidth = 'optimal')
# test_acgm$pvalue is a float but it could be a string
# ('< 0.02' for instance or 'None' if the test can not be evaluated)
# test_acgm$test_statistics is a float but it could be a string
# ('None' if the test can not be evaluated)

# The 3 tests accept H0

########################################################################
# Test (bootstrap) under H1 (uncomment the following code block)

# # We want to test if the link function is f(x)=0.5*cos(x)+1
# # The answer is no (see the definition of data.Y above)
#
# linkfunction.H1=function(x){0.8*cos(x)+1}
#
# plot(xgrid,linkfunction.H0(xgrid),type='l',ylim=c(-1,2))
# lines(xgrid,linkfunction.H1(xgrid),type='l')
#
# data.X.H1 = data.X.H0
# data.Y.H1 = linkfunction.H1(data.X.H1)+rnorm(n,mean=0,sd=0.3)
# h.opt.df = df.bandwidth.selection.linkfunction(data.X.H1, data.Y.H1,linkfunction.H1)
#
# cond_cdf.H1=function(x,y)
# {
#   out=matrix(0,nrow=length(x),ncol=length(y))
#   for (i in 1:length(x)){
#     x0=x[i]
#     out[i,]=pnorm(y-linkfunction.H1(x0),0,0.3)
#   }
#   out
# }
#
# test_df.H1 = df.test.bootstrap(data.X,data.Y,cond_cdf.H1,
#                                0.05,h.opt.df,bootstrap=c(20,'Mammen'),
#                                integration.step = 0.1)
# test_acgm.H1 = acgm.test.bootstrap(data.X,data.Y,linkfunction.H1,
#                                    0.05,bandwidth='optimal',bootstrap=c(20,'Mammen'),
#                                    integration.step = 0.1,verbose=FALSE)
# test_vkgmss.H1 = vkgmss.test.bootstrap(data.X,data.Y,linkfunction.H1,
#                                        0.05,bandwidth='optimal',bootstrap=c(20,'Mammen'),
#                                        verbose=FALSE)
#
# # From only 25 points, only Van Keilegom et al. test rejects H0
# # while Ducharme and Ferrigno and Alcala et al. tests accept H0

########################################################################
}
