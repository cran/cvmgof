\name{kernel.function.gauss}
\alias{kernel.function.gauss}

\title{Gaussian kernel}
\description{This function computes the Gaussian kernel used to estimate the conditional distribution function.}

\usage{
kernel.function.gauss(u)
}

\arguments{
  \item{u}{a numeric vector.}
}

\references{
J. Fan and I. Gijbels. Local polynomial modelling and its applications. \emph{Chapman & Hall},
London, 1996.
}

\author{Romain Azais, Sandie Ferrigno and Marie-Jose Martinez}

\examples{
x = runif(10)  #simulating a uniform random sample
kernel.function.gauss(x)
}
