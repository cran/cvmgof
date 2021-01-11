\name{kernel.function.epan}
\alias{kernel.function.epan}

\title{Epanechnikov kernel}
\description{This function computes the Epanechnikov kernel used to estimate the conditional distribution function.}

\usage{
kernel.function.epan(u)
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
kernel.function.epan(x)
}
