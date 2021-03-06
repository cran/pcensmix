% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcensmix.R
\name{run_pcensmix}
\alias{run_pcensmix}
\title{Generating Progressively Type-II Censored Mixture Data and Fitting a Model}
\usage{
run_pcensmix(N, r, p, param, repetition = 100)
}
\arguments{
\item{N}{population size.}

\item{r}{required number of failures to observe.}

\item{p}{a parameter controlling the amount of censoring. The action of
censoring individuals after each failure occurs with probabilty \code{p}
from a binomial distribution at each stage. If \code{p = 0}, there will be
no censoring.}

\item{param}{a numeric vector; used as starting values for the EM and
simulating a new data to replace in case of happening singularity in the
likelihood.}

\item{repetition}{the required number of repetition of the algorithm-- default
is 100.}
}
\value{
It returns the parameter estimates given by
 \code{\link[pcensmix]{pcensmixSim}} with the desired number of repetitions.
 In each repetition it generates a new normal mixture progressive Type-II
 censored dataset from the same true parameter values and fits a model.
}
\description{
This function implements an algorithm using the
\code{\link[pcensmix]{mixgen}}, \code{\link[pcensmix]{pcgen}} and
\code{\link[pcensmix]{pcensmixSim}} functions to generate data and fit a model
using EM algorithm with a specified number of iterations.
}
\examples{
\dontrun{

## Example 1: with very well separated mixture components
set.seed(3)
f1 <- run_pcensmix(N = 160, r = 120, p = 0.3, param = c(10, 2, 25, 4, 0.3), repetition = 100)
colMeans(f1)

## Example 2.
set.seed(150)
f2 <- run_pcensmix(N = 160, r = 130, p = 0.35, param = c(10, 2, 17, 4, 0.3), repetition = 100)
colMeans(f2)

## Example 3.
set.seed(20)
f3 <- run_pcensmix(N = 160, r = 130, p = 0.3, param = c(20, 6, 22, 12, 0.6), repetition = 100)
colMeans(f3)

}

}
\seealso{
\code{\link[pcensmix]{pcgen}}, \code{\link[pcensmix]{pcensmixSim}}, \code{\link[pcensmix]{mixgen}}
}
\author{
Lida Fallah, John Hinde

Maintainer: Lida Fallah <l.fallah22@gmail.com>
}
