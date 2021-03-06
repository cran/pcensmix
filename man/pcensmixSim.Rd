% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcensmix.R
\name{pcensmixSim}
\alias{pcensmixSim}
\alias{pcensmixSim.pcgen}
\title{Fitting a Normal Mixture Model to a Simulated Progressive Type-II Censored
Data Using EM Algorithm}
\usage{
pcensmixSim(Pdat, ...)

\method{pcensmixSim}{pcgen}(Pdat, r, p, param, iteration = 1e+05,
  INERiter = 20, ...)
}
\arguments{
\item{Pdat}{an object of class \code{"pcgen"} created by function
\code{\link[pcensmix]{pcgen}} or a two-column matrix (or data.frame)
with first column
giving a vector of censored version of a two-component mixed normal data,
and the other one indicating the censoring status associated with them (1
if not censored, otherwise zero).}

\item{...}{additinal arguments to pass by.}

\item{r}{desired number of failures to observe.}

\item{p}{a parameter controlling the amount of censoring. The action of
censoring individuals after each failure occurs with probabilty \code{p}
from binomial distribution at each stage. If \code{p = 0}, there will be no
censoring.}

\item{param}{a numeric vector; used as starting values for the EM and
simulating a new data to replace in case of happening singularity in the
likelihood.}

\item{iteration}{the maximum number of required iteration for the EM
algorithm until convergence-- default value is 1e+05.}

\item{INERiter}{the maximum number of required iteration for the second EM
algorithm-- default is 20.}
}
\value{
\code{pcensmixSim} gives an object of class \code{data.frame}
  containing the following components: \item{muhat1,sigmahat1}{component one
  parameter estimates (\eqn{\hat\mu_1}{\hat{\mu_1}},
  \eqn{\hat\sigma_1}{\hat{\sigma_1}} )} \item{muhat2,sigmahat2}{component two
  parameter estimates (\eqn{\hat\mu_2}{\hat{\mu_2}},
  \eqn{\hat{\sigma_2}}{\hat{\sigma_2} )}} \item{pihat}{estimation of mixture
  proportion \eqn{\hat\pi}{\hat{\pi}}} \item{se.muhat1,se.sigmahat1}{standard
  errors of \eqn{\hat\mu_1}{\hat{\mu_1}} and \eqn{\hat{\sigma_1}}}
  \item{se.muhat2,se.sigmahat2}{standard errors of
  \eqn{\hat\mu_2}{\hat{\mu_2}} and \eqn{\hat\sigma_2}{\hat{\sigma_2}}}
  \item{se.pihat}{standard error of \eqn{\hat\pi}{\hat{\pi}}}
  \item{no.fails.comp1,no.fails.comp2}{number of failures from each mixture
  component} \item{no.cens.comp1,no.cens.comp2}{number of censored
  observations from each mixture component} \item{ll}{log-likelihood value}
  \item{datachange_flag}{\code{TRUE} if data has been replaced by a newly
  generated one}
}
\description{
{}

This function fits a normal mixture model
to progressive Type-II censored mixture data by dealing with the two aspects of
missing data, latent mixture components and the censored data, using a maximum
likelihood estimation through a constrained two-layer EM algorithm.
}
\details{
This function fits a two-component normal mixture model to simulated
  progressive Type-II censored data with density function \deqn{\pi
  (\frac{1}{\sigma_1})\,  \phi\!\! \left[\frac{(z - \mu_1)}{\sigma_1}\right]
  + (1 - \pi) (\frac{1}{\sigma_2})\,  \phi\!\! \left[\frac{(z -
  \mu_2)}{\sigma_2}\right]}{\pi (1/ \sigma_1) \phi[ (z - \mu_1) / \sigma_1]
  + (1 - \pi) (1/ \sigma_2) \phi[ (z - \mu_2) / \sigma_2]} where \eqn{\phi}
  is the standard normal density.

  It uses a constrained two-layer EM algorithm to deal with the two forms of
  missing data: the censored survival times and the mixture component labels.
  Given the EM algorithm is at a particular iteration: (i) first, in the
  E-step it obtains the mixture component indicator estimates given the
  current parameter estimates and the observed data. (ii) Next, for
  re-estimation of the unknown parameters, a new EM algorithm is nested in
  the M-step of the initial EM algorithm to deal with the estimation of the
  missing censored survival times and consequently building the maximum
  likelihood equations. These steps are repeated until the model converges.
}
\note{
\itemize{ \item In fitting the model, to overcome the problem of
  singularity and model non-identifiability that might happen in some cases
  depending on the true means and standard deviations of the components, we
  use the constraint proposed by Hathaway (1985). Based on this, the
  ratios of the standard deviations are considered to be greater than a
  pre-fixed constant. We consider the constraints \eqn{\sigma_1/\sigma_2>0.1}
  and \eqn{\sigma_2/\sigma_1>0.1} which lead to obtain a parameter space with
  no singularities. In case this conditions doesn't hold, the data will be
  replaced by a new simulated one and \code{datachange_flag} will appear as
  \code{TRUE} in the output. \item See \code{\link[pcensmix]{pcgen}} for the
  definition of censored version of data.}
}
\examples{
\dontrun{
set.seed(100)

Pdat<- pcgen(r = 60, p = 0.3, data = mixgen(N = 100, dist1 = 'norm',
                 dist2 = 'norm', control = list(12, 2, 14, 5, 0.35)))
pcensmixSim(Pdat, r = 60, p = 0.3, param=c(12, 2, 14, 5, 0.35))}

}
\seealso{
\code{\link[pcensmix]{pcgen}}, \code{\link[pcensmix]{run_pcensmix}}
}
\author{
Lida Fallah, John Hinde

  Maintainer: Lida Fallah <l.fallah22@gmail.com>
}
