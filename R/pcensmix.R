#' Generate and fit a model to censored mixture data
#'
#' @docType package
#' @description This package provides functions to generate two-component
#'   mixture data from various different mixture distribution, generate
#'   progressive Type-II censored data in a mixture structure and fit a normal
#'   mixture model using a constrained EM algorithm. In addition, it can create a
#'   progressive Type-II censored version of a given real dataset and fit a
#'   normal mixture model to.  Main functions are \link[pcensmix]{pcgen},
#'   \link[pcensmix]{pcensmixSim} and \link[pcensmix]{pcensmixR}. Example datasets are included for
#'   illustration.
#'
#' @details
#' Package: pcensmix
#'
#' Type: Package
#'
#' Version: 1.2-1
#'
#' Date: 2017-07-24
#'
#' License: GPL (>= 2)
#'
#' @author Lida Fallah <l.fallah22@gmail.com> and John Hinde
#'
#' Maintainer: Lida Fallah <l.fallah22@gmail.com>
#'
#'
#' @references Aitkin, M., Francis, B., Hinde, J. and Darnell, R., (2009).
#' Statistical Modelling in R. Oxford: Oxford University Press.
#'
#' Balakrishnan, N. and Aggarwala, R., (2000). Progressive Censoring: Theory,
#' Methods, and Applications. Springer Science & Business Media.
#'
#' Hathaway, R.J., (1985). A constrained formulation of maximum-likelihood
#' estimation for normal mixture distributions. The Annals of Statistics,
#' 795-800.
#'
#' McLachlan, G. and Krishnan, T., (2007). The EM Algorithm and Extensions. John
#' Wiley & Sons.
#'
#' McLachlan, G. and Peel, D., (2004). Finite Mixture Models. John Wiley & Sons.
#'
#' @name pcensmix-package
#' @keywords package
NULL

#' Failure times of insulating fluid.
#'
#' A dataset containing measurements of 19 failure times (in minutes) for an
#' insulating fluid between two electrodes subject to a voltage of 34 KV, see
#' Nelson (1982).
#'
#' @docType data
#' @usage insulate
#'
#' @format A vector of length 19.
#'
#' @references
#'
#' Nelson, W., (1982). Applied Life Data Analysis. Wiley, New York.
"insulate"


#' Blood pressure of mine workers.
#'
#' A dataframe containing the blood pressure of a population of mine workers in
#' Ghana.
#'
#' @docType data
#' @usage blood
#'
#' @format A data frame with 495 rows and 3 variables
#'
#'   \itemize{
#'
#'   \item Systolic.BP. Systolic blood pressure \item Diastolic.BP. Diastolic
#'   blood pressure \item Amplitude. Amplitude
#'
#'   }
#'
#' @references
#'
#' Boehning, D., (2000). Computer-assisted Analysis of Mixtures and
#' Applications: Meta-analysis, Disease mapping and Others. CRC press.
#'
#' Gunga, H.C., Forson, K., Amegby, N. and Kirsch, K., (1991). Lebensbedingungen
#' und gesundheitszustand von berg-und fabrikarbeitern im tropischen regenwald
#' von Ghana. Arbeitsmedizin Sozialmedizin Praventivmediz in, 17-25.
"blood"


#' Generating Mixture Datasets
#'
#' This function generates two-component mixture data from a various different
#' mixture distributions.
#'
#' @param N population size.
#' @param dist1,dist2 respective distributions of the first and second mixture
#'   components to be simulated from. For Normal, Log-normal, Weibull, Gamma,
#'   Cauchy and Beta distributions, they must be provided as \code{'norm'},
#'   \code{'lnorm'}, \code{'weibull'}, \code{'gamma'}, \code{'cauchy'} and
#'   \code{'beta'} respectively. The Exponential distribution can be included as
#'   a Gamma distribution with scale parameter one.
#' @param control a list of parameters for controlling the simulation process.
#'   This includes parameters of the first and second mixture component
#'   distributions and the mixture proportion \eqn{\pi} which should be provided
#'   in the order mentioned, i.e., parameters for the first component come
#'   first, those of the second component come after and then the value of the
#'   mixing proportion. All values should be provided in a list of numerics.
#'   Note parameters for each component distribution should be added in the
#'   order as required in their generator functions, see
#'   \code{\link[stats]{rnorm}}, \code{\link[stats]{rlnorm}},
#'   \code{\link[stats]{rweibull}}, \code{\link[stats]{rgamma}},
#'   \code{\link[stats]{rcauchy}} and \code{\link[stats]{rbeta}}.
#'
#'
#' @return An object of class \code{data.frame} containing the following
#'   information: \item{z}{mixture data} \item{label}{component indicator}
#'
#' @details It generates a two-component mixture dataset from a density function
#'   \deqn{\pi f_1 + (1 - \pi) f_2,}{%
#'   \pi f_1 + (1 - \pi) f_2,} where \eqn{f_1} and \eqn{f_2} are the first
#'   and the second mixture component distributions respectively.
#'
#'
#' @examples
#' ## Generate a sample from a two component Normal-Weibull mixture distribution
#' ## with mixing components as N(12, 2) and Weibull(15, 4), mixing proportion 0.3
#' ## and size of N = 20.
#'
#' mixture<- mixgen(N = 20, dist1 = 'norm', dist2 = 'weibull', control = list(12, 2, 15, 4, 0.3))
#'
#' @seealso \code{\link[pcensmix]{pcgen}}
#' @author Lida Fallah, John Hinde
#'
#' Maintainer: Lida Fallah <l.fallah22@gmail.com>
#' @export
#' @importFrom stats dnorm pnorm rbinom rgamma rlnorm rnorm runif rweibull
#'   rcauchy rbeta var
#' @importFrom utils tail


mixgen<-function(N, dist1 , dist2 , control){

  z<- rep(-1, N)
  label<- rep(-1,N)
  d<-runif(N)

  pi<-tail(control, 1)


  for(i in 1:N){
    if (d[i]<pi){
      z[i]<- eval(call(paste('r', dist1, sep=''), 1, control[[1]], control[[2]]))
      label[i]<-1
    }
    else{
      z[i]<- eval(call(paste('r', dist1, sep=''), 1, control[[3]], control[[4]]))
      label[i]<-0
    }
  }

  data.frame(z, label)
}



#' Creating a Progressively Type-II Censored Version of a Given Dataset
#'
#' This function implements an algorithm for generating a progressive Type-II
#' censored version of a specified dataset.
#'
#' @param r desired number of failures to observe.
#' @param p a parameter controlling the amount of censoring. The action of
#'   censoring individuals after each failure occurs with probabilty \code{p}
#'   from binomial distribution at each stage. If \code{p = 0}, no censoring
#'   will happen.
#' @param data a numeric vector of a real dataset (mixture/not mixture) or an
#'   object of class \code{data.frame} generated by
#'   \code{\link[pcensmix]{mixgen}}.
#'
#' @return An object of class \code{"pcgen"} containing the following
#'   information: \item{original_data}{original mixture data}
#'   \item{label}{component membership indicator for the original simulated
#'   mixture data. This will not be returned if a real data has been used.}
#'   \item{censored_version_of_data}{progressive Type-II censored version of
#'   data, i.e., each observation is equal to the actual observed survival time
#'   in the event of failure and is equal to the latest observe failure time if
#'   it is associated to an unobserved censored observation. Notice that they
#'   are order statistics.} \item{component_indicator}{component indicator
#'   associated with the censored_verison_of_data. This will not be returned if
#'   a real data has been used.} \item{censoring_indicator}{censoring indicator
#'   associated with the censored_verison_of_data.}
#'
#' @details It creates a progressive Type-II censored version of a given real
#'   dataset or a simulated dataset from \code{\link[pcensmix]{mixgen}}. The
#'   output of this function can be passed as an argument to
#'   \code{\link[pcensmix]{pcensmixR}} or \code{\link[pcensmix]{pcensmixSim}} for the purpose of fitting a
#'   normal mixture model to the progressively censored dataset.
#' @note See \code{\link[pcensmix]{print.pcgen}} for printing data of class
#' \code{"pcgen"}.
#' @examples
#' ## 1. Generate a progressive Type-II censored data from a simulated mixture data with
#' ## allowing for censoring with controlling parameters p = 0.3 and r = 12.
#' set.seed(0)
#' mixture <- mixgen(N = 20, dist1 = 'norm', dist2 = 'weibull', control = list(12, 2, 15, 4, 0.3))
#' Pdat0 <- pcgen(r = 12, p = 0.3, data = mixture)
#' print(Pdat0)
#'
#'
#' ## 2. Examples of generating a progresively Type-II censored data
#'
#' set.seed(0)
#' Pdat1 <- pcgen(r = 6, p = 0.3, data = insulate)
#' print(Pdat1)
#'
#' set.seed(100)
#' Pdat2 <- pcgen(r = 260, p = 0.35, data = blood$Systolic.BP)
#' print(Pdat2)
#'
#'
#' @seealso \code{\link[pcensmix]{mixgen}}, \code{\link[pcensmix]{print.pcgen}}.
#'
#' @author Lida Fallah, John Hinde
#'
#' Maintainer: Lida Fallah <l.fallah22@gmail.com>
#' @export
#' @importFrom stats dnorm pnorm rbinom rgamma rlnorm rnorm runif rweibull var


pcgen<-function(r, p, data){

  dat <- data

  if(is.data.frame(dat)){
    N <- nrow(dat)
    z <- dat[,1] ; label <- dat[,2]}

  if(is.vector(dat)){
    N <- length(dat)
    z <- dat ; label <- rep(NA, length(z))
  }

  sort.z = sort(z)

  if(r!=N){
    lastremov <- -1
    while(lastremov < 0){
      Rstar <- rbinom((r-1),1,p)
      lastremov <- N-r-sum(Rstar)
    }
  }else
  {Rstar = rep(0, (r-1)); lastremov = 0}

  R <- c(Rstar, lastremov)

  times <- sort.z
  Cstar = label[order(z)]
  W = Z = NULL

  C <- NULL
  Z <- c(Z, times[1])
  W <- c(W, 1)

  index <- 1:length(times)

  for(i in 1:length(Rstar)){

    times <- times[!index %in% index[1]]

    index <- index[-1] ;

    C = c(C, Cstar[1])
    Cstar <- Cstar[-1]


    samp <- sample(index, R[i], replace=FALSE)
    times <- times[!index %in% samp]


    C = c(C, Cstar[index %in% samp]);  Cstar <- Cstar[!index %in% samp]
    index <- index[!index %in% samp]

    W = c(W, rep(0,length(samp)), 1)
    Z <- c( Z, rep(Z[length(Z)], length(samp)), times[1] )
  }

  times <- times[!index %in% index[1]]
  C <- c( C, Cstar)
  W <- c(W, rep(0, length(times)))
  Z <- c(Z, rep(Z[length(Z)], length(times)))

  Pdat <- data.frame('original_data'=z, 'label'=label, 'censored_version_of_data'=Z,
                     'component_indicator'=C, 'censoring_indicator'=W)

  if(is.vector(dat)){
    Pdat <- Pdat[, -c(2,4)]
  }

  class(Pdat) <-  "pcgen"

  invisible(Pdat)
}

#' Print Method for pcgen Objects
#'
#' This function prints the progressive censored data generated by the S3
#' class \code{\link{pcgen}}.
#'
#' @param x object of class \code{pcgen}.
#' @param ... optional arguments to pass by.
#'
#' @return This function uses the generic function \code{\link[base]{print}}
#' to print the dataset of class \code{"pcgen"} in a nice format.
#' @examples
#' ## Generate a two component normal mixture data,
#' Pdat <- pcgen(r = 80, p = 0.3, data = mixgen(N = 100, dist1 = 'norm',
#'               dist2 = 'norm', control = list(12, 2, 14, 4, 0.3)))
#' # and print it.
#' print(Pdat)
#'
#' @export

print.pcgen <- function(x, ...) {

  print(data.frame("original_data"= x$original_data, "censored_version_of_data" = x$censored_version_of_data, "censoring_indicator" = x$censoring_indicator), ...)
#
#   invisible(x)
}



#' Internal function for pcensmix
#'
#' .f_u function is for internal usage only and is not intended to be called
#' by the user.
#'
#' @keywords internal
#'
#' @author Lida Fallah, John Hinde
#'
#' Maintainer: Lida Fallah <l.fallah22@gmail.com>
#' @name pcensmix-internal
#' @aliases pcensmix-internal

.f_u<-function(eps, Pdat, muhat1, muhat2, sigmahat1, sigmahat2, pi_hat, muhat1_est, muhat2_est, sigmahat1_est, sigmahat2_est, pi_hat_est, ll, W, Z, INERiter, N){


  ExpU<-W*( pi_hat*dnorm(Z, muhat1 , sigmahat1)/(pi_hat*dnorm(Z, muhat1,sigmahat1)+(1-pi_hat)*dnorm(Z, muhat2, sigmahat2)) ) + (1-W)*( pi_hat*(1-pnorm(Z, muhat1 , sigmahat1))/(pi_hat*(1-pnorm(Z, muhat1,sigmahat1))+(1-pi_hat)*(1-pnorm(Z, muhat2, sigmahat2))) )


  if((sigmahat1/sigmahat2)< 0.1 & (sigmahat2/sigmahat1)< 0.1){
    ExpU<- ((1-eps))*ExpU + eps
  }


  for (k in 1:INERiter){

    x_tilde<-rep(0,length(W))
    xSq_tilde<-rep(0,length(W))
    y_tilde<-rep(0,length(W))
    ySq_tilde<-rep(0,length(W))

    for(j in 1:length(W)){
      if(W[j]==0 & (Z[j]-muhat1) < 6*sigmahat1){
        x_tilde[j]<- W[j]*Z[j] + (1-W[j])*(muhat1+sigmahat1*dnorm((Z[j]-muhat1)/sigmahat1, 0, 1)/(1-pnorm((Z[j]-muhat1)/sigmahat1, 0, 1))  )

      }else{x_tilde[j]<- W[j]*Z[j]}  }

    for(j in 1:length(W)){
      if(W[j]==0 & (Z[j]-muhat1) < 6*sigmahat1){
        xSq_tilde[j]<- W[j]*Z[j]^2 + (1-W[j])*( muhat1^2 + sigmahat1^2+ sigmahat1*(muhat1+Z[j])*dnorm((Z[j]-muhat1)/sigmahat1, 0, 1)/(1-pnorm((Z[j]-muhat1)/sigmahat1, 0, 1)) )
      }else{ xSq_tilde[j]<- W[j]*Z[j]^2}   }


    for(j in 1:length(W)){
      if(W[j]==0 & (Z[j]-muhat2) < 6*sigmahat2){
        y_tilde[j]<- W[j]*Z[j] + (1-W[j])*(muhat2+sigmahat2*dnorm((Z[j]-muhat2)/sigmahat2, 0, 1)/(1-pnorm((Z[j]-muhat2)/sigmahat2, 0, 1))  )
      }else{y_tilde[j]<- W[j]*Z[j] }     }

    for(j in 1:length(W)){
      if(W[j]==0 & (Z[j]-muhat2) < 6*sigmahat2){
        ySq_tilde[j]<- W[j]*Z[j]^2 + (1-W[j])*( muhat2^2 + sigmahat2^2+ sigmahat2*(muhat2+Z[j])*dnorm((Z[j]-muhat2)/sigmahat2, 0, 1)/(1-pnorm((Z[j]-muhat2)/sigmahat2, 0, 1)) )
      }else{ySq_tilde[j]<- W[j]*Z[j]^2}  }

    muhat1 <- sum(x_tilde*ExpU)/sum(ExpU)

    muhat2<- sum(y_tilde*(1-ExpU))/sum(1-ExpU)

    sigmahat1<- sqrt( sum((xSq_tilde -2*muhat1*x_tilde+ muhat1^2)*ExpU)/sum(ExpU) )

    sigmahat2<- sqrt( sum((ySq_tilde -2*muhat2*y_tilde + muhat2^2)*(1-ExpU))/sum(1-ExpU) )

    pi_hat<- sum(ExpU)/N

  }

  muhat1_est<-c(muhat1_est, muhat1)
  muhat2_est<-c(muhat2_est, muhat2)
  sigmahat1_est<- c(sigmahat1_est, sigmahat1)
  sigmahat2_est<- c(sigmahat2_est, sigmahat2)
  pi_hat_est<-c(pi_hat_est, pi_hat)


  l<- sum(log( (pi_hat*dnorm(Z, muhat1 , sigmahat1)+ (W*ExpU<=1e-12))^(W*ExpU))  + log(((1-pi_hat)*dnorm(Z, muhat2 , sigmahat2)+ (W*(1-ExpU)<=1e-12))^(W*(1-ExpU)))  + log((pi_hat*dnorm(x_tilde, muhat1 , sigmahat1)+ ((1-W)*(ExpU)<=1e-12))^((1-W)*ExpU) ) + log(( (1-pi_hat)*dnorm(y_tilde, muhat2 , sigmahat2)+ ((1-W)*(1-ExpU)<=1e-12) )^((1-W)*(1-ExpU)) ) )

  ll<-c(ll, l)


  list('muhat1'=muhat1,'muhat2'=muhat2,'sigmahat1'= sigmahat1,'sigmahat2'=sigmahat2,'pi_hat'= pi_hat, 'ExpU'=ExpU, 'x_tilde'=x_tilde, 'y_tilde'=y_tilde, 'xSq_tilde'=xSq_tilde, 'ySq_tilde'=ySq_tilde, 'l'=l)
}


#' Fitting a Normal Mixture Model to a Real Progressive Type-II Censored Mixture
#' Data Using EM Algorithm
#'
#' {}
#'
#'
#' @param Pdat an object of class \code{"pcgen"} created by function
#'   \code{\link[pcensmix]{pcgen}} or a two-column matrix (or data.frame)
#'   with first column
#'   giving a vector of censored version of a two-component mixed normal data,
#'   and the other one indicating the censoring status associated with them (1
#'   if not censored, otherwise zero).
#' @param ... additinal arguments to pass by.
#'
#' @export


pcensmixR<-function (Pdat, ...)
{
  UseMethod("pcensmixR")
}




#' Fitting a Normal Mixture Model to a Real Progressive Type-II Censored Mixture
#' Data Using EM Algorithm
#'
#' This function uses a two-layer EM algorithm to fit a mixture model to
#' progressive Type-II censored mixture data by estimating the latent mixture
#' components and the censored data.
#'
#'
#' @param start a numeric vector; used as starting values for the EM algorithm.
#' @param iteration the maximum number of required iteration for the EM
#'   algorithm until convergence-- default value is 1e+05.
#' @param INERiter the maximum number of required iteration for the second EM
#'   algorithm-- default is 20.
#' @param warn logical. shows warning messages if \code{TRUE}, if there is any--
#' default is \code{FALSE}.
#'
#'
#' @return \code{pcensmixR} gives an object of class \code{data.frame}
#'   containing the following components:
#'   \item{muhat1,sigmahat1}{component one parameter
#'   estimates (\eqn{\hat\mu_1}{\hat{\mu_1}}, \eqn{\hat\sigma_1}{\hat{\sigma_1}} )}
#'   \item{muhat2,sigmahat2}{component two parameter
#'   estimates (\eqn{\hat\mu_2}{\hat{\mu_2}}, \eqn{\hat\sigma_2}{\hat{\sigma_2}} )}
#'   \item{pihat}{estimation of mixture proportion \eqn{\hat\pi}{\hat{\pi}}}
#'   \item{se.muhat1,se.sigmahat1}{standard errors of
#'   \eqn{\hat\mu_1}{\hat{\mu_1}} and \eqn{\hat\sigma_1}{\hat{\sigma_1}}}
#'   \item{se.muhat2,se.sigmahat2}{standard errors of
#'   \eqn{\hat\mu_2}{\hat{\mu_2}} and \eqn{\hat\sigma_2}{\hat{\sigma_2}}}
#'   \item{se.pihat}{standard error of \eqn{\hat\pi}{\hat{\pi}}}
#'   \item{no.fails.comp1,no.fails.comp2}{number of failures from each mixture
#'   component} \item{no.cens.comp1,no.cens.comp2}{number of censored
#'   observations from each mixture component} \item{ll}{log-likelihood value}
#'
#' @details This function fits a two-component normal mixture model to
#'   a given progressive Type-II censored data.
#'
#'   It uses a two-layer EM algorithm for fitting the model. Generally speaking,
#'   the first layer estimates the mixture component latent variables, in the E-step, by finding their conditional expected
#'   values given the current parameter estimates and the data; and the second layer consists
#'   of another EM algorithm to estimate the missing censored data and eventually the parameters of interest.
#'   The layers are repeated until convergence achieved.
#'
#'
#'
#' @note See \code{\link[pcensmix]{pcgen}} for the definition of censored
#'   version of data.
#' @examples
#' ## Example 1: fit a mixture model to 'insulate' data
#' set.seed(107)
#' Pdat<- pcgen(r = 15, p = 0.6, data = insulate)
#' pcensmixR(Pdat, start = c(5, 3, 35, 20, 0.6))
#'
#' \dontrun{
#' ## Example 2: fit a mixture model to 'Systolic blood pressure' data
#' set.seed(1010)
#' pcensmixR(Pdat = pcgen(360, 0.35, blood$Systolic.BP),
#'                start = c(120, 15, 150, 20, 0.6))}
#'
#' @rdname pcensmixR
#' @export
#' @seealso \code{\link{pcgen}}, \code{\link[pcensmix]{pcensmixSim}}
#' @author Lida Fallah, John Hinde
#'
#' Maintainer: Lida Fallah <l.fallah22@gmail.com>
#' @importFrom stats dnorm pnorm rbinom rgamma rlnorm rnorm runif rweibull rcauchy
#' rbeta var
#' @importFrom utils tail


pcensmixR.pcgen<-function(Pdat, start, iteration = 1e+05, INERiter = 20, warn = FALSE, ...){

  if (!warn)
    options(warn = -1)

  if( any(c("matrix", "data.frame") %in% class(Pdat)) ){
    N=nrow(Pdat)
    W=Pdat[,2]; Z=Pdat[,1]
  }
  else if ( "pcgen" %in% class(Pdat) )
  {
    N=length(Pdat$original_data); W=Pdat$censoring_indicator ; Z=Pdat$censored_version_of_data
  }

  if(!any(c("pcgen", "matrix", "data.frame") %in% class(Pdat))){
    stop("Pdat must be of class 'pcgen','data.frame', or 'matrix'.")
  }

  W=Pdat$censoring_indicator ; Z=Pdat$censored_version_of_data
  mu_1 = start[1]; sigma_1 = start[2]; mu_2 = start[3]; sigma_2 = start[4]; pi = start[5]

  tryCatch({

    muhat1 = mu_1; muhat2=mu_2; muhat1_est=NULL; muhat2_est = NULL; pi_hat=pi; pi_hat_est=NULL
    sigmahat1 = sigma_1; sigmahat2=sigma_2; sigmahat1_est=NULL; sigmahat2_est=NULL; ll=NULL

    eps<- 0.95

    for (i in 1:iteration){

      f2_u=.f_u(eps, Pdat, muhat1, muhat2, sigmahat1, sigmahat2, pi_hat, muhat1_est, muhat2_est, sigmahat1_est, sigmahat2_est, pi_hat_est, ll, W, Z, INERiter, N)
      muhat1<-f2_u$muhat1
      muhat2<-f2_u$muhat2
      sigmahat1<-f2_u$sigmahat1
      sigmahat2<-f2_u$sigmahat2
      pi_hat<-f2_u$pi_hat
      l<-f2_u$l


      ExpU<-c(f2_u$ExpU)
      x_tilde<-c(f2_u$x_tilde)
      y_tilde<-c(f2_u$y_tilde)
      xSq_tilde<-c(f2_u$xSq_tilde)
      ySq_tilde<-c(f2_u$ySq_tilde)



      muhat1_est=c(muhat1_est, muhat1)
      muhat2_est=c(muhat2_est, muhat2)
      sigmahat1_est<- c(sigmahat1_est, sigmahat1)
      sigmahat2_est<- c(sigmahat2_est, sigmahat2)
      pi_hat_est<-c(pi_hat_est, pi_hat)
      ll<-c(ll, l)

      if(i>1){

        if( (abs(sigmahat2_est[i]-sigmahat2_est[(i-1)])<1e-3) & (abs(muhat1_est[i]-muhat1_est[(i-1)])<1e-3) & (abs(sigmahat1_est[i]-sigmahat1_est[(i-1)])<1e-3) & (abs(muhat2_est[i]-muhat2_est[(i-1)])<1e-3) & (abs(pi_hat_est[i]-pi_hat_est[(i-1)])<1e-3) & (abs(ll[i]-ll[(i-1)])<1e-5) )
          break

      }
    }

    B <- matrix(NA, ncol=5, nrow=5)
    B[1,1] = sum(ExpU/sigmahat1^2)
    B[1,2] = B[2,1] = sum(2*ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^3)
    B[2,2] = - sum( ExpU*(W*( 1/(sigmahat1^2) - 3*(Z - muhat1)^2/sigmahat1^4) + (1-W)*( 1/(sigmahat1^2) - 3*(xSq_tilde -2*muhat1*x_tilde + muhat1^2)/sigmahat1^4)) )
    B[3,3] = sum((1-ExpU)/sigmahat2^2)
    B[3,4] = B[4,3]=sum(2*(1-ExpU)*(W*(Z-muhat2)+(1-W)*(y_tilde-muhat2))/sigmahat2^3 )
    B[4,4] =  sum( -(1-ExpU)*(W*( 1/(sigmahat2^2) - 3*(Z - muhat2)^2/sigmahat2^4) + (1-W)*( 1/(sigmahat2^2) - 3*(ySq_tilde -2*muhat2*y_tilde + muhat2^2)/sigmahat2^4)) )
    B[5,5] = sum(ExpU/pi_hat^2 + (1-ExpU)/(1-pi_hat)^2)
    B[is.na(B)] <- 0

    E_S.St <- matrix(NA, ncol=5, nrow=5)
    g1 <- (xSq_tilde -2*muhat1*x_tilde + muhat1^2)
    E_S.St[1,1] = sum(ExpU*(W^2*(x_tilde-muhat1)^2 + (1-W)^2*g1)/sigmahat1^4)

    xCub_tilde <- rep(0,length(W))
    for(j in 1:length(W)){
      if(W[j]==0 & (Z[j]-muhat1) < 6*sigmahat1){
        xCub_tilde[j]<- W[j]*(Z[j]^3) + (1-W[j])*( muhat1^3 + 3*sigmahat1^2*muhat1 + (2*sigmahat1^3 + sigmahat1*(Z[j]^2 + muhat1^2 + muhat1*Z[j]) )*dnorm((Z[j]-muhat1)/sigmahat1, 0, 1)/(1-pnorm((Z[j]-muhat1)/sigmahat1, 0, 1)) )
      }else{ xCub_tilde[j]<- W[j]*(Z[j]^3) }    }

    g2 <- (xCub_tilde-muhat1^3-3*muhat1*xSq_tilde + 3*muhat1^2*x_tilde)

    E_S.St[1,2] = E_S.St[2,1] = sum(ExpU*(W^2*(Z-muhat1)*((-1/(sigmahat1) + (Z-muhat1)^2/sigmahat1^3)) + (1-W)^2*(-(x_tilde-muhat1)/sigmahat1 + g2/sigmahat1^3))/sigmahat1^2)

    xforth_tilde <- rep(0,length(W))

    for(j in 1:length(W)){
      if(W[j]==0 & (Z[j]-muhat1) < 6*sigmahat1){
        xforth_tilde[j]<-W[j]*(Z[j]^4) + (1-W[j])*( ( 3*sigmahat1^4 + 6*muhat1^2*sigmahat1^2 + muhat1^4 ) + (7*sigmahat1^3*muhat1 + muhat1^3*sigmahat1 + sigmahat1^3*Z[j]+ muhat1^2*sigmahat1*Z[j] +muhat1*sigmahat1*Z[j]^2 + sigmahat1*Z[j]^3)*dnorm((Z[j]-muhat1)/sigmahat1, 0, 1)/(1-pnorm((Z[j]-muhat1)/sigmahat1, 0, 1)) )
      }else{xforth_tilde[j]<-W[j]*(Z[j]^4)}   }


    g3 <- xforth_tilde - 4*muhat1*xCub_tilde + 6*xSq_tilde*muhat1^2 -4*x_tilde*muhat1^3 + muhat1^4


    E_S.St[2,2] = sum( ExpU*( W^2*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3)^2  + (1-W)^2*(1/sigmahat1^2 + g3/sigmahat1^6 - 2*g1/sigmahat1^4) ) )

    g4 <- (ySq_tilde -2*muhat2*y_tilde + muhat2^2)
    E_S.St[3,3] = sum((1-ExpU)*(W^2*(Z-muhat2)^2 + (1-W)^2*g4)/sigmahat2^4)

    yforth_tilde <- rep(0,length(W))

    for(j in 1:length(W)){

      if(W[j]==0 & (Z[j]-muhat2)< 6*sigmahat2){
        yforth_tilde[j]<-W[j]*(Z[j]^4) + (1-W[j])*( ( 3*sigmahat2^4 + 6*muhat2^2*sigmahat2^2 + muhat2^4 ) + (5*sigmahat2^3*muhat2 + muhat2^3*sigmahat2 + 3*sigmahat2^3*Z[j]+ muhat2^2*sigmahat2*Z[j] +muhat2*sigmahat2*Z[j]^2 + sigmahat2*Z[j]^3)*dnorm((Z[j]-muhat2)/sigmahat2, 0, 1)/(1-pnorm((Z[j]-muhat2)/sigmahat2, 0, 1)) )
      }else{yforth_tilde[j]<-W[j]*(Z[j]^4)}    }

    yCub_tilde <- rep(0,length(W))

    for(j in 1:length(W)){
      if(W[j]==0 & (Z[j]-muhat2) < 6*sigmahat2){
        yCub_tilde[j]<- W[j]*(Z[j]^3) + (1-W[j])*( muhat2^3 + 3*sigmahat2^2*muhat2 + (2*sigmahat2^3 + sigmahat2*(Z[j]^2 + muhat2^2 + muhat2*Z[j]) )*dnorm((Z[j]-muhat2)/sigmahat2, 0, 1)/(1-pnorm((Z[j]-muhat2)/sigmahat2, 0, 1)) )
      }else{yCub_tilde[j]<- W[j]*(Z[j]^3) }     }

    g5 <- yforth_tilde - 4*muhat2*yCub_tilde + 6*ySq_tilde*muhat2^2 -4*y_tilde*muhat2^3 + muhat2^4
    E_S.St[4,4] = sum( (1-ExpU)*( W^2*(-1/sigmahat2 + (Z-muhat2)^2/(sigmahat2^3))^2  + (1-W)^2*(1/(sigmahat2^2) + g5/(sigmahat2^6) - 2*g4/(sigmahat2^4)) ) )

    g6 <- (yCub_tilde-muhat2^3-3*muhat2*ySq_tilde + 3*muhat2^2*y_tilde)

    E_S.St[3,4] = E_S.St[4,3]=sum((1-ExpU)*(W^2*(Z-muhat2)*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)^2*(-(y_tilde-muhat2)/sigmahat2 + g6/sigmahat2^3))/sigmahat2^2)
    E_S.St[1,5] = E_S.St[5,1]=sum(ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/(pi_hat*sigmahat1^2))
    E_S.St[2,5] = E_S.St[5,2]=sum((ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3)) )/pi_hat)
    E_S.St[3,5] = E_S.St[5,3]=sum((-(1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/((1-pi_hat)*sigmahat2^2))
    E_S.St[4,5] = E_S.St[5,4]=sum(-(1-ExpU)*( W*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3))/(1-pi_hat) )


    E_S.St[5,5] = sum(ExpU/pi_hat^2 + (1-ExpU)/(1-pi_hat)^2)

    E_S.St[is.na(E_S.St)] <- 0

    E_S.E_St <- matrix(NA, ncol=5, nrow=5)
    E_S.E_St[1,1] = sum( (ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^2)^2  )

    E_S.E_St[1,2] = E_S.E_St[2,1] = sum( (ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^2)*( ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3) ) ) )

    E_S.E_St[2,2] = sum( (ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3)) )^2 )

    E_S.E_St[3,1] = E_S.E_St[1,3] = sum( (ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^2)*( ((1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/sigmahat2^2 ) )

    E_S.E_St[3,3] = sum( (((1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/sigmahat2^2)^2 )

    E_S.E_St[2,3] = E_S.E_St[3,2] = sum( (ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3)) )*(((1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/sigmahat2^2)  )

    E_S.E_St[2,4] = E_S.E_St[4,2] = sum( (ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3)) )*((1-ExpU)*( W*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3)) ) )

    E_S.E_St[1,4] = E_S.E_St[4,1] = sum( (ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^2)*((1-ExpU)*( W*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3)) ) )

    E_S.E_St[3,4] = E_S.E_St[4,3] = sum( (((1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/sigmahat2^2)*((1-ExpU)*( W*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3)) )  )

    E_S.E_St[4,4] = sum( ((1-ExpU)*( W*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3)) )^2 )

    E_S.E_St[5,5] = sum((ExpU/pi_hat - (1-ExpU)/(1-pi_hat))^2)

    E_S.E_St[1,5] = E_S.E_St[5,1] = sum( (ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^2)*(ExpU/pi_hat - (1-ExpU)/(1-pi_hat)) )

    E_S.E_St[2,5] = E_S.E_St[5,2] = sum( (ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3)) )*(ExpU/pi_hat - (1-ExpU)/(1-pi_hat)) )

    E_S.E_St[3,5] = E_S.E_St[5,3] = sum( (((1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/sigmahat2^2)*(ExpU/pi_hat - (1-ExpU)/(1-pi_hat)) )

    E_S.E_St[4,5] = E_S.E_St[5,4] = sum( ((1-ExpU)*( W*(-1/sigmahat2 + (y_tilde-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3)) )*(ExpU/pi_hat - (1-ExpU)/(1-pi_hat)) )

    I_y <- B - E_S.St + E_S.E_St

    se.lu_mixtpcens.em <- sqrt(diag(solve(I_y)))

    n1_mixtpcens <- sum(W*ExpU)
    n2_mixtpcens <- sum(W*(1-ExpU))
    r1_mixtpcens <- sum((1-W)*ExpU)
    r2_mixtpcens <- sum((1-W)*(1-ExpU))

    variability_comp1 <- var((W*ExpU))
    variability_comp2 <- var((W*(1-ExpU)))

    res4<-data.frame('muhat1'=muhat1, 'sigmahat1'=sigmahat1, 'muhat2'=muhat2, 'sigmahat2'=sigmahat2, 'pihat'=pi_hat, 'se.muhat1'=se.lu_mixtpcens.em[1],  'se.sigmahat1'=se.lu_mixtpcens.em[2], 'se.muhat2'=se.lu_mixtpcens.em[3],'se.sigmahat2'=se.lu_mixtpcens.em[4] , 'se.pihat'= se.lu_mixtpcens.em[5], 'no.fails.comp1'=n1_mixtpcens, 'no.fails.comp2'=n2_mixtpcens, 'no.cens.comp1'=r1_mixtpcens, 'no.cens.comp2'=r2_mixtpcens, 'variability_comp1'=variability_comp1, 'variability_comp2'=variability_comp2, 'll'=l)
    res<- t(res4)

  }, error=function(err) {

    if(length(err) != 0)
    message("ERROR: singularity in the likelihood... Try refitting with different starting values.")
    else
     return(res)

  })

  if(exists("res"))
    return(res)
}



#' Fitting a Normal Mixture Model to a Simulated Progressive Type-II Censored
#' Data Using EM Algorithm
#'
#' {}
#'
#' @param Pdat an object of class \code{"pcgen"} created by function
#'   \code{\link[pcensmix]{pcgen}} or a two-column matrix (or data.frame)
#'   with first column
#'   giving a vector of censored version of a two-component mixed normal data,
#'   and the other one indicating the censoring status associated with them (1
#'   if not censored, otherwise zero).
#' @param ... additinal arguments to pass by.
#'
#' @export


pcensmixSim<-function (Pdat, ...)
{
  UseMethod("pcensmixSim")
}




#' Fitting a Normal Mixture Model to a Simulated Progressive Type-II Censored
#' Data Using EM Algorithm
#'
#' This function fits a normal mixture model
#' to progressive Type-II censored mixture data by dealing with the two aspects of
#' missing data, latent mixture components and the censored data, using a maximum
#' likelihood estimation through a constrained two-layer EM algorithm.
#'
#' @param r desired number of failures to observe.
#' @param p a parameter controlling the amount of censoring. The action of
#'   censoring individuals after each failure occurs with probabilty \code{p}
#'   from binomial distribution at each stage. If \code{p = 0}, there will be no
#'   censoring.
#' @param param a numeric vector; used as starting values for the EM and
#'   simulating a new data to replace in case of happening singularity in the
#'   likelihood.
#' @param iteration the maximum number of required iteration for the EM
#'   algorithm until convergence-- default value is 1e+05.
#' @param INERiter the maximum number of required iteration for the second EM
#'   algorithm-- default is 20.
#'
#' @return \code{pcensmixSim} gives an object of class \code{data.frame}
#'   containing the following components: \item{muhat1,sigmahat1}{component one
#'   parameter estimates (\eqn{\hat\mu_1}{\hat{\mu_1}},
#'   \eqn{\hat\sigma_1}{\hat{\sigma_1}} )} \item{muhat2,sigmahat2}{component two
#'   parameter estimates (\eqn{\hat\mu_2}{\hat{\mu_2}},
#'   \eqn{\hat{\sigma_2}}{\hat{\sigma_2} )}} \item{pihat}{estimation of mixture
#'   proportion \eqn{\hat\pi}{\hat{\pi}}} \item{se.muhat1,se.sigmahat1}{standard
#'   errors of \eqn{\hat\mu_1}{\hat{\mu_1}} and \eqn{\hat{\sigma_1}}}
#'   \item{se.muhat2,se.sigmahat2}{standard errors of
#'   \eqn{\hat\mu_2}{\hat{\mu_2}} and \eqn{\hat\sigma_2}{\hat{\sigma_2}}}
#'   \item{se.pihat}{standard error of \eqn{\hat\pi}{\hat{\pi}}}
#'   \item{no.fails.comp1,no.fails.comp2}{number of failures from each mixture
#'   component} \item{no.cens.comp1,no.cens.comp2}{number of censored
#'   observations from each mixture component} \item{ll}{log-likelihood value}
#'   \item{datachange_flag}{\code{TRUE} if data has been replaced by a newly
#'   generated one}
#'
#' @details This function fits a two-component normal mixture model to simulated
#'   progressive Type-II censored data with density function \deqn{\pi
#'   (\frac{1}{\sigma_1})\,  \phi\!\! \left[\frac{(z - \mu_1)}{\sigma_1}\right]
#'   + (1 - \pi) (\frac{1}{\sigma_2})\,  \phi\!\! \left[\frac{(z -
#'   \mu_2)}{\sigma_2}\right]}{\pi (1/ \sigma_1) \phi[ (z - \mu_1) / \sigma_1]
#'   + (1 - \pi) (1/ \sigma_2) \phi[ (z - \mu_2) / \sigma_2]} where \eqn{\phi}
#'   is the standard normal density.
#'
#'   It uses a constrained two-layer EM algorithm to deal with the two forms of
#'   missing data: the censored survival times and the mixture component labels.
#'   Given the EM algorithm is at a particular iteration: (i) first, in the
#'   E-step it obtains the mixture component indicator estimates given the
#'   current parameter estimates and the observed data. (ii) Next, for
#'   re-estimation of the unknown parameters, a new EM algorithm is nested in
#'   the M-step of the initial EM algorithm to deal with the estimation of the
#'   missing censored survival times and consequently building the maximum
#'   likelihood equations. These steps are repeated until the model converges.
#' @note \itemize{ \item In fitting the model, to overcome the problem of
#'   singularity and model non-identifiability that might happen in some cases
#'   depending on the true means and standard deviations of the components, we
#'   use the constraint proposed by Hathaway (1985). Based on this, the
#'   ratios of the standard deviations are considered to be greater than a
#'   pre-fixed constant. We consider the constraints \eqn{\sigma_1/\sigma_2>0.1}
#'   and \eqn{\sigma_2/\sigma_1>0.1} which lead to obtain a parameter space with
#'   no singularities. In case this conditions doesn't hold, the data will be
#'   replaced by a new simulated one and \code{datachange_flag} will appear as
#'   \code{TRUE} in the output. \item See \code{\link[pcensmix]{pcgen}} for the
#'   definition of censored version of data.}
#' @examples
#' \dontrun{
#' set.seed(100)
#'
#' Pdat<- pcgen(r = 60, p = 0.3, data = mixgen(N = 100, dist1 = 'norm',
#'                  dist2 = 'norm', control = list(12, 2, 14, 5, 0.35)))
#' pcensmixSim(Pdat, r = 60, p = 0.3, param=c(12, 2, 14, 5, 0.35))}
#'
#' @rdname pcensmixSim
#' @export
#' @seealso \code{\link[pcensmix]{pcgen}}, \code{\link[pcensmix]{run_pcensmix}}
#' @author Lida Fallah, John Hinde
#'
#'   Maintainer: Lida Fallah <l.fallah22@gmail.com>
#' @importFrom stats dnorm pnorm rbinom rgamma rlnorm rnorm runif rweibull
#'   rcauchy rbeta var
#' @importFrom utils tail


pcensmixSim.pcgen<-function(Pdat, r, p, param, iteration = 1e+05, INERiter = 20, ...){

  if( any(c("matrix", "data.frame") %in% class(Pdat)) ){
    N = nrow(Pdat)
    W = Pdat[,2]; Z = Pdat[,1]
  }
  else if ( "pcgen" %in% class(Pdat) )
  {
    N = length(Pdat$original_data); W = Pdat$censoring_indicator ; Z = Pdat$censored_version_of_data
  }

  if(!any(c("pcgen", "matrix", "data.frame") %in% class(Pdat))){
    stop("Pdat must be of class 'pcgen','data.frame', or 'matrix'.")
  }

  W = Pdat$censoring_indicator ; Z = Pdat$censored_version_of_data

  mu_1 = param[1]; sigma_1 = param[2]; mu_2 = param[3]; sigma_2 = param[4]; pi = param[5]


  muhat1 = mu_1; muhat2=mu_2; muhat1_est=NULL; muhat2_est = NULL; pi_hat = pi; pi_hat_est = NULL
  sigmahat1 = sigma_1; sigmahat2 = sigma_2; sigmahat1_est = NULL; sigmahat2_est = NULL; ll = NULL

  datachange <- NULL
  eps <- 0.95

  for (i in 1:iteration){

    f2_u=.f_u(eps, Pdat, muhat1, muhat2, sigmahat1, sigmahat2, pi_hat, muhat1_est, muhat2_est, sigmahat1_est, sigmahat2_est, pi_hat_est, ll, W, Z, INERiter, N)
    muhat1 <- f2_u$muhat1
    muhat2 <- f2_u$muhat2
    sigmahat1 <- f2_u$sigmahat1
    sigmahat2 <- f2_u$sigmahat2
    pi_hat <- f2_u$pi_hat
    l <- f2_u$l


    ExpU <- c(f2_u$ExpU)
    x_tilde <- c(f2_u$x_tilde)
    y_tilde <- c(f2_u$y_tilde)
    xSq_tilde <- c(f2_u$xSq_tilde)
    ySq_tilde <- c(f2_u$ySq_tilde)


    if((sigmahat1/sigmahat2)> 0.1 & (sigmahat2/sigmahat1)> 0.1){
      datachange <- c(datachange, 0)
    }
    datachange2 <- datachange

    while((sigmahat1/sigmahat2)< 0.1 | (sigmahat2/sigmahat1)< 0.1){

      datachange <- c(datachange2, 1)

      muhat1 = mu_1; muhat2 = mu_2; pi_hat = pi
      sigmahat1 = sigma_1; sigmahat2 = sigma_2

      eps <- log(eps, base = i)


      mixt <- mixgen(N, dist1 = 'norm', dist2 = 'norm', control = list(mu_1, sigma_1, mu_2, sigma_2, pi) )
      Pdat <- pcgen(r, p, mixt)

      W = Pdat$censoring_indicator; Z = Pdat$censored_version_of_data

      f2_u = .f_u(eps, Pdat, muhat1, muhat2, sigmahat1, sigmahat2, pi_hat, muhat1_est, muhat2_est, sigmahat1_est, sigmahat2_est, pi_hat_est, ll, W, Z, INERiter, N)

      muhat1 <- f2_u$muhat1
      muhat2 <- f2_u$muhat2
      sigmahat1 <- f2_u$sigmahat1
      sigmahat2 <- f2_u$sigmahat2
      pi_hat <- f2_u$pi_hat

      l <- f2_u$l

      ExpU <- c(f2_u$ExpU)
      x_tilde <- c(f2_u$x_tilde)
      y_tilde <- c(f2_u$y_tilde)
      xSq_tilde <- c(f2_u$xSq_tilde)
      ySq_tilde <- c(f2_u$ySq_tilde)
    }

    muhat1_est = c(muhat1_est, muhat1)
    muhat2_est = c(muhat2_est, muhat2)
    sigmahat1_est <- c(sigmahat1_est, sigmahat1)
    sigmahat2_est <- c(sigmahat2_est, sigmahat2)
    pi_hat_est <- c(pi_hat_est, pi_hat)
    ll <- c(ll, l)

    if(i>1 & datachange[i] == 0){

      if( (abs(sigmahat2_est[i]-sigmahat2_est[(i-1)])<1e-3) & (abs(muhat1_est[i]-muhat1_est[(i-1)])<1e-3) & (abs(sigmahat1_est[i]-sigmahat1_est[(i-1)])<1e-3) & (abs(muhat2_est[i]-muhat2_est[(i-1)])<1e-3) & (abs(pi_hat_est[i]-pi_hat_est[(i-1)])<1e-3) & (abs(ll[i]-ll[(i-1)])<1e-5) )
        break

    }
  }

  datachange_flag<-any(datachange == 1)

  B <- matrix(NA, ncol=5, nrow=5)
  B[1,1] = sum(ExpU/sigmahat1^2)
  B[1,2] = B[2,1]=sum(2*ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^3)
  B[2,2] = - sum( ExpU*(W*( 1/(sigmahat1^2) - 3*(Z - muhat1)^2/sigmahat1^4) + (1-W)*( 1/(sigmahat1^2) - 3*(xSq_tilde -2*muhat1*x_tilde + muhat1^2)/sigmahat1^4)) )
  B[3,3] = sum((1-ExpU)/sigmahat2^2)
  B[3,4] = B[4,3] = sum(2*(1-ExpU)*(W*(Z-muhat2)+(1-W)*(y_tilde-muhat2))/sigmahat2^3 )
  B[4,4] = sum( -(1-ExpU)*(W*( 1/(sigmahat2^2) - 3*(Z - muhat2)^2/sigmahat2^4) + (1-W)*( 1/(sigmahat2^2) - 3*(ySq_tilde -2*muhat2*y_tilde + muhat2^2)/sigmahat2^4)) )
  B[5,5] = sum(ExpU/pi_hat^2 + (1-ExpU)/(1-pi_hat)^2)
  B[is.na(B)] <- 0

  E_S.St <- matrix(NA, ncol = 5, nrow = 5)
  g1 <- (xSq_tilde -2*muhat1*x_tilde + muhat1^2)
  E_S.St[1,1] = sum(ExpU*(W^2*(x_tilde-muhat1)^2 + (1-W)^2*g1)/sigmahat1^4)

  xCub_tilde <- rep(0,length(W))
  for(j in 1:length(W)){
    if(W[j] == 0 & (Z[j]-muhat1) < 6*sigmahat1){
      xCub_tilde[j] <- W[j]*(Z[j]^3) + (1-W[j])*( muhat1^3 + 3*sigmahat1^2*muhat1 + (2*sigmahat1^3 + sigmahat1*(Z[j]^2 + muhat1^2 + muhat1*Z[j]) )*dnorm((Z[j]-muhat1)/sigmahat1, 0, 1)/(1-pnorm((Z[j]-muhat1)/sigmahat1, 0, 1)) )
    }else{ xCub_tilde[j]<- W[j]*(Z[j]^3) }    }

  g2 <- (xCub_tilde-muhat1^3-3*muhat1*xSq_tilde + 3*muhat1^2*x_tilde)

  E_S.St[1,2] = E_S.St[2,1] = sum(ExpU*(W^2*(Z-muhat1)*((-1/(sigmahat1) + (Z-muhat1)^2/sigmahat1^3)) + (1-W)^2*(-(x_tilde-muhat1)/sigmahat1 + g2/sigmahat1^3))/sigmahat1^2)

  xforth_tilde <- rep(0,length(W))
  for(j in 1:length(W)){

    if(W[j] == 0 & (Z[j]-muhat1) < 6*sigmahat1){
      xforth_tilde[j]<-W[j]*(Z[j]^4) + (1-W[j])*( ( 3*sigmahat1^4 + 6*muhat1^2*sigmahat1^2 + muhat1^4 ) + (7*sigmahat1^3*muhat1 + muhat1^3*sigmahat1 + sigmahat1^3*Z[j]+ muhat1^2*sigmahat1*Z[j] +muhat1*sigmahat1*Z[j]^2 + sigmahat1*Z[j]^3)*dnorm((Z[j]-muhat1)/sigmahat1, 0, 1)/(1-pnorm((Z[j]-muhat1)/sigmahat1, 0, 1)) )
    }else{xforth_tilde[j] <- W[j]*(Z[j]^4)}   }

  g3 <- xforth_tilde - 4*muhat1*xCub_tilde + 6*xSq_tilde*muhat1^2 -4*x_tilde*muhat1^3 + muhat1^4


  E_S.St[2,2] = sum( ExpU*( W^2*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3)^2  + (1-W)^2*(1/sigmahat1^2 + g3/sigmahat1^6 - 2*g1/sigmahat1^4) ) )

  g4 <- (ySq_tilde -2*muhat2*y_tilde + muhat2^2)
  E_S.St[3,3] = sum((1-ExpU)*(W^2*(Z-muhat2)^2 + (1-W)^2*g4)/sigmahat2^4)

  yforth_tilde <- rep(0,length(W))

  for(j in 1:length(W)){
    if(W[j] == 0 & (Z[j]-muhat2)< 6*sigmahat2){
      yforth_tilde[j] <- W[j]*(Z[j]^4) + (1-W[j])*( ( 3*sigmahat2^4 + 6*muhat2^2*sigmahat2^2 + muhat2^4 ) + (7*sigmahat2^3*muhat2 + muhat2^3*sigmahat2 + sigmahat2^3*Z[j]+ muhat2^2*sigmahat2*Z[j] +muhat2*sigmahat2*Z[j]^2 + sigmahat2*Z[j]^3)*dnorm((Z[j]-muhat2)/sigmahat2, 0, 1)/(1-pnorm((Z[j]-muhat2)/sigmahat2, 0, 1)) )
    }else{yforth_tilde[j] <- W[j]*(Z[j]^4)}    }

  yCub_tilde <- rep(0,length(W))

  for(j in 1:length(W)){
    if(W[j]==0 & (Z[j]-muhat2) < 6*sigmahat2){
      yCub_tilde[j] <- W[j]*(Z[j]^3) + (1-W[j])*( muhat2^3 + 3*sigmahat2^2*muhat2 + (2*sigmahat2^3 + sigmahat2*(Z[j]^2 + muhat2^2 + muhat2*Z[j]) )*dnorm((Z[j]-muhat2)/sigmahat2, 0, 1)/(1-pnorm((Z[j]-muhat2)/sigmahat2, 0, 1)) )
    }else{yCub_tilde[j]<- W[j]*(Z[j]^3) }     }

  g5 <- yforth_tilde - 4*muhat2*yCub_tilde + 6*ySq_tilde*muhat2^2 -4*y_tilde*muhat2^3 + muhat2^4
  E_S.St[4,4] = sum( (1-ExpU)*( W^2*(-1/sigmahat2 + (Z-muhat2)^2/(sigmahat2^3))^2  + (1-W)^2*(1/(sigmahat2^2) + g5/(sigmahat2^6) - 2*g4/(sigmahat2^4)) ) )

  g6 <- (yCub_tilde-muhat2^3-3*muhat2*ySq_tilde + 3*muhat2^2*y_tilde)

  E_S.St[3,4] = E_S.St[4,3] = sum((1-ExpU)*(W^2*(Z-muhat2)*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)^2*(-(y_tilde-muhat2)/sigmahat2 + g6/sigmahat2^3))/sigmahat2^2)
  E_S.St[1,5] = E_S.St[5,1] = sum(ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/(pi_hat*sigmahat1^2))
  E_S.St[2,5] = E_S.St[5,2] = sum((ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3)) )/pi_hat)
  E_S.St[3,5] = E_S.St[5,3] = sum((-(1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/((1-pi_hat)*sigmahat2^2))
  E_S.St[4,5] = E_S.St[5,4] = sum(-(1-ExpU)*( W*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3))/(1-pi_hat) )


  E_S.St[5,5] = sum(ExpU/pi_hat^2 + (1-ExpU)/(1-pi_hat)^2)

  E_S.St[is.na(E_S.St)] <- 0

  E_S.E_St <- matrix(NA, ncol=5, nrow=5)
  E_S.E_St[1,1] = sum( (ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^2)^2  )

  E_S.E_St[1,2] = E_S.E_St[2,1] = sum( (ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^2)*( ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3) ) ) )

  E_S.E_St[2,2] = sum( (ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3)) )^2 )

  E_S.E_St[3,1] = E_S.E_St[1,3] = sum( (ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^2)*( ((1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/sigmahat2^2 ) )

  E_S.E_St[3,3] = sum( (((1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/sigmahat2^2)^2 )

  E_S.E_St[2,3] = E_S.E_St[3,2] = sum( (ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3)) )*(((1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/sigmahat2^2)  )

  E_S.E_St[2,4] = E_S.E_St[4,2] = sum( (ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3)) )*((1-ExpU)*( W*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3)) ) )

  E_S.E_St[1,4] = E_S.E_St[4,1] = sum( (ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^2)*((1-ExpU)*( W*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3)) ) )

  E_S.E_St[3,4] = E_S.E_St[4,3] = sum( (((1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/sigmahat2^2)*((1-ExpU)*( W*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3)) )  )

  E_S.E_St[4,4] = sum( ((1-ExpU)*( W*(-1/sigmahat2 + (Z-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3)) )^2 )

  E_S.E_St[5,5] = sum((ExpU/pi_hat - (1-ExpU)/(1-pi_hat))^2)

  E_S.E_St[1,5] = E_S.E_St[5,1] = sum( (ExpU*(W*(Z-muhat1)+(1-W)*(x_tilde-muhat1))/sigmahat1^2)*(ExpU/pi_hat - (1-ExpU)/(1-pi_hat)) )

  E_S.E_St[2,5] = E_S.E_St[5,2] = sum( (ExpU*( W*(-1/sigmahat1 + (Z-muhat1)^2/sigmahat1^3) + (1-W)*(-1/sigmahat1 + g1/sigmahat1^3)) )*(ExpU/pi_hat - (1-ExpU)/(1-pi_hat)) )

  E_S.E_St[3,5] = E_S.E_St[5,3] = sum( (((1-ExpU)*( W*(Z-muhat2)+(1-W)*(y_tilde-muhat2) ) )/sigmahat2^2)*(ExpU/pi_hat - (1-ExpU)/(1-pi_hat)) )

  E_S.E_St[4,5] = E_S.E_St[5,4] = sum( ((1-ExpU)*( W*(-1/sigmahat2 + (y_tilde-muhat2)^2/sigmahat2^3) + (1-W)*(-1/sigmahat2 + g4/sigmahat2^3)) )*(ExpU/pi_hat - (1-ExpU)/(1-pi_hat)) )

  I_y <- B - E_S.St + E_S.E_St

  se.lu_mixtpcens.em <- sqrt(diag(solve(I_y)))

  n1_mixtpcens <- sum(W*ExpU)
  n2_mixtpcens <- sum(W*(1-ExpU))
  r1_mixtpcens <- sum((1-W)*ExpU)
  r2_mixtpcens <- sum((1-W)*(1-ExpU))

  variability_comp1 <- var((W*ExpU))
  variability_comp2 <- var((W*(1-ExpU)))


  res4 <- data.frame('muhat1'=muhat1, 'sigmahat1'=sigmahat1, 'muhat2'=muhat2, 'sigmahat2'=sigmahat2, 'pihat'=pi_hat, 'se.muhat1'=se.lu_mixtpcens.em[1],  'se.sigmahat1'=se.lu_mixtpcens.em[2], 'se.muhat2'=se.lu_mixtpcens.em[3],'se.sigmahat2'=se.lu_mixtpcens.em[4] , 'se.pihat'= se.lu_mixtpcens.em[5],  'no.fails.comp1'=n1_mixtpcens, 'no.fails.comp2'=n2_mixtpcens, 'no.cens.comp1'=r1_mixtpcens, 'no.cens.comp2'=r2_mixtpcens, 'variability_comp1'=variability_comp1, 'variability_comp2'=variability_comp2, 'll'=l, 'datachange_flag'=datachange_flag)

  t(res4)
}



#'Generating Progressively Type-II Censored Mixture Data and Fitting a Model
#'
#'This function implements an algorithm using the
#'\code{\link[pcensmix]{mixgen}}, \code{\link[pcensmix]{pcgen}} and
#'\code{\link[pcensmix]{pcensmixSim}} functions to generate data and fit a model
#'using EM algorithm with a specified number of iterations.
#'
#'@param N population size.
#'@param r required number of failures to observe.
#'@param p a parameter controlling the amount of censoring. The action of
#'  censoring individuals after each failure occurs with probabilty \code{p}
#'  from a binomial distribution at each stage. If \code{p = 0}, there will be
#'  no censoring.
#'@param param a numeric vector; used as starting values for the EM and
#'  simulating a new data to replace in case of happening singularity in the
#'  likelihood.
#'@param repetition the required number of repetition of the algorithm-- default
#'  is 100.
#'
#'@return It returns the parameter estimates given by
#'  \code{\link[pcensmix]{pcensmixSim}} with the desired number of repetitions.
#'  In each repetition it generates a new normal mixture progressive Type-II
#'  censored dataset from the same true parameter values and fits a model.
#'
#'
#' @examples
#' \dontrun{
#'
#' ## Example 1: with very well separated mixture components
#' set.seed(3)
#' f1 <- run_pcensmix(N = 160, r = 120, p = 0.3, param = c(10, 2, 25, 4, 0.3), repetition = 100)
#' colMeans(f1)
#'
#' ## Example 2.
#' set.seed(150)
#' f2 <- run_pcensmix(N = 160, r = 130, p = 0.35, param = c(10, 2, 17, 4, 0.3), repetition = 100)
#' colMeans(f2)
#'
#' ## Example 3.
#' set.seed(20)
#' f3 <- run_pcensmix(N = 160, r = 130, p = 0.3, param = c(20, 6, 22, 12, 0.6), repetition = 100)
#' colMeans(f3)
#'
#' }
#'
#'@export
#'@seealso \code{\link[pcensmix]{pcgen}}, \code{\link[pcensmix]{pcensmixSim}}, \code{\link[pcensmix]{mixgen}}
#'@author Lida Fallah, John Hinde
#'
#'Maintainer: Lida Fallah <l.fallah22@gmail.com>
#'@importFrom stats dnorm pnorm rbinom rgamma rlnorm rnorm runif rweibull
#'  rcauchy rbeta var
#'@importFrom utils tail

run_pcensmix<-function(N, r, p, param, repetition = 100){

  mu_1 = param[1]; sigma_1 = param[2]; mu_2 = param[3]; sigma_2 = param[4]; pi = param[5]

  PTIICres<-NULL;

  for(k in 1:repetition){


    mixt<-mixgen(N, dist1 = 'norm', dist2 = 'norm', control = list(mu_1, sigma_1, mu_2, sigma_2, pi) )
    Pdat<-pcgen(r, p, mixt)

    PTIIres<-t(pcensmixSim.pcgen(Pdat, r, p, param, iteration=100000, INERiter = 20))

    PTIICres<-rbind(PTIICres, PTIIres)

  }

  PTIICres
}

