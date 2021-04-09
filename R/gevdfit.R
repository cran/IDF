# This file contains the functions:
# - gev.d.fit, gev.d.init for fitting
# - gev.d.diag for diagnostic plots
# - gev.d.params for calculation of parameters
# and the documentation of the example data

#### gev.d.fit ####

#' @title Maximum-likelihood Fitting of the duration-dependent GEV Distribution
#' @description Modified \code{\link[ismev]{gev.fit}} function for Maximum-likelihood fitting
#' for the duration-dependent generalized extreme
#' value distribution, following Koutsoyiannis et al. (1998), including generalized linear
#' modeling of each parameter.
#' @param xdat A vector containing maxima for different durations.
#' This can be obtained from \code{\link{IDF.agg}}.
#' @param ds A vector of aggregation levels corresponding to the maxima in xdat.
#' 1/60 corresponds to 1 minute, 1 corresponds to 1 hour.
#' @param ydat A matrix of covariates for generalized linear modeling of the parameters
#' (or NULL (the default) for stationary fitting). The number of rows should be the same as the
#' length of xdat.
#' @param  mutl,sigma0l,xil,thetal,etal,taul,eta2l Numeric vectors of integers, giving the columns of ydat that contain
#'  covariates for generalized linear modeling of the parameters (or NULL (the default)
#'  if the corresponding parameter is stationary).
#'  Parameters are: modified location, scale offset, shape, duration offset, duration exponent, respectively.
#' @param mutlink,sigma0link,xilink,thetalink,etalink,eta2link,taulink Link functions for generalized linear
#' modeling of the parameters, created with \code{\link{make.link}}. The default is \code{make.link("identity")}.
#' @param init.vals list, giving initial values for all or some parameters
#' (order: mut, sigma0, xi, theta, eta, eta2, tau). If one of those parameters shall not be used (see theta_zero, eta2_zero, tau_zero),
#' the number of parameters has to be reduced accordingly. If some or all given values in init.vals are NA or 
#' no init.vals at all is declared (the default), initial parameters are obtained
#' internally by fitting the GEV separately for each duration and applying a linear model to obtain the
#' duration dependency of the location and shape parameter.
#' Initial values for covariate parameters are assumed as 0 if not given.
#' @param theta_zero Logical value, indicating whether theta should be estimated (FALSE, the default) or
#' should stay zero.
#' @param tau_zero,eta2_zero Logical values, indicating whether tau,eta2 should be estimated (TRUE, the default).
#' @param show Logical; if TRUE (the default), print details of the fit.
#' @param method The optimization method used in \code{\link{optim}}.
#' @param maxit The maximum number of iterations.
#' @param ... Other control parameters for the optimization.
#' @return A list containing the following components.
#' A subset of these components are printed after the fit.
#' If \code{show} is TRUE, then assuming that successful convergence is indicated,
#' the components nllh, mle and se are always printed.
#' \item{nllh}{single numeric giving the negative log-likelihood value}
#' \item{mle}{numeric vector giving the MLE's for the modified location, scale_0, shape,
#' duration offset and duration exponent, resp. If requested, contains also second duration exponent and intensity-offset}
#' \item{se}{numeric vector giving the standard errors for the MLE's (in the same order)}
#' \item{trans}{A logical indicator for a non-stationary fit.}
#' \item{model}{A list with components mutl, sigma0l, xil, thetal and etal. If requested, contains also eta2l and taul}
#' \item{link}{A character vector giving inverse link functions.}
#' \item{conv}{The convergence code, taken from the list returned by \code{\link{optim}}.
#' A zero indicates successful convergence.}
#' \item{data}{data is standardized to standard Gumbel.}
#' \item{cov}{The covariance matrix.}
#' \item{vals}{Parameter values for every data point.}
#' \item{init.vals}{Initial values that were used.}
#' \item{ds}{Durations for every data point.}
#' @details For details on the d-GEV and the parameter definitions, see \link{IDF-package}.
#' @seealso \code{\link{IDF-package}}, \code{\link{IDF.agg}}, \code{\link{gev.fit}}, \code{\link{optim}}
#' @export
#' @importFrom stats optim
#' @importFrom stats make.link
#'
#' @examples
#' # sampled random data from d-gev with covariates
#' # GEV parameters:
#' # mut = 4 + 0.2*cov1 +0.5*cov2
#' # sigma0 = 2+0.5*cov1
#' # xi = 0.5
#' # theta = 0
#' # eta = 0.5
#' # eta2 = 0.5
#' # tau = 0
#'
#' data('example',package ='IDF')
#'
#' gev.d.fit(xdat=example$dat,ds = example$d,ydat=as.matrix(example[,c('cov1','cov2')])
#' ,mutl=c(1,2),sigma0l=1)

gev.d.fit<-
  function(xdat, ds, ydat = NULL, mutl = NULL, sigma0l = NULL, xil = NULL, thetal = NULL, etal = NULL, taul = NULL, eta2l = NULL,
           mutlink = make.link("identity"), sigma0link = make.link("identity"), xilink = make.link("identity"),
           thetalink = make.link("identity"), etalink = make.link("identity"), taulink = make.link("identity"), eta2link = make.link("identity"),
           init.vals = NULL, theta_zero = FALSE, tau_zero=TRUE, eta2_zero=TRUE,
           show = TRUE, method = "Nelder-Mead", maxit = 10000, ...)
  {
    if (length(xdat) != length(ds)) {
      stop(paste0('The length of xdat is ',length(xdat),', but the length of ds is ',length(ds),'.'))
    }
    z <- list()
    # number of parameters (betas) to estimate for each parameter:
    npmu <- length(mutl) + 1
    npsc <- length(sigma0l) + 1
    npsh <- length(xil) + 1
    npth <- ifelse(!theta_zero,length(thetal) + 1,0)
    npet <- length(etal) + 1
    npta <- ifelse(!tau_zero,  length(taul)   + 1,0)
    npe2 <- ifelse(!eta2_zero,  length(eta2l)   + 1,0)
    z$trans <- FALSE  # indicates if fit is non-stationary
    z$model <- list(mutl, sigma0l, xil, thetal, etal, eta2l, taul)
    z$link <- list(mutlink=mutlink, sigma0link=sigma0link, xilink=xilink, thetalink=thetalink, etalink=etalink, eta2link=eta2link, taulink=taulink)

    # test for NA values:
    if(any(is.na(xdat))) stop('xdat contains NA values. NA values need to be removed first.')
    # test for finite values:
    if(any(is.infinite(xdat))) stop('xdat contains non finite values. Inf and -Inf need to be removed first.')

    # test if covariates matrix is given correctly
    npar <- max(sapply(z$model,function(x){return(ifelse(is.null(x),0,max(x)))}))
    if(any(npar>ncol(ydat),npar>0 & is.null(ydat)))stop("Not enough columns in covariates matrix 'ydat'.")

    # initial values
    init.necessary.length = 4 + as.numeric(!theta_zero) + as.numeric(!eta2_zero) + as.numeric(!tau_zero)  # mut, sigma0, xi, theta, eta, eta2, tau
    if (is.null(init.vals)) {init.vals = as.list(rep(NA,init.necessary.length))
    }else{init.vals = as.list(init.vals)}
    
    if(length(init.vals)!=init.necessary.length | !is.list(init.vals)) {
      print(paste0('Parameter init.vals is not used, because it is no list of length ',init.necessary.length,'.'))
      init.vals <- gev.d.init(xdat,ds,z$link)
      
    }else{ # length of given values is correct

      # name given initial values
      names1=c('mu','sigma','xi')                 # standard set of parameters
      if (!theta_zero){names1=c(names1,'theta')}  # add theta (in case)
      names1=c(names1,'eta')                      # add eta   (always)
      if (!eta2_zero){names1=c(names1,'eta2')}    # add eta2  (in case)
      if (!tau_zero){names1=c(names1,'tau')}      # add tau   (in case)
      names(init.vals) <- names1
      # add missing initial value (fixed internal number of parameters: 7)
      if (theta_zero) init.vals$theta = 0
      if (eta2_zero) init.vals$eta2 = init.vals$eta
      if (tau_zero) init.vals$tau = 0
      init.vals=init.vals[c("mu","sigma","xi","theta","eta","eta2","tau")]
      iv=init.vals
      init.vals = list(mu=iv$mu,sigma=iv$sigma,xi=iv$xi,theta=iv$theta,eta=iv$eta,eta2=iv$eta2,tau=iv$tau)
      
      if(!any(is.na(init.vals))){ #all initial values are given
        # do nothing
      }else if(any(!is.na(init.vals))) { #some initial values are given
        #if (!eta2_zero) print("autmoatic inital value setting not implemented yet for multiscaling (eta2_zero=FALSE)")
        init.vals.user <- init.vals
        init.vals <- gev.d.init(xdat,ds,z$link) #calculate init.vals using gev.d.init
        for (i in 1:length(init.vals)){ #overwrite the calculated initial values with the ones given by the user
          if(!is.na(init.vals.user[[names(init.vals.user)[i]]])) {
            init.vals[[names(init.vals.user)[i]]]<-init.vals.user[[names(init.vals.user)[i]]]
          } 
        }
      }else{ #no initial values are given
        #if (!eta2_zero) print("autmoatic inital value setting not implemented yet for multiscaling (eta2_zero=FALSE)")
        init.vals <- gev.d.init(xdat,ds,z$link)
      }
    } 

    # generate covariates matrices:
    if (is.null(mutl)) { #stationary
      mumat <- as.matrix(rep(1, length(xdat)))
      muinit <- init.vals$mu
    }else { #non stationary
      z$trans <- TRUE
      mumat <- cbind(rep(1, length(xdat)), ydat[, mutl])
      muinit <- c(init.vals$mu, rep(0, length(mutl)))[1:npmu] #fill with 0s to length npmu
    }
    if (is.null(sigma0l)) {
      sigmat <- as.matrix(rep(1, length(xdat)))
      siginit <- init.vals$sigma
    }else {
      z$trans <- TRUE
      sigmat <- cbind(rep(1, length(xdat)), ydat[, sigma0l])
      siginit <- c(init.vals$sigma, rep(0, length(sigma0l)))[1:npsc]
    }
    if (is.null(xil)) {
      shmat <- as.matrix(rep(1, length(xdat)))
      shinit <- init.vals$xi
    }else {
      z$trans <- TRUE
      shmat <- cbind(rep(1, length(xdat)), ydat[, xil])
      shinit <- c(init.vals$xi, rep(0, length(xil)))[1:npsh]
    }
    if (is.null(thetal)) {
      thmat <- as.matrix(rep(1, length(xdat)))
      thetainit <- init.vals$theta
    }else {
      z$trans <- TRUE
      thmat <- cbind(rep(1, length(xdat)), ydat[, thetal])
      thetainit <- c(init.vals$theta, rep(0, length(thetal)))[1:npth]
    }
    if (is.null(etal)) {
      etmat <- as.matrix(rep(1, length(xdat)))
      etainit <- init.vals$eta
    }else {
      z$trans <- TRUE
      etmat <- cbind(rep(1, length(xdat)), ydat[, etal])
      etainit <- c(init.vals$eta, rep(0, length(etal)))[1:npet]
    }
    if (is.null(taul)) {
      tamat <- as.matrix(rep(1, length(xdat)))
      tauinit <- init.vals$tau
    }else {
      z$trans <- TRUE
      tamat <- cbind(rep(1, length(xdat)), ydat[, taul])
      tauinit <- c(init.vals$tau, rep(0, length(taul)))[1:npta]
    }
    if (is.null(eta2l)) {
      e2mat <- as.matrix(rep(1, length(xdat)))
      eta2init <- init.vals$eta2
    }else {
      z$trans <- TRUE
      e2mat <- cbind(rep(1, length(xdat)), ydat[, eta2l])
      eta2init <- c(init.vals$eta2, rep(0, length(eta2l)))[1:npe2]
    }
    
    init <- c(muinit,siginit,shinit)
    if (!theta_zero) {init <- c(init,thetainit)} # add theta init (in case)
    init <- c(init,etainit)                      # add eta init   (always)
    if (!eta2_zero)  {init <- c(init,eta2init)}  # add eta2 init  (in case)
    if (!tau_zero)   {init <- c(init,tauinit)}   # add tau init   (in case)
     
    # function to calculate neg log-likelihood:
    gev.lik <- function(a) {
      # computes neg log lik of d-gev model
      mu <- mutlink$linkinv(mumat %*% (a[1:npmu]))
      sigma <- sigma0link$linkinv(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
      xi <- xilink$linkinv(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
      # Next line will set the theta likelihood as non-existent in case user requested it. (same for tau)
      if(!theta_zero) {theta <- thetalink$linkinv(thmat %*% (a[seq(npmu + npsc + npsh + 1, length = npth)]))}
      eta <- etalink$linkinv(etmat %*% (a[seq(npmu + npsc + npsh + npth + 1, length = npet)]))
      if(!eta2_zero)  {eta2  <- eta2link$linkinv( e2mat %*% (a[seq(npmu + npsc + npsh + npth + npet + 1, length = npe2)]))}
      if(!tau_zero)   {tau   <- taulink$linkinv(  tamat %*% (a[seq(npmu + npsc + npsh + npth + npet + npe2 + 1, length = npta)]))}
      
      ifelse(!theta_zero, ds.t <- ds+theta, ds.t <- ds) #Don't use theta if user requested not to have it.
      #ifelse(!tau_zero,   sigma.d <- sigma/(ds.t^eta)+tau, sigma.d <- sigma/(ds.t^eta)) # don't use tau if user requested not to have it.
      if (tau_zero){ # don't use tau if user requested not to have it.
        if (eta2_zero){ # don't use eta2
          sigma.d <- sigma/(ds.t^eta)
          mu.d    <- mu*sigma.d
        }else{ # use eta2 (and no tau)
          sigma.d <-    sigma/(ds.t^eta2)
          mu.d    <- mu*sigma/(ds.t^eta)
        }
      }else{ # use tau
        if (eta2_zero){ # don't use eta2 
          sigma.d <- sigma/(ds.t^eta)+tau
          mu.d    <- mu*sigma.d
        }else{ # use eta2 (and tau)
          sigma.d <-     sigma/(ds.t^eta2)+tau
          mu.d    <- mu*(sigma/(ds.t^eta)+tau)
        }
      }
      #sigma.d <- sigma/(ds.t^eta) # old
      
      y = (xdat - mu.d) / sigma.d # new
      #y = (xdat - mu*sigma.d) / sigma.d # derivation
      #y <- xdat/sigma.d - mu # old
      
      y <- 1 + xi * y


      if(!theta_zero) {if(any(theta < 0)) {return(10^6)} } # check definition condition for theta
      if(any(eta <= 0) || any(sigma.d <= 0) || any(y <= 0)) return(10^6)
      if(!tau_zero)   {if(any(tau < 0))    {return(10^6)} } # check definition condition for tau.
      if(!eta2_zero) {if(any(eta2 < 0))    {return(10^6)} } # check definition condition for eta2.
      
      sum(log(sigma.d)) +                 sum(y^(-1/xi)) +              sum(log(y) * (1/xi + 1))
    }


    # finding minimum of log-likelihood:
    x <- optim(init, gev.lik, hessian = TRUE, method = method,
               control = list(maxit = maxit, ...))

    # saving output parameters:
    z$conv <- x$convergence
    mut <- mutlink$linkinv(mumat %*% (x$par[1:npmu]))
    sc0 <- sigma0link$linkinv(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
    xi <- xilink$linkinv(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
    if(!theta_zero){ #When user does NOT set theta parameter to zero (default)
      theta <- thetalink$linkinv(thmat %*% (x$par[seq(npmu + npsc + npsh + 1, length = npth)]))
    }else{ #When user requests theta_parameter to be zero
      theta <- thetalink$linkinv(thmat %*% (0))
    }
    eta <- etalink$linkinv(etmat %*% (x$par[seq(npmu + npsc + npsh + npth + 1, length = npet)]))
    
    if(!eta2_zero){  #When user wants to use eta2 parameter
      eta2 <- eta2link$linkinv(e2mat %*% (x$par[seq(npmu + npsc + npsh + npth + npet + 1,length = npe2)]))
    }else{ #When user requests not to have eta2 parameter (default)
      eta2 <- eta
    }
    
    if(!tau_zero){  #When user does NOT set tau parameter to zero (not default)
      tau <- taulink$linkinv(tamat %*% (x$par[seq(npmu + npsc + npsh + npth + npet + npe2 + 1,length = npta)]))
    }else{ #When user requests tau parameter to be zero
      tau <- taulink$linkinv(tamat %*% (0))
    }
      
    z$nllh <- x$value
    # normalize data to standard Gumbel:
    #sc.d  <-      sc0/((ds+theta)^eta)+tau  # old
    sc.d  <-      sc0/((ds+theta)^eta2)+tau  # new
    mut.d <- mut*(sc0/((ds+theta)^eta )+tau) # new

    #z$data <-  - log(as.vector((1 + xi * (xdat/sc.d-mut))^(-1/xi)))  # old
    z$data <-  - log(as.vector((1 + xi * ((xdat-mut.d)/sc.d))^(-1/xi))) # new
    z$mle <- x$par
    test <- try(              # catch error
    z$cov <- solve(x$hessian) # invert hessian to get estimation on var-covar-matrix
    ,silent = TRUE )
    if("try-error" %in% class(test)){
      warning("Hessian could not be inverted. NAs were produced.")
      z$cov <- matrix(NA,length(z$mle),length(z$mle))
        }
    z$se <- sqrt(diag(z$cov)) # sqrt(digonal entries) = standart error of mle's
    'if (!theta_zero) {#When theta parameter is returned (default)
      if (!tau_zero){ # when tau is returned
        z$vals <- cbind(mut, sc0, xi, theta, eta, tau)
      }else{ # when tau is not returned
        z$vals <- cbind(mut, sc0, xi, theta, eta)
      }
    } else {#When theta parameter is not returned, asked by user
      if (!tau_zero){ # if tau is returned
        z$vals <- cbind(mut, sc0, xi, eta, tau)
      }else{ # if tau is not returned
        z$vals <- cbind(mut, sc0, xi, eta)
      }
    }'
    z$vals <- cbind(mut, sc0, xi, theta, eta, eta2, tau)
    z$init.vals <- unlist(init.vals)

    'names2 = c("mut","sigma0","xi")               # fixed part for set of names
    if(!theta_zero){names2=c(names2,"theta")}     # add theta (in case)
    names2 = c(names2, "eta")                     # add eta (always)
    if(!tau_zero){names2=c(names2, "tau")}        # add tau (in case)
    colnames(z$vals) <- names2'
    colnames(z$vals) <- c("mut","sigma0","xi","theta","eta","eta2","tau")
    
    z$ds <- ds
    z$theta_zero <- theta_zero #Indicates if theta parameter was set to zero by user.
    z$tau_zero <- tau_zero     #Indicates if tau parameter was set to zero by user. (default)
    z$eta2_zero <- eta2_zero   #Indicates if eta2 parameter was set to zero by user. (default)
    if(show) {
      if(z$trans) { # for nonstationary fit
        print(z[c(2, 4)]) # print model, link (3) , conv
        # print names of link functions:
        cat('$link\n')
        #if(!tau_zero){
        print(c(z$link$mutlink$name,z$link$sigma0link$name,z$link$xilink$name,z$link$thetalink$name,z$link$etalink$name,z$link$eta2link$name,z$link$taulink$name))
        #} else{
        #  print(c(z$link$mutlink$name,z$link$sigma0link$name,z$link$xilink$name,z$link$thetalink$name,z$link$etalink$name))
        #}
        cat('\n')
      }else{print(z[4])} # for stationary fit print only conv
      if(!z$conv){ # if fit converged
        print(z[c(5, 7, 9)]) # print nll, mle, se
      }
    }
    class(z) <- "gev.d.fit"
    invisible(z)
}


#### gev.d.init ####

# function to get initial values for gev.d.fit:
# obtain initial values
# by fitting every duration separately

# possible ways to improve:
# take given initial values into account, if there are any
# xi -> mean vs. median ... how do we improve that?
# mu_tilde -> is not very good for small sample sizes yet
# improved initial value for eta, by fitting both mu~d and sigma~d in log-log scale

#' @title get initial values for gev.d.fit
#' @description obtain initial values by fitting every duration separately
#' @param xdat vector of maxima for different durations
#' @param ds vector of durations belonging to maxima in xdat
#' @param link list of 5, link functions for parameters, created with \code{\link{make.link}}
#' @return list of initial values for mu_tilde, sigma_0, xi, theta, eta, eta2, tau
#' @importFrom stats lm
#' @importFrom stats median
#' @importFrom ismev gev.fit
#' @keywords internal

gev.d.init <- function(xdat,ds,link){
  durs <- unique(ds)
  mles <- matrix(NA, nrow=length(durs), ncol= 3)
  for(i in 1:length(durs)){
    test <- try(fit <- ismev::gev.fit(xdat[ds==durs[i]],show = FALSE),silent = TRUE)
    if("try-error" %in% class(test) | fit$conv!=0){mles[i,] <- rep(NA,3)}else{mles[i,] <- fit$mle}
  }
  if(all(is.na(mles))){stop('Initial values could not be computed for this dataset.')}
  # get values for sig0 and eta (also mu_0) from linear model in log-log scale
  lmsig <- lm(log(mles[,2])~log(durs))
  lmmu <- lm(log(mles[,1])~log(durs))

  # sig0 <- exp Intercept
  siginit <- link$sigma0link$linkfun(exp(lmsig$coefficients[[1]]))
  # eta <- mean of negativ slopes
  etainit <- link$etalink$linkfun(mean(c(-lmsig$coefficients[[2]],-lmmu$coefficients[[2]])))
  eta2init <- etainit + 0.0
  # mean of mu_d/sig_d
  # could try:
  # mu0/sig0 = exp(lmmu$coefficients[[1]])/exp(lmsig$coefficients[[1]])
  muinit <- link$mutlink$linkfun(median(c(mles[,1]/mles[,2]),na.rm = TRUE))
  # mean of shape parameters
  shinit <- link$xilink$linkfun(median(mles[,3],na.rm = TRUE))
  thetainit <- link$thetalink$linkfun(0)
  tauinit <- link$taulink$linkfun(0)

  return(list(mu=muinit,sigma=siginit,xi=shinit,theta=thetainit,eta=etainit, eta2=eta2init, tau=tauinit))
}

#### gev.d.lik ####

#' d-GEV Likelihood
#'
#' Computes (log-) likelihood of d-GEV model
#' @param xdat numeric vector containing observations
#' @param ds numeric vector containing corresponding durations (1/60 corresponds to 1 minute, 1 corresponds to 1 hour)
#' @param mut,sigma0,xi,theta,eta,eta2,tau numeric vectors containing corresponding estimates for each of the parameters
#' @param log Logical; if TRUE, the log likelihood is returned.
#'
#' @return single value containing (log) likelihood
#' @export
#'
#' @examples
#' # compute log-likelihood of observation values not included in fit
#' train.set <- example[example$d!=2,]
#' test.set <- example[example$d==2,]
#' fit <- gev.d.fit(train.set$dat,train.set$d,mutl = c(1,2),sigma0l = 1
#'           ,ydat = as.matrix(train.set[c('cov1','cov2')]))
#' params <- gev.d.params(fit,ydat = as.matrix(test.set[c('cov1','cov2')]))
#' gev.d.lik(xdat = test.set$dat,ds = test.set$d,mut = params[,1],sigma0 = params[,2],xi = params[,3]
#'           ,theta = params[,4],eta = params[,5],log=TRUE)
gev.d.lik <- function(xdat,ds,mut,sigma0,xi,theta,eta,log=FALSE,tau=0,eta2=NULL) {
  if (is.null(eta2)){eta2=eta}
  if(any(xi==0)){stop('Function is not defined for shape parameter of zero.')}
  if(any(! c(length(ds),length(mut),length(sigma0),length(xi),length(theta),length(eta),length(eta2),length(tau)) %in%
         c(1,length(xdat)))){
    stop('Input vectors differ in length, but must have the same length.')
  }

  ds.t <- ds+theta
  sigma.d <-     sigma0/(ds.t^eta2)+tau
  mu.d    <- mut*(sigma0/(ds.t^eta)+tau)
  y = (xdat - mu.d) / sigma.d # new
  y <- 1 + xi * y
  
  #sigma.d <- sigma0/(ds.t^eta) + tau # old
  #y <- xdat/sigma.d - mut            # old
  #y <- 1 + xi * y                    # old
  
  if(log){
    return(-sum(log(sigma.d) + y^(-1/xi) + log(y) * (1/xi + 1)))
  }else{
    return(-prod(sigma.d * exp(y^(-1/xi)) * y ^ (1/xi + 1)))
  }

}

#### gev.d.diag ####

#' Diagnostic Plots for d-gev Models
#'
#' @description  Produces diagnostic plots for d-gev models using
#' the output of the function \code{\link{gev.d.fit}}. Values for different durations can be plotted in
#' different colors of with different symbols.
#' @param fit object returned by \code{\link{gev.d.fit}}
#' @param subset an optional vector specifying a subset of observations to be used in the plot
#' @param cols optional either one value or vector of same length as \code{unique(fit$ds)} to
#' specify the colors of plotting points.
#' The default uses the \code{rainbow} function.
#' @param pch optional either one value or vector of same length as \code{unique(fit$ds)} containing
#' integers or symbols to specify the plotting points.
#' @param which string containing 'both', 'pp' or 'qq' to specify, which plots should be produced.
#' @param mfrow vector specifying layout of plots. If both plots should be produced separately,
#' set to \code{c(1,1)}.
#' @param legend logical indicating if legends should be plotted
#' @param title character vector of length 2, giving the titles for the pp- and the qq-plot
#' @param emp.lab,mod.lab character string containing names for empirical and model axis
#' @param ... additional parameters passed on to the plotting function
#'
#' @export
#' @importFrom graphics plot abline par title
#' @importFrom grDevices rainbow
#'
#' @examples
#' data('example',package ='IDF')
#'
#' fit <- gev.d.fit(xdat=example$dat,ds = example$d,ydat=as.matrix(example[,c('cov1','cov2')])
#'                  ,mutl=c(1,2),sigma0l=1)
#' # diagnostic plots for complete data
#' gev.d.diag(fit,pch=1)
#' # diagnostic plots for subset of data (e.g. one station)
#' gev.d.diag(fit,subset = example$cov1==1,pch=1)
gev.d.diag <- function(fit,subset=NULL,cols=NULL,pch=NULL,which='both',mfrow=c(1,2),legend=TRUE,
                       title=c('Residual Probability Plot','Residual Quantile Plot'),
                       emp.lab='Empirical',mod.lab='Model',...){
  # check parameter:
  if(!is.element(which,c('both','pp','qq'))) stop("Parameter 'which'= ",which,
                                                 " but only 'both','pp' or 'qq' are allowed.")
  # subset data
  df <- data.frame(data=fit$data,ds=fit$ds)
  if(!is.null(subset)){
    if(dim(df)[1]!=length(subset)){stop("Length of 'subset' does not match length of data
                                        'xdat' used for fitting.")}
    df <- df[subset,]
  }
  # get single durations
  durs <- sort(unique(df$ds))
  # rescale durations to assign colors
  df$cval <- sapply(df$ds,function(d){which(durs==d)})

  # sort data
  df <- df[order(df$data),]

  # plotting position
  n <- length(df$data)
  px <- (1:n)/(n + 1)

  # create plots:
  if(which=='both') par(mfrow=mfrow) # 2 subplots
  # colors and symbols
  if(is.null(cols))cols <- rainbow(length(durs))
  if(is.null(pch))pch <- df$cval

  if(which=='both'|which=='pp'){
    # pp
    plot(px, exp( - exp( - df$data)), xlab =
           emp.lab, ylab = mod.lab,col=cols[df$cval],pch=pch,...)
    abline(0, 1, col = 1,lwd=1)
    title(title[1])
    if(legend){legend('bottomright',legend = round(durs,digits = 2),pch=pch,
                      col = cols[1:length(durs)],title = 'Duration [h]',ncol = 2)}
  }
  if(which=='both'|which=='qq'){
    # qq
    plot( - log( - log(px)), df$data, ylab =
            emp.lab, xlab = mod.lab,col=cols[df$cval],pch=pch,...)
    abline(0, 1, col = 1,lwd=1)
    title(title[2])
    if(legend){legend('bottomright',legend = round(durs,digits = 2),pch=pch,
                      col = cols[1:length(durs)],title = 'Duration [h]',ncol = 2)}
  }
  if(which=='both') par(mfrow=c(1,1)) # reset par
}

#### gev.d.params ####

#' Calculate gev(d) parameters from \code{gev.d.fit} output
#'
#' @description function to calculate mut, sigma0, xi, theta, eta, eta2, tau
#' (modified location, scale offset, shape, duration offset, duration exponent, second duration exponent, intensity offset)
#' from results of \code{\link{gev.d.fit}} with covariates or link functions other than identity.
#' @param fit fit object returned by \code{\link{gev.d.fit}} or \code{\link{gev.fit}}
#' @param ydat A matrix containing the covariates in the same order as used in \code{gev.d.fit}.
#' @seealso \code{\link{IDF-package}}
#' @return data.frame containing mu_tilde, sigma0, xi, theta, eta, eta2, tau (or mu, sigma, xi for gev.fit objects)
#' @export
#'
#' @examples
#' data('example',package = 'IDF')
#' fit <- gev.d.fit(example$dat,example$d,ydat = as.matrix(example[,c("cov1","cov2")])
#'                  ,mutl = c(1,2),sigma0l = 1)
#' gev.d.params(fit = fit,ydat = cbind(c(0.9,1),c(0.5,1)))


gev.d.params <- function(fit,ydat=NULL){
  if(!class(fit)%in%c("gev.d.fit","gev.fit"))stop("'fit' must be an object returned by 'gev.d.fit' or 'gev.fit'.")
  if(!is.null(ydat)){
    # check covariates matrix
    if(!is.matrix(ydat))stop("'ydat' must be of class matrix.")
    n.par <- max(sapply(fit$model,function(x){return(ifelse(is.null(x),0,max(x)))}))
    if(n.par>ncol(ydat))stop("Covariates-Matrix 'ydat' has ",ncol(ydat), " columns, but ", n.par," are required.")
  }else{if(!fit$trans){# no model -> no covariates matrix
    ydat <- matrix(1)
    }else{stop("To calculate parameter estimates, covariates matrix 'ydat' must be provided.")}
  }

  # number of parameters
  npmu <- length(fit$model[[1]]) + 1
  npsc <- length(fit$model[[2]]) + 1
  npsh <- length(fit$model[[3]]) + 1
  if(class(fit)=="gev.d.fit"){
    if(!fit$theta_zero){
      npth <- length(fit$model[[4]]) + 1 #Including theta parameter (default)]
    }else{
      npth <- 0
    }#With no theta parameter, asked by user
    npet <- length(fit$model[[5]]) + 1
    if(!fit$eta2_zero){
      npe2 <- length(fit$model[[6]]) + 1 #Including eta2 parameter (not default)]
    }else{
      npe2 <- 0
    }
    if(!fit$tau_zero){
      npta <- length(fit$model[[7]]) + 1 #Including tau parameter (not default)]
    }else{
      npta <- 0
    }#With no tau parameter, asked by user
  }

  # inverse link functions
  if(class(fit)=="gev.d.fit"){
    mulink <- fit$link$mutlink$linkinv
    siglink <- fit$link$sigma0link$linkinv
    shlink <- fit$link$xilink$linkinv
    if(!fit$theta_zero) thetalink <- fit$link$thetalink$linkinv
    etalink <- fit$link$etalink$linkinv
    if(!fit$eta2_zero) eta2link <- fit$link$eta2link$linkinv
    if(!fit$tau_zero) taulink <- fit$link$taulink$linkinv
  }else{
    mulink <- eval(parse(text=fit$link))[[1]]
    siglink <- eval(parse(text=fit$link))[[2]]
    shlink <- eval(parse(text=fit$link))[[3]]
  }

  # covariates matrices
  mumat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[1]]],dim(ydat)[1],npmu-1))
  sigmat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[2]]],dim(ydat)[1],npsc-1))
  shmat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[3]]],dim(ydat)[1],npsh-1))
  if(class(fit)=="gev.d.fit"){
    if(!fit$theta_zero){thmat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[4]]],dim(ydat)[1],npth-1))}
    etmat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[5]]],dim(ydat)[1],npet-1))
    if(!fit$eta2_zero) {e2mat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[6]]],dim(ydat)[1],npe2-1))}
    if(!fit$tau_zero)  {tamat <- cbind(rep(1, dim(ydat)[1]), matrix(ydat[, fit$model[[7]]],dim(ydat)[1],npta-1))}
  }

  # calculate parameters
  mut <- mulink(mumat %*% (fit$mle[1:npmu]))
  sc0 <- siglink(sigmat %*% (fit$mle[seq(npmu + 1, length = npsc)]))
  xi <- shlink(shmat %*% (fit$mle[seq(npmu + npsc + 1, length = npsh)]))
  if(class(fit)=="gev.d.fit"){
    if(!fit$theta_zero){
      theta <- thetalink(thmat %*% (fit$mle[seq(npmu + npsc + npsh + 1, length = npth)]))
    }else{
      theta <- rep(0,dim(ydat)[1])
    }
    eta <- etalink(etmat %*% (fit$mle[seq(npmu + npsc + npsh + npth + 1, length = npet)]))
    if(!fit$eta2_zero){
      eta2 <- eta2link(e2mat %*% (fit$mle[seq(npmu + npsc + npsh + npth + npet + 1, length = npe2)]))
    }else{
      eta2 <- eta #rep(0,dim(ydat)[1])
    }
    if(!fit$tau_zero){
      tau <- taulink(tamat %*% (fit$mle[seq(npmu + npsc + npsh + npth + npet + npe2 + 1, length = npta)]))
    }else{
      tau <- rep(0,dim(ydat)[1])
    }
    return(data.frame(mut=mut,sigma0=sc0,xi=xi,theta=theta,eta=eta, eta2=eta2, tau=tau))
  }else{
    return(data.frame(mu=mut,sig=sc0,xi=xi))
  }
}


#### example data ####
#' Sampled data for duration-dependent GEV
#' 
#' @description
#' Randomly sampled data set used for running the example code, containing:
#' \itemize{
#'   \item \code{$xdat}: 'annual' maxima values
#'   \item \code{$ds}: corresponding durations
#'   \item \code{$cov1}, \code{$cov2}: covariates}
#' d-GEV parameters used for sampling:
#' \itemize{
#'   \item \eqn{\tilde{\mu} = 4 + 0.2 cov_1 +0.5 cov_2}
#'   \item \eqn{\sigma_0 = 2+0.5 cov_1}
#'   \item \eqn{\xi = 0.5}
#'   \item \eqn{\theta = 0}
#'   \item \eqn{\eta = 0.5}
#'   \item \eqn{\eta_2 = 0.5}
#'   \item \eqn{\tau = 0}}
#'
#'
#' @docType data
#' @keywords datasets
#' @name example
#' @usage data('example',package ='IDF')
#' @format A data frame with 330 rows and 4 variables
NULL
