# This file contains:
# -IDF-package description
# -IDF.agg for preparing the data
# -IDF.plot for plotting of IDF curves at a chosen station


#### IDF-package ####

#' Introduction
#' @name IDF-package
#' @docType package
#' @description This package provides functions to estimate IDF relations for given
#' precipitation time series on the basis of a duration-dependent
#' generalized extreme value distribution (d-GEV). 
#' The central function is \code{\link{gev.d.fit}}, which uses the method 
#' of maximum-likelihood estimation for the d-GEV parameters, whereby it is 
#' possible to include generalized linear modeling 
#' for each parameter. This function was implemented on the basis of \code{\link[ismev]{gev.fit}}.
#' For more detailed information on the methods and the application of the package for estimating 
#' IDF curves with spatial covariates, see Ulrich et. al (2020). 
#' @details 
#' * The __d-GEV__ is defined following Koutsoyiannis et al. (1998): 
#' \deqn{G(x)= \exp[-( 1+\xi(x/\sigma(d)- \tilde{\mu}) ) ^{-1/\xi}] } 
#' defined on \eqn{ \{ x: 1+\xi(x/\sigma(d)- \tilde{\mu} > 0) \} },
#' with the duration dependent scale parameter \eqn{\sigma(d)=\sigma_0/(d+\theta)^\eta > 0},
#' modified location parameter \eqn{\tilde{\mu}=\mu/\sigma(d)\in R}
#' and shape parameter \eqn{\xi\in R}, \eqn{\xi\neq 0}.
#' The parameters \eqn{\theta \leq 0} and \eqn{0<\eta<1} are duration offset and duration exponent
#' and describe the slope and curvature in the resulting IDF curves, respectively.
#' * A useful introduction to __Maximum Likelihood Estimation__ for fitting for the 
#' generalized extreme value distribution (GEV) is provided by Coles (2001). It should be noted, however, that this method uses
#' the assumption that block maxima (of different durations or stations) are independent of each other. 
#' @references 
#' * Ulrich, J.; Jurado, O.E.; Peter, M.; Scheibel, M.; 
#' Rust, H.W. Estimating IDF Curves Consistently over Durations with Spatial Covariates. Water 2020, 12, 3119,
#' https://doi.org/10.3390/w12113119
#' * Demetris Koutsoyiannis, Demosthenes Kozonis, Alexandros Manetas,
#' A mathematical framework for studying rainfall intensity-duration-frequency relationships,
#' Journal of Hydrology,
#' Volume 206, Issues 1–2,1998,Pages 118-135,ISSN 0022-1694, https://doi.org/10.1016/S0022-1694(98)00097-3
#' * Coles, S.An Introduction to Statistical Modeling of Extreme Values;   Springer:  New York, NY, USA, 2001,
#' https://doi.org/10.1198/tech.2002.s73
#' @md
#' 
#' @examples 
#' ## Here are a few examples to illustrate the order in which the functions are intended to be used.
#' 
#' ## Step 0: sample 20 years of example hourly 'precipitation' data
# dates <- seq(as.POSIXct("2000-01-01 00:00:00"),as.POSIXct("2019-12-31 23:00:00"),by = 'hour')
# sample.precip <- rgamma(n = length(dates), shape = 0.05, rate = 0.4)
# precip.df <- data.frame(date=dates,RR=sample.precip)
# 
# ## Step 1: get annual maxima
# durations <- 2^(0:6) # accumulation durations [h] 
# ann.max <- IDF.agg(list(precip.df),ds=durations,na.accept = 0.1)
# # plotting the annual maxima in log-log representation
# plot(ann.max$ds,ann.max$xdat,log='xy',xlab = 'Duration [h]',ylab='Intensity [mm/h]')
# 
# ## Step 2: fit d-GEV to annual maxima
# fit <- gev.d.fit(xdat = ann.max$xdat,ds = ann.max$ds,sigma0link = make.link('log'))
# # checking the fit 
# gev.d.diag(fit,pch=1,legend = FALSE)
# # parameter estimates 
# params <- gev.d.params(fit)
# print(params)
# # plotting the probability density for a single duration 
# q.min <- floor(min(ann.max$xdat[ann.max$ds%in%1:2]))
# q.max <- ceiling(max(ann.max$xdat[ann.max$ds%in%1:2]))
# q <- seq(q.min,q.max,0.2)
# plot(range(q),c(0,0.55),type = 'n',xlab = 'Intensity [mm/h]',ylab = 'Density')
# for(d in 1:2){ # d=1h and d=2h
#   hist(ann.max$xdat[ann.max$ds==d],main = paste('d=',d),q.min:q.max
#        ,freq = FALSE,add=TRUE,border = d)   # sampled data
#   lines(q,dgev.d(q,params$mut,params$sigma0,params$xi,params$theta,params$eta,d = d),col=d) # etimated prob. density
# }
# legend('topright',col=1:2,lwd=1,legend = paste('d=',1:2,'h'),title = 'Duration')
# 
# ## Step 3: adding the IDF-curves to the data
# plot(ann.max$ds,ann.max$xdat,log='xy',xlab = 'Duration [h]',ylab='Intensity [mm/h]')
# IDF.plot(durations,params,add=TRUE)
NULL

#### IDF.agg ####

#' Aggregation and annual maxima for chosen durations
#' @description Aggregates several time series for chosen durations and finds annual maxima 
#' (either for the whole year or chosen months). Returns data.frame that can be used for
#' the function \code{\link{gev.d.fit}}.
#'
#' @param data list of data.frames containing time series for every station. 
#' The data.frame must have the columns 'date' and 'RR' unless other names 
#' are specified in the parameter `names`. The column 'date' must contain strings with 
#' standard date format.
#' @param ds numeric vector of aggregation durations in hours. 
#' (Must be multiples of time resolution at all stations.)
#' @param na.accept numeric in [0,1) giving maximum percentage of missing values 
#' for which block max. should still be calculated.
#' @param which.stations optional, subset of stations. Either numeric vector or character vector 
#' containing names of elements in data. If not given, all elements in `data` will be used.
#' @param which.mon optional, subset of months (as list containing values from 0 to 11) of which to calculate the annual maxima from. 
#' @param names optional, character vector of length 2, containing the names of the columns to be used. 
#' @param cl optional, number of cores to be used from \code{\link[parallel]{parLapply}} for parallel computing.
#'
#' @details If data contains stations with different time resolutions that need to be aggregated at
#' different durations, IDF.agg needs to be run separately for the different groups of stations. 
#' Afterwards the results can be joint together using `rbind`.
#'
#' @return data.frame containing the annual intensity maxima [mm/h] in `$xdat`, the corresponding duration in `$ds`,
#' the `$year` and month (`$mon`) in which the maxima occurred 
#' and the station id or name in `$station`.
#' 
#' @seealso \code{\link{pgev.d}}
#' 
#' @export
#' @importFrom parallel parLapply makeCluster stopCluster
#' @importFrom pbapply pblapply 
#' @importFrom RcppRoll roll_sum
#' @importFrom fastmatch ctapply
#'
#' @examples
#' dates <- as.Date("2019-01-01")+0:729
#' x <- rgamma(n = 730, shape = 0.4, rate = 0.5)
#' df <- data.frame(date=dates,RR=x)
#' 
#' # get annual maxima
#' IDF.agg(list('Sample'= df),ds=c(24,48),na.accept = 0.01)
#' 
#' ##      xdat ds year  mon station
#' ## 0.2853811 24 2019 0:11  Sample
#' ## 0.5673122 24 2020 0:11  Sample
#' ## 0.1598448 48 2019 0:11  Sample
#' ## 0.3112713 48 2020 0:11  Sample
#' 
#' # get monthly maxima for each month of june, july and august
#' IDF.agg(list('Sample'=df),ds=c(24,48),na.accept = 0.01,which.mon = list(5,6,7))
#' 
#' # get maxima for time range from june to august
#' IDF.agg(list('Sample'=df),ds=c(24,48),na.accept = 0.01,which.mon = list(5:7))
#' 
    IDF.agg <- function(data,ds,na.accept = 0,
                        which.stations = NULL,which.mon = list(0:11),names = c('date','RR'),cl = 1){
      
      if(!inherits(data, "list"))stop("Argument 'data' must be a list, instead it is a: ", class(data))
      
      # function 2: aggregate station data over durations and find annual maxima:                                
      agg.station <- function(station){
        data.s <- data[[station]]
        if(!is.data.frame(data.s)){stop("Elements of 'data' must be data.frames. But element "
                                        ,station," contains: ", class(data.s))}
        if(sum(is.element(names[1:2],names(data.s)))!=2){stop('Dataframe of station ', station 
                                                              ,' does not contain $', names[1]
                                                              ,' or $', names[2], '.')}
        dtime<-as.numeric((data.s[,names[1]][2]-data.s[,names[1]][1]),units="hours")
        
        if(any((ds/dtime)%%1 > 10e-8)){
          stop('At least one of the given aggregation durations is not multiple of the time resolution = '
               ,dtime,'hours at station ',station,'.')}
        
        # function 1: aggregate over single durations and find annual maxima:
        agg.ts <- function(ds){
          runsum = RcppRoll::roll_sum(data.s[,names[2]],round(ds/dtime),fill=NA,align='right')
          #runmean <- rollapplyr(as.zoo(data.s[,names[2]]),ds/dtime,FUN=sum,fill =NA,align='right')
          runsum <- runsum/ds #intensity per hour
          max.subset <- lapply(1:length(which.mon),function(m.i){
            subset <- is.element(as.POSIXlt(data.s[,names[1]])$mon,which.mon[[m.i]])
            max <- fastmatch::ctapply(runsum[subset],(as.POSIXlt(data.s[,names[1]])$year+1900)[subset],
                          function(vec){
                            n.na <- sum(is.na(vec))
                            max <- ifelse(n.na <= na.accept*length(vec),max(vec,na.rm = TRUE),NA)
                            return(max)})
            df <- data.frame(xdat=max,ds=ds,year=as.numeric(names(max)),mon=deparse(which.mon[[m.i]]),
                             station= station,
                             stringsAsFactors = FALSE)
            return(df)})
          df <- do.call(rbind,max.subset)  
          return(df) # maxima for single durations
        }
        # call function 1 in lapply to aggregate over all durations at single station
        clust <- parallel::makeCluster(cl)
        data.agg <- parallel::parLapply(cl = clust,ds,agg.ts)  
        parallel::stopCluster(clust)
        df <- do.call(rbind,data.agg)
        return(df) # maxima for all durations at one station
      }
      # which stations should be used?
      if(is.null(which.stations))which.stations <- if(is.null(names(data))){1:length(data)}else{names(data)}
      # call function 2 in lapply to aggregate over all durations at all stations
      station.list <- pbapply::pblapply(which.stations,agg.station)
      
      return(do.call('rbind',station.list))
    }

#### IDF.plot ####

#' Plotting of IDF curves at a chosen station
#'
#' @param durations vector of durations for which to calculate the quantiles. 
#' @param fitparams vector containing parameters mut, sigma0, xi, theta, eta
#' (modified location, scale offset, shape, duration offset, duration exponent) for chosen station
#' as obtained from \code{\link{gev.d.fit}}
#' (or \code{\link{gev.d.params}} for model with covariates).
#' @param probs vector of non-exceedance probabilities for which to plot IDF curves (p = 1-1/(Return Period))
#' @param cols vector of colors for IDF curves. Should have same length as \code{probs}
#' @param add logical indicating if plot should be added to existing plot, default is FALSE
#' @param legend logical indicating if legend should be plotted (TRUE, the default)
#' @param ... additional parameters passed on to the \code{plot} function
#'
#' @export
#' @importFrom grDevices rgb
#' @importFrom graphics axis box lines plot points 
#' @examples
#' data('example',package = 'IDF')
#' # fit d-gev
#' fit <- gev.d.fit(example$dat,example$d,ydat = as.matrix(example[,c("cov1","cov2")])
#'                  ,mutl = c(1,2),sigma0l = 1)
#' # get parameters for cov1 = 1, cov2 = 1
#' par <- gev.d.params(fit = fit, ydat = matrix(1,1,2))
#' # plot quantiles
#' IDF.plot(durations = seq(0.5,35,0.2),fitparams = par)
#' # add data points
#' points(example[example$cov1==1,]$d,example[example$cov1==1,]$dat)
IDF.plot <- function(durations,fitparams,probs=c(0.5,0.9,0.99),
                     cols=4:2,add=FALSE,
                     legend=TRUE,...){
  
  # if cols is to short, make longer    
  if(length(cols)!=length(probs))cols <- rep_len(cols,length.out=length(probs))
  
  ## calculate IDF values for given probability and durations
  qs <- lapply(durations,qgev.d,p=probs,mut=fitparams[1],sigma0=fitparams[2],xi=fitparams[3],
         theta=fitparams[4],eta=fitparams[5])
  idf.array <- matrix(unlist(qs),length(probs),length(durations)) # array[probs,durs]
  if(!add){ #new plot
    ## initialize plot window
    # check if limits were passed
    if(is.element('ylim',names(list(...)))){
      ylim <- list(...)[['ylim']]
    }else{ylim <- range(idf.array,na.rm=TRUE)}
    if(is.element('xlim',names(list(...)))){
      xlim <- list(...)[['xlim']]
    }else{xlim <- range(durations,na.rm=TRUE)}
    if(is.element('main',names(list(...)))){
      main <- list(...)[['main']]
    }else{main <- ''}
    
    # empty plot
    plot(NA,xlim=xlim,ylim=ylim,xlab="Duration [h]",ylab="Intensity [mm/h]",log="xy",main=main)
  }
  
  ## plot IDF curves
  for(i in 1:length(probs)){
    lines(durations,idf.array[i,],col=cols[i],...)
  }
  
  if(legend){## plot legend
    # check if lwd, lty were passed
    if(is.element('lwd',names(list(...)))){
      lwd <- list(...)[['lwd']]
    }else{lwd <- 1}
    if(is.element('lty',names(list(...)))){
      lty <- list(...)[['lty']]
    }else{lty <- 1}
    
    legend(x="topright",title = 'p-quantile',legend=probs,
           col=cols,lty=lty,lwd=lwd)
  }
}
