##################################################
## IDF package 
## Authors: Sarah Joedicke, Carola Detring, Christoph Ritschel
## Update: 15.09.2017  
###################################################

###############################################
############# Read Data function ##############
###############################################

#' @title Reading precipitation data 
#' @description The function \code{IDF.read} reads a file in table format and creates a \code{data.frame} from it
#' and adds some attributes (station information, aggregation time, data source). The only data values used are: 
#' date, precipitation
#' The \code{data.frame} will have the following format:
#' | year | mon | day | hour | min | RR |
#' |------+-----+-----+------+-----+----+
#' |      |     |     |      |     |    |
#' @usage IDF.read(file, type) 
#' @param file a \code{character string} naming the file from which the data is to be read. 
#' @param type a \code{character string} defining the type of data to be read: either "stadtmessnetz" or "webwerdis", depending on if the data comes from the Stadtmessnetz Berlin
#' or WebWerdis. If type = "webwerdis", the data will be read, then sorted, formatted and missing lines added, 
#' while if type = "stadtmessnetz", the data will just be read and formatted. 
#' Both source types have a different layout in the original file.
#' @return Liste a \code{data.frame} of date and time information and precipitation values for each time step
#' @details This function is designed to prepare a data file for doing an estimation on IDF parameters in function \code{IDF.fit}.
#' The time given in the data is the end time, so the precipitation was measured up to that time.  
#' @seealso read.table, IDF.fit
#' @author Sarah Joedicke \email{sarah.joedicke@@fu-berlin.de}
#' @author Christoph Ritschel \email{christoph.ritschel@@met.fu-berlin.de}
IDF.read <- function(file,type){
  
  if(type != "stadtmessnetz" && type != "webwerdis") {
    
    cat("Warning: wrong type declared for input file")
    stop()
  }
  
  if (type == "stadtmessnetz") {
    
    Tab_MN <- read.csv2(file)  #STADTMESSNETZ
    new_time <- strptime(Tab_MN$Zeitstempel,format="%d.%m.%Y %H:%M")   #STADTMESSNETZ date vector
  }
  
  # Da die Stadtmessnetzdaten (bisher) konstistent aussehen, wird auf das Erstellen einer neuen Tabelle mit sicher allen
  # Zeiten verzichtet, da die Minutendaten sehr gross sind. Sollte es inkonsistente Tabellen geben, sollte man diese seperat behandeln,
  # sonst wird viel Rechenzeit fuer die kompletten Tabellen verschwendet. 
  
  if (type == "webwerdis") {
    Tab <- read.table(file,header=TRUE,sep=";")   #WEBWERDIS
    Tab_kurz <- Tab[,c("Date","precipitation")]
    
    ## Sort table in output format
    time <- strptime(Tab_kurz$Date,format="%Y-%m-%d T %H:%M:%S")
    Tab_sort <- Tab_kurz[order(as.character(time)),]
    time_sort <- strptime(Tab_sort$Date,format="%Y-%m-%d T %H:%M:%S")
    Tab_sort$Date <- as.character(time_sort)
    
    ## If dates are missing, add lines containing NA preicipitation measurments for these time steps. 
    h_diff <- as.numeric(difftime(format(time_sort[length(time_sort)],"%Y-%m-%d T %H:%M:%S") , 
                                  format(time_sort[1],"%Y-%m-%d T %H:%M:%S"),units="hours")) #h_diff is the difference in time steps
    new_time <- seq(time_sort[1], length = h_diff+1, by = "hour")
    new_tab <- data.frame(Date=as.character(new_time), precipitation=NA)  # predefine table with NAs and every time steps
    
    Tab_na <- (merge(Tab_sort, new_tab, "Date", all.y=TRUE))[,1:2]
  }
  
  new_timect <- as.POSIXct(new_time)
  
  J <- as.numeric(format(new_timect,'%Y'))
  M <- as.numeric(format(new_timect,'%m'))
  d <- as.numeric(format(new_timect,'%d'))
  h <- as.numeric(format(new_timect,'%H'))
  m <- as.numeric(format(new_timect,'%M'))
  
  if (type == "webwerdis") Tab_end <- data.frame(J,M,d,h,m,Tab_na$precipitation.x) #WEBWERDIS
  if (type == "stadtmessnetz") Tab_end <- data.frame(J,M,d,h,m,Tab_MN[,2]) #STADTMESSNETZ
  
  ## Name table attributes: 
  
  colnames(Tab_end) <- c("year","mon","day","hour","min","RR")
  attr(Tab_end,"accumulation time (min)") <- as.numeric(difftime(new_timect[2],new_timect[1], units="mins"))
 # Liste <- list(t1=Tab_end)
  Liste <- Tab_end 
 
  if (type == "webwerdis"){
    # WEBWERDIS:
    attr(Liste,"StationName") <- as.character(Tab$Stationname[1])
    attr(Liste,"StationID") <- "NA"
    attr(Liste,"Long (deg N)")  <- Tab$Longitude[1]
    attr(Liste,"Lat (deg E)") <- Tab$Latitude[1]
    attr(Liste,"Heigth (m)")   <- Tab$StationHeight[1]
    attr(Liste,"Source") <- "Web-WERDIS"
  } #Listen-Attribute benennen
  
  if (type == "stadtmessnetz"){
    # STADTMESSNETZ:
    attr(Liste,"StationName") <- colnames(Tab_MN)[2]
    attr(Liste,"StationID") <- "NA"
    attr(Liste,"Long (deg N)")  <- "NA"
    attr(Liste,"Lat (deg E)") <- "NA"
    attr(Liste,"Height (m)")   <- "NA"
    attr(Liste,"Source") <- "Stadtmessnetz"
  } #Listen-Attribute benennen
  
  cat(paste("read.data of", file , "done \n"))
  str(Liste)   # optional; so sieht man beim Einlesen, womit man es zu tun hat und ob alles geklappt hat
  
  return(Liste)
} 
# End of function IDF.read
####################################################################################################################

##### Aggregation ###

#' \code{TS.acc} accumulates a given time series \code{x} at a given accumulation level \code{acc.val}. Minimum value
#' for acc.val is 2 [unit time]. Option for using moving sum is given.
#' @title Accumulation of a time series
#' @param x \code{vector} of a time series
#' @param acc.val \code{value} specifying the accumulation level, minimum value is 2
#' @param moving.sum \code{logical} 'TRUE' means moving sum will be applied
#' @return x.acc \code{TS.acc} returns a \code{vector} of an accumulated time series 
#' @usage TS.acc(x,acc.val,moving.sum="FALSE")
#' @examples
#' TS <- rgamma(n=1000,shape=1)
#' acc.2 <- TS.acc(TS,acc.val=2)
#' \donttest{
#' acc.24 <- TS.acc(TS,acc.val=24,moving.sum=TRUE)
#' }
#' @author Christoph Ritschel \email{christoph.ritschel@@met.fu-berlin.de}
#' @author Carola Detring \email{carola.detring@@met.fu-berlin.de}
TS.acc <- function(x,acc.val=2,moving.sum="FALSE") {
  
  ## check for input value of acc.val
  if(acc.val<1) cat(paste("Warning: accumulation value acc.val too small for accumulation of the time series \n"))
  
  if(moving.sum){
   
    x.acc <- as.numeric(filter(x,filter=rep(1,acc.val),method="convolution",sides=1))
     
  }else{
  
  l.new <- length(x)%/%acc.val ## calculate new length of accumulated time series
  l.rest <- length(x)%%acc.val ## calculate values left over
  if(l.rest==0) {
    x.acc <- apply(matrix(x,nrow=l.new,byrow=T),1,sum) 
  }else{
    x.acc <- apply(matrix(x[1:(length(x)-l.rest)],nrow=l.new,byrow=T),1,sum)   
    #cat(paste("Warning: ",l.rest,"time steps left and not used for accumulation \n"))
  }
  
  }
  
  ## return accumulated time series
  return(x.acc)

} # End of function TS.acc
#####################################################################################


#######################
## Fitting Functions ##
#######################

#'@title Density function of modified generalized extreme value distribution
#'@description The function \code{dgev.d} is a modified version of the function \code{dgev} for different durations \code{d} developed by Koutsoyiannis et al. (1998).
#'@param q Vector of quantiles
#'@param mu location value
#'@param sigma scale value
#'@param xi shape value
#'@param theta value defining the curvature of the IDF
#'@param eta value defining the slope of the IDF
#'@param d vector of durations
#'@param log \code{logical} option to use logarithmic parameter values, default=FALSE
#'@seealso \code{\link[evd]{dgev}}
#'@return dgev.d gives the density function
#' @author Christoph Ritschel \email{christoph.ritschel@@met.fu-berlin.de}

dgev.d <- function(q,mu=0,sigma=1,xi=0,theta=0,eta=1,d=1,log=FALSE) {
  sigma.d <- sigma/(d+theta)^eta
  ##problem if sigma.d is NaN (d+theta) negative and eta smaller than 1 --> cant calculate root of negative value 
  sigma.d[which(is.nan(sigma.d))] <- Inf
  dens <- dgev(q,loc=mu*sigma.d,scale=sigma.d,shape=xi,log=log)
  dens[which(is.nan(dens))] <- NA
  return(dens)
}


#'@title Quantile function of modified generalized extreme value distribution
#'@description The function \code{qgev.d} is a modified version of the function \code{qgev} for different durations \code{d} developed by Koutsoyiannis et al. (1998).
#'@param p Vector of probabilities
#'@param mu location value
#'@param sigma scale value
#'@param xi shape value
#'@param theta value defining the curvature of the IDF
#'@param eta value defining the slope of the IDF
#'@param d vector of durations
#'@param lower.tail \code{logical} if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#'@seealso \code{\link[evd]{qgev}}
#'@return qgev.d gives the quantile function
#' @author Christoph Ritschel \email{christoph.ritschel@@met.fu-berlin.de}
qgev.d <- function(p,mu=0,sigma=1,xi=0,theta=0,eta=1,d=1,lower.tail=TRUE) {
  
  sigma.d <- sigma/(d+theta)^eta
  ##problem if sigma.d is NaN (d+theta) negative and eta smaller than 1 --> cant calculate root of negative value 
  sigma.d[which(is.nan(sigma.d))] <- Inf
  quant <- qgev(p,loc=mu*sigma.d,scale=sigma.d,shape=xi,lower.tail=lower.tail)
  quant[is.infinite(quant)] <- NA
  return(quant)
}

#'@title Random generation for the modified generalized extreme value distribution
#'@description The function \code{rgev.d} is a modified version of the function \code{rgev} for different durations \code{d} developed by Koutsoyiannis et al. (1998).
#'@param n Number of observations
#'@param mu location value
#'@param sigma scale value
#'@param xi shape value
#'@param theta value defining the curvature of the IDF
#'@param eta value defining the slope of the IDF
#'@param d vector of durations
#'@seealso \code{\link[evd]{rgev}}
#'@return rgev.d generates random derivates
#' @author Christoph Ritschel \email{christoph.ritschel@@met.fu-berlin.de}
rgev.d <- function(n,mu=0,sigma=1,xi=0,theta=0,eta=1,d=1) {
  ## gumbel
  sigma.d <- sigma/(d+theta)^eta
  ##problem if sigma.d is NaN (d+theta) negative and eta smaller than 1 --> cant calculate root of negative value 
  sigma.d[which(is.nan(sigma.d))] <- Inf
  x <- rgev(n, loc=mu*sigma.d,scale=sigma.d,shape=xi)
  x[which(is.nan(x))] <- NA
  return(x)
  
  }

#######################################################################
#' @title Negativ log-likelihood of modified GEV
#' @description The function \code{IDF.nll} calculates the negative log-likelihood for a given set of model parameters
#' \code{mu,sigma,xi,theta,eta}, given observations \code{x} and given durations \code{d}. Options for the usage of
#' logartihmic values \code{use.log} and a debugging function \code{DEBUG} are available.
#'@param mu location value
#'@param sigma scale value
#'@param xi shape value
#'@param theta value defining the curvature of the IDF
#'@param eta value defining the slope of the IDF
#'@param x vector of observations at different durations d
#'@param d vector of durations
#'@param use.log \code{logical} value for usage of logarithmic values, default is \code{FALSE}
#'@param DEBUG \code{logical} value for usage of debugging, if \code{TRUE} the input parameters and the value of negative
#'log-likelihood are printed on console.
#'@return retruns weightes negative log-likelihood by number of observatons uesd
#' @author Christoph Ritschel \email{christoph.ritschel@@met.fu-berlin.de}

IDF.nll <- function(mu=0,sigma=1,xi=0,theta=0,eta=1,x,d,use.log=F,DEBUG=F) {
  ## mu is the mu~ from Koutsoyiannis
  
  if(use.log){
    ## ensure that critical parameters are positive
    sigma <- exp(sigma)
    theta <- exp(theta)
    eta <- exp(eta)
  }
  
  sigma.d <- sigma/((d+theta)^eta) 
  if(DEBUG) debug.values <- c(mu,sigma,xi,theta,eta)
  
  if(sum(is.nan(sigma.d))==0) {
  
  ## Weibull und Frechet
  if(xi!=0){
    C <- 1 + xi * (x/sigma.d - mu )
    nll <- switch((sum(C<0,na.rm=T)>0)+1,
                  sum(log(sigma.d),na.rm=T)+(1+1/xi)*sum(log(C),na.rm=T)+sum((C)^(-1/xi),na.rm=T),
                  NA)
    #       + penalty*(sum(C[C<0]^2))
    ## Gumbel
  }else if(xi==0){# & sigma<1 & eta<1) 
    Y <- x/sigma.d-mu
    nll <- -(-sum(log(sigma.d),na.rm=T)-sum((Y),na.rm=T)-sum(exp(-Y),na.rm=T))
  }
  }else{ nll <- NA}
    
  if(DEBUG){ 
    cat(debug.values,nll,"\n")
    options(digits.secs=6)
    ##    debug.values <- c(debug.values,nll,as.character(Sys.time()))
    ##    write(debug.values,file="optim.log",append=TRUE,ncolumns=length(debug.values))
    ##    cat(debug.values,nll,sum(A<0),"\n")
  }
  
  return(nll/length(x))
  
  } # end of function IDF.nll
######################################################################################################

#' @title Fitting function to optimize IDF model parameters
#' @description The function \code{fit.fun} fits IDF model parameters \code{mu,sigma,xi,theta,eta} to a set of given observations \code{obs}, 
#' typically a series of yearly maxima at different durations \code{d}. Options for using logarithmic parameter values and debugging
#' are given. Also the \code{optim} parameters \code{method} and \code{upper,lower} can be defined.
#' @param obs vector of yearly intensity maxima at different durations. Order: Y1D1, Y2D1,...,YnD1,Y1D2,...YnD2,Y1D3,...,YnDk
#' @param dur vector of durations with same length as \code{obs}. Order: n x D1, n x D2, ... n x Dk 
#' @param mu location value
#' @param sigma scale value
#' @param xi shape value
#' @param theta value defining the curvature of the IDF
#' @param eta value defining the slope of the IDF
#' @param use.log \code{logical} value for usage of logarithmic values, default is \code{FALSE}
#' @param DEBUG \code{logical} value for usage of debugging, if \code{TRUE} the input parameters and the value of negative
#' log-likelihood are printed on console for each iteration during optimization.
#' @param method \code{character} defining the method to be used in \code{optim}, preferences are: "Nelder-Mead", "BFGS", "L-BFGS-B"e
#' @param lower \code{vector} specifying the lower boundary of parameters for "L-BFGS-B" method
#' @param upper \code{vector} specifying the upper boundary of parameters for "L-BFGS-B" method
#' @return $min value of negative log-likelihood at optimization minimum
#' @return $par vector of IDF parameters at optimization minimum
#' @author Christoph Ritschel \email{christoph.ritschel@@met.fu-berlin.de}

fit.fun <- function(obs,dur,mu=1,sigma=1,xi=0.5,theta=1,eta=1,use.log=F,DEBUG=F,method="Nelder-Mead",upper=Inf,lower=-Inf) {

  use.log=use.log
  
  if(use.log) {
    if(sigma<=0){sigma <- 1E-10}
    if(theta<=0){theta <- 1E-10}
    if(eta<=0){eta <- 1E-10}
    sigma <- log(sigma)
    theta <- log(theta)
    eta <- log(eta)
    
    if(method=="L-BFGS-B") {
    upper[2] <- log(upper[2])
    upper[4] <- log(upper[4])
    upper[5] <- log(upper[5])
    
    lower[2] <- log(lower[2])
    lower[4] <- log(lower[4])
    lower[5] <- log(lower[5])
    }
    
  }
  
  ## check initial value of negative log-Likelihood function
  nll <- IDF.nll(mu,sigma,xi,theta,eta,x=obs,d=dur,use.log=use.log,DEBUG=DEBUG)
  
  ## if initial value is acceptable...
  if(!is.infinite(nll)&!is.na(nll)) {
    
    
    if(method=="L-BFGS-B") {
      
      ## problem: optimization algrorithm often has difficulities concerning infinite or NA-difference values betweeen iterations
      ## solution: ignore this error message using functon tryCatch and return NULL if there was an error during optimization
      fit <- tryCatch(mle(IDF.nll,start=list(mu=mu,sigma=sigma,xi=xi,theta=theta,eta=eta),fixed=list(x=obs,d=dur,use.log=use.log,DEBUG=DEBUG),
                          control=list(trace=0,maxit=5000),
                          method=method,upper=upper,lower=lower), error=function(e) NULL)#,
      #upper=upper,lower=lower)
      
    }else{
      
      ## problem: optimization algrorithm often has difficulities concerning infinite or NA-difference values betweeen iterations
      ## solution: ignore this error message using functon tryCatch and return NULL if there was an error during optimization
      fit <- tryCatch(mle(IDF.nll,start=list(mu=mu,sigma=sigma,xi=xi,theta=theta,eta=eta),fixed=list(x=obs,d=dur,use.log=use.log,DEBUG=DEBUG),
                          control=list(trace=0,maxit=5000),
                          method=method), error=function(e) NULL)#,
      #upper=upper,lower=lower)
      
      
      
    }
    
    ## if there was no error
    if(!is.null(fit)) {
      fit.min <- fit@min
      fit.par <- fit@coef
    }else { ## else return NA
      fit.min <- NA
      fit.par <- rep(NA,5)  
    } ## end if error
    
  }else { ## else retunr NA
    
    fit.min <- NA
    fit.par <- rep(NA,5)  
    
  } ## end if initial value..
  
  if(use.log){
    fit.par[2] <- exp(fit.par[2])
    fit.par[4] <- exp(fit.par[4])
    fit.par[5] <- exp(fit.par[5])
  }
  names(fit.par) <- c("mu","sigma","xi","theta","eta")
  
  return(list("min"=fit.min,"par"=fit.par))
  
} ## end of function fit.fun
##################################################################################

#' @title Fitting IDF model parameters to observations at different durations
#' @description The function \code{IDF.fit} fits the IDF model parameters \code{mu,sigma,xi,eta,theta}
#' to a data.frame of observations \code{data} with temporal inforamtion (at least years) and values of precipitation
#' at a given temporal resoultion. This precipitation time series gets aggregated at given aggregation levels.
#' \code{agg.lev} and yearly maxima of intensity are caluclated for a specific month or the whole year/dataset. 
#' The starting values of the IDF model parameters can be determined by the user as well as specific options to use
#' during optimization. Logartihmic transformation, debugging, the optimization method, and an option to plot the
#' IDF curves.
#' @param data a \code{data,frame}, preferably generated by function \code{IDF.read}. It should at least contain a \code{$RR} and \code{$year} element for the 
#' function tow work properly.
#' @param agg.lev a vector of aggregation levels used to fit the IDF curves.
#' @param month \code{integer} value specifying the month to be used for estimating the IDF parameters. Type "all" for all months or if
#' the whole time series should be fitted.
#' @param moving.sum \code{logical} specifying if moving sum filtering should be applied for time series aggregation.
#' @param theta.init inital value defining the curvature of the IDF, default is zero, it is not recommended to change it
#' @param use.log \code{logical} value for usage of logarithmic values, default is \code{FALSE}
#' @param DEBUG \code{logical} value for usage of debugging, if \code{TRUE} the input parameters and the value of negative
# 'log-likelihood are printed on console for each iteration during optimization.
#' @param method \code{character} defining the method to be used in \code{optim}, preferences are: "Nelder-Mead", "BFGS", "L-BFGS-B"e
#' @param lower \code{vector} specifying the lower boundary of parameters for "L-BFGS-B" method
#' @param upper \code{vector} specifying the upper boundary of parameters for "L-BFGS-B" method
#' @param plot \code{logical} option of creating a plot of IDF curves with estimated parameters.
#' @param probs a vector of probabilities for which the IDF curves are calculated
#' @param cols a vector of colors for the seperate IDF curves, needs same length as \code{probs}
#' @param station.name \code{character} overall naming of the IDF plot, e.g. name of location or model name
#' @param data.name \code{character} naming the data points, e.g. obs or model name
#' @return $ints vector of sorted intensities for selected aggregation levels
#' @return $durs vector of sorted aggregation levels
#' @return $min minimum value of negative log-likelihood during optimization
#' @return $par vector of estimated IDF model parameters mu,sigma,xi,theta,eta at minimum value of negative log-likelihood.
#' @examples 
#' RR <- rgamma(10*30*24,shape=1)
#' year <- sort(rep(1:(10),30*24))
#' data <- data.frame(RR,year)
#' fit <- IDF.fit(data)
#' pars <- fit$par 
#' @author Christoph Ritschel \email{christoph.ritschel@@met.fu-berlin.de}

IDF.fit <- function(data,agg.lev=c(2,3,6,12,24,48,72,96),month="all",moving.sum="FALSE",theta.init=0,
                    use.log=FALSE,DEBUG=FALSE,method="Nelder-Mead",upper=Inf,lower=-Inf,plot=FALSE,
                    probs=c(0.5,0.9,0.99),cols=c(rgb(1,0,0,1),rgb(0,1,0,1),rgb(0,0,1,1)),station.name="Berlin",data.name="obs") {
  
  #########################################################################
  ### Calculate extreme values for each year and each aggregation level ###
  #########################################################################
  
  RR <- data$RR ## get precipitation time series from data.frame
  years <- unique(data$year) # get years from data.frame
  n.y <- length(years) # number of years
  n.a <- length(agg.lev) # number of aggregation times
  
  ## initilise arrays 
  agg.1 <- array(NA,dim=c(n.y)) 
  ints <- array(NA,dim=c(n.y*n.a))
  
  ###loop over years
  for(y in 1:n.y) {
    
    if(month[1]=="all") { 
      index <- which(data$year==years[y])
    }else if(is.integer(month) | is.numeric(month)) {
      index <- which(data$year==years[y] & data$mon >= min(month) & data$mon <= max(month))    
    }
    if(length(index)>0) {
    RR.year <- RR[index]
    agg.1[y] <- max(RR.year,na.rm=T) 
    
    ###loop over agg.lev
    for(a in 1:n.a) {
      
      ints[y+((a-1)*n.y)] <-  max(TS.acc(RR.year,agg.lev[a],moving.sum=moving.sum),na.rm=T)/agg.lev[a]
      
    } # end for all aggregation times
  } # end if lenght
  } # end for all years
  
  ## vector of all intensities
  int.vec <- c(agg.1,ints)
  
  ## vector of all durations (single)
  d.all <- c(1,agg.lev)
  ## long vector of all durations (repeated for each year to have same length as intensity vector)
  durs <- rep(d.all,each=n.y)
  
  ###############################################
  ### Estimate Parameters for single duration ###
  ###############################################

  ## Fit a generalized extreme value distribution to the maximum intensities of each year for a single 
  ## aggregation level and write the estimated parameters in an array for further analyisis. 

  ints.all <- matrix(int.vec,nrow=n.y) ## sort intensities in a matrix, rows are years, columns are aggregation levels
  pars <- array(NA,dim=c(3,length(d.all)))
  
  ## In case of NA values the optimization fails, therefore years with NA values need to be removed.
  ints.all <- matrix(ints.all[rowSums(!is.na(ints.all)) == length(d.all)],ncol=length(d.all))
  
  if(length(ints.all)==0) {
    cat("Warning: optimization did not converge and no parameters were estimated. Time Series contains too many NA values. \n")
  }else{
  
  ## loop over all aggregation levels
  for(d in 1:length(d.all)) {
    
    #fit <- fit.fun.emp(obs=ints.all[,d],mu=mu,sigma=sigma,xi=xi,use.log=use.log,
    #                   DEBUG=DEBUG,method=method,upper=upper,lower=lower)
    if(method=="L-BFGS-B"){
    fit <- gev.fit(xdat=ints.all[,d],method=method,upper=upper,lower=lower,show=F)
    }else{
      fit <- gev.fit(xdat=ints.all[,d],method=method,show=F)
    }
    
    pars[,d] <- fit$mle
    
  } ## end loop over aggregation levels
  
  #############################################################
  ### Derive starting parameters for duration-dependent GEV ###
  #############################################################
  
  ## Fit a linear model to the individual sigmas for individual aggregation times in a log-log environment
  ## The slope coefficient is an estimate for the slope in the duration-dependent GEV, namely parameter eta
  ## The intersection is an estimation of the starting parameter sigma
  ## Parameter mu is estimated as mean value of individual mus divided by indiviudal sigmas
  ## The initial value for xi will be the mean of all individual xi, since it is approximately independent of duration
  formel <- lm(log(pars[2,]) ~ log(d.all))
  sigma.est <- as.numeric(exp(formel$coefficients[1]))
  mu.est <- mean(pars[1,]/pars[2,])
  eta.est <- as.numeric(-formel$coefficients[2])
  
  xi.est <- max(0,mean(pars[3,],na.rm=T))

  ######################################################
  ### Estimate parameters for duration-dependent GEV ###
  ######################################################
  
  fit <- fit.fun(obs=int.vec,dur=durs,mu=mu.est,sigma=sigma.est,xi=xi.est,theta=theta.init,eta=eta.est,use.log=use.log,
                 DEBUG=DEBUG,method=method,upper=upper,lower=lower)
  
  if(plot&& !is.na(fit$min)) {
    ds <- sort(rep(d.all,length(int.vec)/length(d.all)))
    IDF.plot(pars=fit$par,probs,st.name=station.name,dt.name=data.name,ints=int.vec,ds=durs)
  }
  
  }
  
  if(!plot && is.na(fit$min)) {
    cat("Warning: optimization did not converge and no parameters were estimated. \n")
  }
  
  if(plot && is.na(fit$min)) {
    cat("Warning: optimization did not converge and no parameters were estimated. Plot not possible. \n")
  }
  
  return(list("ints"=int.vec,"durs"=durs,"min"=fit$min,"par"=fit$par))
  
} ## End of function IDF.fit
######################################################################################################################

########################################################################################################
#' @title Plotting IDF curves
#' @description The function \code{IDF.plot} plots a set of IDF curves with given IDF model parameters \code{pars} for
#' several probability levels \code{probs} at given durations \code{dur}. The colors of the curves can be defined with
#' parameter \code{cols} (need to have same length as \code{probs}). The \code{station.name} will be printed in the legend.
#' @param pars a vector of IDF model parameters mu,sigma,xi,eta,theta
#' @param probs a vector of probabilities for which the IDF curves are calculated
#' @param dur a vector of durations at which the IDF curves are calculated
#' @param cols a vector of colors for the seperate IDF curves, needs same length as \code{probs}
#' @param st.name \code{character} overall naming of the IDF plot, e.g. name of location or model name
#' @param dt.name \code{character} naming the data points, e.g. obs or model name
#' @param ints \code{vector} of observational intensities (surted by durations)
#' @param ds \code{vector} of durations (same length as intensities)
#' @examples 
#' RR <- rgamma(10*30*24,shape=1)
#' year <- sort(rep(1:(10),30*24))
#' data <- data.frame(RR,year)
#' fit <- IDF.fit(data)
#' param <- fit$par
#' IDF.plot(pars=param,st.name="example",dt.name="rgamma")
#' @author Christoph Ritschel \email{christoph.ritschel@@met.fu-berlin.de}

IDF.plot <- function(pars,probs=c(0.5,0.9,0.99),dur=c(0.5,1,2,3,6,12,24,48,72,96),cols=c(rgb(1,0,0,1),rgb(0,1,0,1),rgb(0,0,1,1)),st.name="Berlin-Dahlem",dt.name="obs",ints=NA,ds=NA) {
  
## initialize array for IDF values at different durations and for different probabilities
  idf.array <- array(NA,dim=c(length(dur),length(probs)))
  
  ## loop over probabilities
  for(i in 1:length(probs)) {
    
    ## calculate IDF values for given probability at all durations
    idf.array[,i] <- qgev.d(probs[i],mu=pars[1],sigma=pars[2],xi=pars[3],theta=pars[4],eta=pars[5],d=dur)
    
  } ## end of loop over probs
  
  ## initiialize plot window with limits of IDF values
  plot(NA,axes=F,xlim=c(min(dur,na.rm=T),max(dur,na.rm=T)),ylim=c(min(idf.array[,1],na.rm=T),max(idf.array[,3],na.rm=T)),xlab="duration [h]",ylab="intensity [mm/h]",log="xy")
  axis(1,at=dur,labels=dur)
  axis(2)  
  points(ds,ints,pch=16,col=rgb(0,0,0,0.5))
  
  ## loop over probabilities
  ## plot IDF curve
  for(i in 1:length(probs)) {
    points(dur,idf.array[,i],type="l",col=cols[i],lwd=1.5)
  } 
  
  legend.text.2 <- "quantile"
  
  ## plot legend
  legend(x="topright",legend=c(st.name,dt.name,paste(probs,legend.text.2,sep=" ")),
         col=c(1,rgb(0,0,0,0.5),cols),lty=c(NA,NA,rep(1,length(cols))),pch=c(NA,16,rep(NA,length(cols))))
  
} ## end of function IDF.plot
###################################################################################


