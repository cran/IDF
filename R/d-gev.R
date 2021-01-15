# This file contains the functions  dgev.d, pgev.d, qgev.d, rgev.d for the duration-dependent-gev.

#### dgev.d() ####

#' d-GEV probability density function
#'
#' @description Probability density function of duration-dependent GEV distribution
#' @param q vector of quantiles
#' @param mut,sigma0,xi numeric value, giving modified location \eqn{\tilde{\mu}}, scale offset \eqn{\tilde{\sigma_0}} and 
#' shape parameter \eqn{\xi}.
#' @param theta numeric value, giving duration offset \eqn{\theta} (defining curvature of the IDF curve)
#' @param eta numeric value, giving duration exponent \eqn{\eta} (defining slope of the IDF curve)
#' @param d positive numeric value, giving duration
#' @param ... additional parameters passed to \code{\link[evd]{dgev}}
#' 
#' @details For details on the d-GEV and the parameter definitions, see \link{IDF-package}.
#' 
#' @return list containing vectors of density values for given quantiles.
#' The first element of the list are the density values for the first given duration etc.
#' 
#' @seealso \code{\link{pgev.d}}, \code{\link{qgev.d}}, \code{\link{rgev.d}}
#' 
#' @export
#' @importFrom evd dgev 
#'
#' @examples
#' x <- seq(4,20,0.1)
#' # calculate probability density for one duration
#' dgev.d(q=x,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.1,d=1)
#' 
#' # calculate probability density for different durations
#' ds <- 1:4
#' dens <- lapply(ds,dgev.d,q=x,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.1)
#' 
#' plot(x,dens[[1]],type='l',ylim = c(0,0.21),ylab = 'Probability Density')
#' for(i in 2:4){
#'   lines(x,dens[[i]],lty=i)
#' }
#' legend('topright',title = 'Duration',legend = 1:4,lty=1:4)
dgev.d <- function(q,mut,sigma0,xi,theta,eta,d,...) {
  if(any(c(length(mut),length(sigma0),length(xi),length(theta),length(eta))>1)){
    message('One of the parameters mut, sigma0, xi, theta, eta is a vector. ',
            'This is not intended and might cause an error.')}
  if (d<=0) {stop('The duration d has to be positive.')}
  if(any(d+theta<=0)){
    warning('Some shape parameters are negative,resulting from a negativ theta '
            ,theta, ' this will prododuce NAs.')}
  # if denominator is negative NAs will be returned
  if(d+theta<=0){return(rep(NA,length(q)))}else{
    sigma.d <-sigma0/(d+theta)^eta
    return(dgev(q,loc=mut*sigma.d,scale=sigma.d,shape=xi,...))}
}


#### pgev.d() ####

#' d-GEV cumulative distribution function
#'
#' @description Cumulative probability distribution function of duration-dependent GEV distribution
#' @param q vector of quantiles
#' @param mut,sigma0,xi numeric value, giving modified location, modified scale and shape parameter
#' @param theta numeric value, giving duration offset (defining curvature of the IDF curve)
#' @param eta numeric value, giving duration exponent (defining slope of the IDF curve)
#' @param d positive numeric value, giving duration
#' @param ... additional parameters passed to \code{\link[evd]{pgev}}
#' 
#' @details The duration dependent GEV distribution is defined after 
#' [Koutsoyiannis et al., 1998]:
#' \deqn{G(x)= \exp[-\left( 1+\xi(x/\sigma(d)-\mu_t) \right)^{-1/\xi}] } 
#' with the duration dependent scale \eqn{\sigma(d)=\sigma_0/(d+\theta)^\eta} and 
#' modified location parameter \eqn{\mu_t=\mu/\sigma(d)}.
#' 
#' @details For details on the d-GEV and the parameter definitions, see \link{IDF-package}.
#' 
#' @return list containing vectors of probability values for given quantiles. 
#' The first element of the list are the probability values for the first given duration etc.
#' 
#' @seealso \code{\link{dgev.d}}, \code{\link{qgev.d}}, \code{\link{rgev.d}}
#' 
#' @export
#' @importFrom evd pgev 
#'
#' @examples
#' x <- seq(4,20,0.1)
#' prob <- pgev.d(q=x,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.1,d=1)
pgev.d <- function(q,mut,sigma0,xi,theta,eta,d,...) {
  if(any(c(length(mut),length(sigma0),length(xi),length(theta),length(eta))>1)){
    message('One of the parameters mut, sigma0, xi, theta, eta is a vector. ',
            'This is not intended and might cause an error.')}
  if (d<=0) {stop('The duration d has to be positive.')}
  if(any(d+theta<=0)){
    warning('Some shape parameters are negative,resulting from a negativ theta '
            ,theta, ' this will prododuce NAs.')}
  # if denominator is negative NAs will be returned
  if(d+theta<=0){return(rep(NA,length(q)))}else{
    sigma.d <-sigma0/(d+theta)^eta
    return(pgev(q,loc=mut*sigma.d,scale=sigma.d,shape=xi,...))}
}


#### qgev.d() ####

#' d-GEV quantile function 
#'
#' @description Quantile function of duration-dependent GEV distribution (inverse of the cumulative probability distribution function)
#' @param p vector of probabilities
#' @param mut,sigma0,xi numeric value, giving modified location, modified scale and shape parameter
#' @param theta numeric value, giving duration offset (defining curvature of the IDF curve)
#' @param eta numeric value, giving duration exponent (defining slope of the IDF curve)
#' @param d positive numeric value, giving duration
#' @param ... additional parameters passed to \code{\link[evd]{qgev}}
#' 
#' @details The duration dependent GEV distribution is defined after 
#' [Koutsoyiannis et al., 1998]:
#' \deqn{ G(x)= \exp[-\left( 1+\xi(x/\sigma(d)-\mu_t) \right)^{-1/\xi}] } 
#' with the duration dependent scale \eqn{\sigma(d)=\sigma_0/(d+\theta)^\eta} and 
#' modified location parameter \eqn{\mu_t=\mu/\sigma(d)}.
#' 
#' @return list containing vectors of quantile values for given probabilities. 
#' The first element of the list are the q. values for the first given duration etc.
#' 
#' @details For details on the d-GEV and the parameter definitions, see \link{IDF-package}.
#' 
#' @seealso \code{\link{pgev.d}}, \code{\link{dgev.d}}, \code{\link{rgev.d}}
#' 
#' @export
#' @importFrom evd qgev 
#'
#' @examples
#' p <- c(0.5,0.9,0.99)
#' # calulate quantiles for one duration
#' qgev.d(p=p,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.3,d=1)
#' 
#' # calculate quantiles for sequence of durations
#' ds <- 2^seq(0,4,0.1)
#' qs <- lapply(ds,qgev.d,p=p,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.3)
#' qs <- simplify2array(qs)
#' 
#' plot(ds,qs[1,],ylim=c(3,20),type='l',log = 'xy',ylab='Intensity',xlab = 'Duration')
#' for(i in 2:3){
#'   lines(ds,qs[i,],lty=i)
#' }
#' legend('topright',title = 'p-quantile',
#'        legend = p,lty=1:3,bty = 'n')
qgev.d <- function(p,mut,sigma0,xi,theta,eta,d,...) {
  if(any(c(length(mut),length(sigma0),length(xi),length(theta),length(eta))>1)){
    message('One of the parameters mut, sigma0, xi, theta, eta is a vector. ',
            'This is not intended and might cause an error.')}
  if (d<=0) {stop('The duration d has to be positive.')}
  if(any(d+theta<=0)){
    warning('Some shape parameters are negative, resulting from a negativ theta '
            ,theta, ' this will prododuce NAs.')}
  # if denominator is negative NAs will be returned
  if(d+theta<=0){return(rep(NA,length(p)))}else{
    sigma.d <-sigma0/(d+theta)^eta
    return(qgev(p,loc=as.numeric(mut*sigma.d)
                ,scale=as.numeric(sigma.d),shape=as.numeric(xi),...))}
}


#### rgev.d() ####

#' Generation of random variables from d-GEV
#'
#' @description Generation of random variables following duration-dependent GEV.
#' @param n number of random variables per duration
#' @param mut,sigma0,xi numeric value, giving modified location, modified scale and shape parameter
#' @param theta numeric value, giving duration offset (defining curvature of the IDF curve)
#' @param eta numeric value, giving duration exponent (defining slope of the IDF curve)
#' @param d positive numeric value, giving duration
#' 
#' @details For details on the d-GEV and the parameter definitions, see \link{IDF-package}
#' 
#' @return list containing vectors of random variables.  
#' The first element of the list are the random values for the first given duration etc.
#' Note that the random variables for different durations are nor ordered (contrary to precipitation maxima of different durations).
#' 
#' @seealso \code{\link{pgev.d}}, \code{\link{qgev.d}}, \code{\link{dgev.d}}
#' 
#' @export
#' @importFrom evd rgev 
#'
#' @examples
#' # random sample for one duration
#' rgev.d(n=100,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.3,d=1)
#' 
#' # compare randomn samples for different durations
#' ds <- c(1,4)
#' samp <- lapply(ds,rgev.d,n=100,mut=4,sigma0=2,xi=0,theta=0.1,eta=0.3)
#' 
#' hist(samp[[1]],breaks = 10,col=rgb(1,0,0,0.5),freq = FALSE
#'      ,ylim=c(0,0.3),xlim=c(3,20),xlab='x',main = 'Random d-GEV samples')
#' hist(samp[[2]],breaks = 10,add=TRUE,col=rgb(0,0,1,0.5),freq = FALSE)
#' legend('topright',fill = c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),
#' legend = paste('d=',1:2,'h'),title = 'Duration')
rgev.d <- function(n,mut,sigma0,xi,theta,eta,d) {
  if(any(c(length(mut),length(sigma0),length(xi),length(theta),length(eta))>1)){
    message('One of the parameters mut, sigma0, xi, theta, eta is a vector. ',
            'This is not intended and might cause an error.')}
  if (d<=0) {stop('The duration d has to be positive.')}
  if(any(d+theta<=0)){
    warning('Some shape parameters are negative, resulting from a negative theta '
            ,theta, ' this will prododuce NAs.')}
  # if denominator is negative NAs will be returned
  if(d+theta<=0){return(rep(NA,n))}else{
    sigma.d <-sigma0/(d+theta)^eta
    return(rgev(n,loc=mut*sigma.d,scale=sigma.d,shape=xi))}
}


