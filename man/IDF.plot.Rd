% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IDF.R
\name{IDF.plot}
\alias{IDF.plot}
\title{Plotting of IDF curves at a chosen station}
\usage{
IDF.plot(
  durations,
  fitparams,
  probs = c(0.5, 0.9, 0.99),
  cols = 4:2,
  add = FALSE,
  legend = TRUE,
  ...
)
}
\arguments{
\item{durations}{vector of durations for which to calculate the quantiles.}

\item{fitparams}{vector containing parameters mut, sigma0, xi, theta, eta
(modified location, scale offset, shape, duration offset, duration exponent) for chosen station
as obtained from \code{\link{gev.d.fit}}
(or \code{\link{gev.d.params}} for model with covariates).}

\item{probs}{vector of non-exceedance probabilities for which to plot IDF curves (p = 1-1/(Return Period))}

\item{cols}{vector of colors for IDF curves. Should have same length as \code{probs}}

\item{add}{logical indicating if plot should be added to existing plot, default is FALSE}

\item{legend}{logical indicating if legend should be plotted (TRUE, the default)}

\item{...}{additional parameters passed on to the \code{plot} function}
}
\description{
Plotting of IDF curves at a chosen station
}
\examples{
data('example',package = 'IDF')
# fit d-gev
fit <- gev.d.fit(example$dat,example$d,ydat = as.matrix(example[,c("cov1","cov2")])
                 ,mutl = c(1,2),sigma0l = 1)
# get parameters for cov1 = 1, cov2 = 1
par <- gev.d.params(fit = fit, ydat = matrix(1,1,2))
# plot quantiles
IDF.plot(durations = seq(0.5,35,0.2),fitparams = par)
# add data points
points(example[example$cov1==1,]$d,example[example$cov1==1,]$dat)
}
