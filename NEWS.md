# IDF 2.1.1

## New features
The function gev.d.diag, which creates diagnostic plots, now provides the option to 
add (95% bootstrapped) confidence intervals to the plot.

We implemented the analytic likelihood gradient into the function gev.d.fit. When using the 
optimization method 'BFGS', it is now possible to obtain estimates, comparable to those of the 
method 'Nelder-Mead' (default), but with considerably reduced computation time. This is especially 
recommended, when using a model with covariates.

## Bug fixes
We fixed a bug, occurring when using different plotting symbols for different durations in the function gev.d.diag.

We changed the definition of the parameter eta2 in the documentation, so that it is in accordance 
with the corresponding publication.

# IDF 2.1.0

## New features:

The package now enables the usage of multiscaling and flattening in IDF curves. All functions are adapted to these new features. With default arguments, the new features are not used.

## Bug fixes

- implemented option to run IDF.agg() without parallelization
- corrected sign in gev.d.lik
   - before it was returning the negative (log) likelihood

# IDF 2.0.0

The package was extended to allow generalized linear modeling of the d-GEV parameters. This can be used to model for example spatial variations of the parameters. 
The package was extensively revised and restructured and some functions were removed as they were too specific. 


# IDF 1.0.0

R package for maximum likelihood fitting of duration-dependent generalized extreme value distribution (d-GEV). 
Additional functions for processing data (obtaining annual maxima) plotting IDF curves.