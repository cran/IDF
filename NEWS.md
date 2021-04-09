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