# gamlssx 1.0.2

## New features

* The `gamlss` package is now imported, so that it is not necessary to load `gamlss` before calling `fitGEV()`. 

## Bug fixes and minor improvements

* Added a description for the `fitGEV()` function, which was missing.
* Fixed CRAN check NOTE by adding dependency on R >= 4.1.0 in DESCRIPTION because package code in `gamlssx-internal.R` uses the function shorthand `\(...)`.
