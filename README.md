
<!-- README.md is generated from README.Rmd. Please edit that file -->
DEPM
====

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/jonathansmart/DEPM.svg?branch=master)](https://travis-ci.org/jonathansmart/DEPM) <!-- badges: end -->

Installation
------------

You can install the released version of DEPM from [Github](https://github.com/jonathansmart/DEPM) with:

``` r
# install.packages("devtools")
devtools::install_github("jonathansmart/DEPM")
```

The DEPM package is a tool to estimate Spawning Stock Biomass for fish populations using Daily Egg Production Methods (DEPM). Two methods are available in this package:

1.  DEPM - the traditional DEPM which uses mean weight and mean fecundity to scale the number of females to total biomass
2.  DEPMWt - a new approach to DEPM presented in McGarvey et al. (in review), where female weight and fecundity are grouped in weight bins. This approach is useful for demersal finfish whose weight distributions are often not unimodal and therefore shouldn't be described using a normal distribution.

The package contains functions that estimate the individual DEPM parameters: daily egg production (*P*<sub>0</sub>), sex ratio (*R*), spawning fraction (*S*), fecundity (both mean and weight specific) and weight (both mean and in weight bins). The only parameter that is not estimated in this package is spawning area (*A*), which is typically estimated using GIS methods and must be provided to the analysis.

How-to guides
=============

Two vignettes are available that wil guide users through both approaches and can be accessed using:

``` r
browseVignettes("DEPM") # access all package vignettes

vignette("DEPMWt_Estimation") # How-to guide for the DEPMWt approach

# Vignette for the Traditional DEPM approach is forthcoming.
```
