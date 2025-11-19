<!-- README.md is generated from README.Rmd. Please edit that file -->

# dispersalENM

<!-- badges: start -->
<!-- badges: end -->

The goal of dispersalENM is to construct a species distribution model.

## Installation

You can install the development version of Dispersal like so:

``` r
devtools::install_github(duanmengge/dispersalENM)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
#Load dispersalENM package
library(dispersalENM)

#load data
occ_data <- readRDS(system.file("extdata", "occ_data.rda",
package = "dispersalENM"))
trainenv <-  terra::rast(list.files(
  system.file("extdata", "env/trainenv", package = "dispersalENM"),
  pattern = "tif", full.names = TRUE))

#Species occurrence data filtering
occourence <-occ_filter(occs = occ_data[,c("lon","lat")],
                        env_1 = trainenv[[1]])

#Correlation analysis of environmental data
data <- env_correlation(occs = occourence, env = trainenv)

#Select the optimal parameter among the alternative parameters of MaxEnt model
model_tuning_results <- model_tuning(occourence, trainenv, nbg=10000,
                                       f_c=c("L","Q"), r_m=c(1, 1, 3))

#MaxEnt model
model <- maxent_model(occourence, trainenv, fc=c("L"), nbg=10000,
                      rm=1, doclamp="false")
```

Predict based on models and climate data.
``` r
dir_out <- getwd()
model <- readRDS(system.file("extdata", "model.rda",
                                 package = "dispersalENM"))
pre_infor <- readRDS(system.file("extdata", "pre_infor.rda",
                                 package = "dispersalENM"))
trainenv <- terra::rast(list.files(
system.file("extdata", "env/trainenv", package = "dispersalENM"),
pattern = "tif", full.names = TRUE))
file.copy(system.file("extdata", "env/predictenv", package = "dispersalENM"), dir_out,
                     recursive = TRUE)
vars <- c("tasmin", "tasmax", "pre")
prediction(model, trainenv, vars, pre_infor, dir_out)
```
