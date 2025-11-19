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

Using the dispersal ability of species to constrain future potential distribution
``` r
#Change to actual data frame
ALL_infor <- readRDS(system.file("extdata", "ALL_infor.rda",
                                 package = "dispersalENM"))
occourence <- readRDS(system.file("extdata", "occourence.rda",
                                 package = "dispersalENM"))
#Result storage path
dirout <- getwd()

#Copy the predicted raster files from the system
file.copy(system.file("extdata", "tif", package = "dispersalENM"), dirout,
                     recursive = TRUE)
file.copy(system.file("extdata", "barrier", package = "dispersalENM"), dirout,
                      recursive = TRUE)
#Dispersal
dispersal(all_infor = ALL_infor, occ = occourence, out_dir=dirout, n = 8,
figure = TRUE)
```
Statistical analysis of the results
``` r
#Change to actual data frame
folder_path <- sprintf("%s/%s", getwd(), "dispersal_tif")
out_dir <- getwd()
statistics(folder_path = folder_path, out_dir = out_dir)
```
