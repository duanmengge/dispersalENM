# occs: longitude, latitude, data.frame
# env: All environmental variable factors, SpatRaster
#' Title
#' @title Correlation analysis of environmental data
#' @description Perform correlation analysis on environmental data to avoid
#' multicollinearity.
#' @usage env_correlation(occs, env)
#' @param occs Data.frame. species occurrence data.
#' @param env SpatRaster. environmental data.
#'
#' @return vector
#' @export
#'
#' @examples
#' occourence <- readRDS(system.file("extdata", "occourence.rda",
#'                                  package = "dispersalENM"))
#' trainenv <-  terra::rast(list.files(
#'   system.file("extdata", "env/trainenv", package = "dispersalENM"),
#'   pattern = "tif", full.names = TRUE))
#' data <- env_correlation(occs = occourence, env = trainenv)
#'

env_correlation <- function(occs, env) {
  if (!all(c("x", "y") %in% names(occs))) {
    names(occs)[grep("^lon", names(occs))] <- "x"
    names(occs)[grep("^lat", names(occs))] <- "y"
  }
  predictors_p <- data.frame(terra::extract(env, occs[, c("x", "y")],
                                            ID = FALSE))
  per <- as.data.frame(stats::cor(predictors_p, method = "pearson"))
  return(per)
}
