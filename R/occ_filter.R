# 1. 数据清理
# occs: longitude, latitude, data.frame
# env_1: a environmental variable factor, SpatRaster
#' occfilter function
#'
#' @title Species occurrence data filtering
#' @description To avoid bias, process species occurrence data. 1. Remove data
#' that is missing x or y; 2. Only one duplicate data is retained; 3. Only one
#' piece of data is retained within each grid.
#' @usage occ_filter(occs, env_1)
#' @param occs Data.frame. species occurrence data.
#' @param env_1 SpatRaster. One of the environmental data.
#'
#' @return data frame
#' @export
#'
#' @examples
#' occ_data <- readRDS(system.file("extdata", "occ_data.rda",
#' package = "dispersalENM"))
#' trainenv <-  terra::rast(list.files(
#'   system.file("extdata", "env/trainenv", package = "dispersalENM"),
#'   pattern = "tif", full.names = TRUE))
#' occourence <-occ_filter(occs = occ_data[,c("lon","lat")],
#'                         env_1 = trainenv[[1]])
#'

occ_filter <- function(occs, env_1) {
  if (all(c("x", "y") %in% names(occs))) {
    print(colnames(occs))
  } else {
    names(occs)[grep("^lon", names(occs))] <- "x"
    names(occs)[grep("^lat", names(occs))] <- "y"
    print(colnames(occs))
  }
  mask <- env_1
  terra::values(mask)[!is.na(terra::values(mask))] <-
    seq_along(terra::values(mask)[!is.na(terra::values(mask))])
  occs <- terra::extract(mask, occs[, c("x", "y")], xy = TRUE, ID = FALSE)
  colnames(occs) <- c("index", "x", "y")
  occs <- occs[!is.na(occs), ]
  occs_unique <- unique(occs$index)
  points <- terra::as.data.frame(mask, xy = TRUE)
  colnames(points)[3] <- "mask"
  selected_points <- points[points$mask %in% occs_unique, ]
  selected_points <- selected_points[, 1:2]
  return(selected_points)
}
