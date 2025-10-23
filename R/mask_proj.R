#' mask_proj function
#'
#' @title base map
#' @description Create a base map of continuous numbers from 0 to n.
#' @usage mask_proj(env_1)
#' @param env_1 SpatRaster. One of the environmental data.
#'
#' @return SpatRaster
#' @export
#'

mask_proj <- function(env_1) {
  mask <- env_1
  terra::values(mask) <- seq_along(terra::values(mask))
  mask <- terra::mask(mask, env_1)
  return(mask)
}
