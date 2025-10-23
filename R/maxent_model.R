#' MaxEntmodel function
#'
#' @title MaxEnt model
#' @description Building the MaxEnt model.
#' @usage maxent_model(occs, env, fc = c("L", "Q", "P"),
#'                    nbg = 10000, rm = 1, doclamp = "false")
#' @param occs Data.frame. species occurrence data.
#' @param env SpatRaster. environmental data.
#' @param nbg Numeric. The quantity of background data. Default is 10000.
#' @param fc Feature combination. the option including "L"(linear),
#' "Q"(quadratic), "P"(product), "T"(threshold), "H"(hinge). the default is
#' c("L","Q","P").
#' @param rm Regularization multiplier. The default is 1.
#' @param doclamp Logical. If "true", using Jackknife method. The default is
#' "false".
#' @return model
#' @export
#'
#' @importFrom rJava .jinit .jnew .jcall
#'
#' @examples
#' occourence <- readRDS(system.file("extdata", "occourence.rda",
#'                                  package = "dispersalENM"))
#' trainenv <-  terra::rast(list.files(
#' system.file("extdata", "env/trainenv", package = "dispersalENM"),
#' pattern = "tif", full.names = TRUE))
#' model <- maxent_model(occourence, trainenv, fc=c("L"), nbg=10000,
#'                       rm=1, doclamp="false")
#'

maxent_model <- function(occs, env, fc = c("L", "Q", "P"),
                         nbg = 10000, rm = 1, doclamp = "false") {
  featuretypes <- c(
    "linear", "quadratic", "product", "threshold",
    "hinge"
  )
  featuresimple <- c("L", "Q", "P", "T", "H")
  for (i in seq_along(featuretypes)) {
    if (length(grep(featuresimple[i], fc)) > 0) {
      assign(featuretypes[i], "true")
    } else {
      assign(featuretypes[i], "false")
    }
  }
  occs_z <- terra::extract(env, occs, ID = FALSE, xy = TRUE)
  env_data <- terra::as.data.frame(env, xy = TRUE, na.rm = TRUE)
  bg_z <- env_data[sample(nrow(env_data), nbg), ]
  occs_z$value <- 1
  bg_z$value <- 0
  raster_stack <- rbind(occs_z, bg_z)
  p <- raster_stack[, c("value")]
  x <- raster_stack[, -which(names(raster_stack) %in% c("x", "y", "value"))]
  model <- dismo::maxent(x, p,
    args = c(
      sprintf("betamultiplier=%.2f", rm),
      "autofeature=false",
      sprintf("linear=%s", linear),
      sprintf("quadratic=%s", quadratic),
      sprintf("product=%s", product),
      sprintf("threshold=%s", threshold),
      sprintf("hinge=%s", hinge),
      sprintf("doclamp=%s", doclamp)
    )
  )
  return(model)
}
