#' Prediction function
#'
#' @title Prediction
#' @description Predict based on models and climate data.
#' @usage prediction(model, trainenv, vars, pre_infor, dir_out, Format=".tif$")
#' @param model ecological niche model
#' @param trainenv SpatRaster. The environmental data which is used in model.
#' @param vars Environment variables are identified in the environment data.
#' @param pre_infor data frame. the information of future environmental data.
#' @param dir_out Path of results.
#' @param Format Format of climate data is GeoTIFF(.tif).
#'
#' @return Return the SpatRaster of the corresponding time within the given
#' storage path.
#' @export
#'
#' @examples
#' \dontrun{
#' dir_out <- getwd()
#' model <- readRDS(system.file("extdata", "model.rda",
#'                                  package = "dispersalENM"))
#' pre_infor <- readRDS(system.file("extdata", "pre_infor.rda",
#'                                  package = "dispersalENM"))
#' trainenv <- terra::rast(list.files(
#' system.file("extdata", "env/trainenv", package = "dispersalENM"),
#' pattern = "tif", full.names = TRUE))
#' file.copy(system.file("extdata", "env/predictenv", package = "dispersalENM"), dir_out,
#'                      recursive = TRUE)
#' vars <- c("tasmin", "tasmax", "pre")
#' prediction(model, trainenv, vars, pre_infor, dir_out)
#' }




prediction <- function(model, trainenv, vars, pre_infor, dir_out, Format=".tif$") {
  dir.create(sprintf("%s/results", dir_out))
  dir.create(sprintf("%s/results/tif", dir_out))
  print("current")
  prediction <- terra::predict(trainenv, model, na.rm = TRUE)
  terra::writeRaster(prediction,
    filename = sprintf("%s/results/tif/current.tif", dir_out),
    overwrite = TRUE
  )

  df <- pre_infor[,-which(names(pre_infor) %in% c("EnvPath", "vars"))]
  for (i in c(1:nrow(pre_infor))) {
    time <- pre_infor$Time[i]
    print(time)
    dir_env_future <- pre_infor$EnvPath[i]
    all_files <- list.files(dir_env_future, pattern = Format)
    for(j in (1:ncol(df))) {
      pattern <- df[i,j]
      if(j == 1){
        tif_files <- all_files[grepl(pattern, all_files)]
      }else{
        tif_files <- tif_files[grepl(pattern, tif_files)]
      }
    }
    match_string  <-  paste0("(", paste(vars, collapse = "|"), ")")
    filtered_files <- tif_files[stringr::str_extract(tif_files, match_string) %in% vars]
    sorted_files <- filtered_files[order(match(
      stringr::str_extract(filtered_files, match_string), vars))]
    stacked_rasters <- rast(sprintf("%s/%s", dir_env_future, sorted_files))
    names(stacked_rasters) <- names(trainenv)
    pred_future <- terra::predict(stacked_rasters, model, na.rm = TRUE)
    output_file <- sprintf("%s/results/tif/%s_%s.tif", dir_out,
                           pre_infor$EnvScene[i], time)
    terra::writeRaster(pred_future, filename = output_file, overwrite = TRUE)
  }
}
