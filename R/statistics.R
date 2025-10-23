#' @title Statistical analysis of the results
#' @description Calculate the ratio of potential distribution area (compared to
#' the study area) within the time interval, gain, refuge, and loss of
#' potential distribution area.
#' @usage statistics(folder_path, out_dir)
#' @param folder_path species dispersal file("tif" format) path.
#' @param out_dir The complete path to the directory where the raster is written
#' as a. tif file. The diffusion results and exposure threshold results at each
#' time will be saved. if NULL, The result will be output in the working
#' directory.
#' @return Return the SpatRaster of the corresponding time within the given
#' storage path.
#' @export
#'
#' @examples
#' #Change to actual data frame
#' folder_path <- sprintf("%s/%s", getwd(), "dispersal_tif")
#' out_dir <- getwd()
#' statistics(folder_path = folder_path, out_dir = out_dir)
#'
statistics <- function(folder_path, out_dir) {
  dis_crs <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

  if (is.null(out_dir)) {
    out_dir <- getwd()
  }else {
    out_dir <- out_dir
  }
  out_path <- sprintf("%s/statistics_tif", out_dir)
  dir.create(out_path)

  filenames <- list.files(folder_path, pattern = "\\.tif$")
  numbers <- gsub("\\D", "", filenames)
  stat_Data <- data.frame()
  for(n in sort(numbers)){
    print(n)
    if(n == min(numbers)){
      name <- grep(n, filenames, value = TRUE)
      data <- terra::rast(sprintf("%s/%s", folder_path, name))
      KK <- terra::ifel(is.na(data), NA,
                        terra::ifel(data == 2, 4, 0))
      points <- terra::as.data.frame(KK, xy=T)
      colnames(points) <- c("x","y","value")
      All <- nrow(points)
      Keep <- nrow(points[which(points$value == 4),])
      stat_data <- data.frame(Year = n,
                              PDA = scales::percent((0 + Keep)/All,
                                                    accuracy=0.01),
                              Gain = scales::percent(0, accuracy=0.01),
                              Refuge = scales::percent(Keep/All, accuracy=0.01),
                              Loss = scales::percent(0, accuracy=0.01))
      stat_Data <- rbind(stat_Data, stat_data)
      terra::writeRaster(KK, sprintf("%s/%s_statistics.tif", out_path, n),
                         overwrite = TRUE)
      n_pre <- n
      next
    }
    name1 <- grep(n_pre, filenames, value = TRUE)
    name2 <- grep(n, filenames, value = TRUE)
    data1 <- terra::rast(sprintf("%s/%s", folder_path, name1))
    data2 <- terra::rast(sprintf("%s/%s", folder_path, name2))
  #Decre 3, Keep 4，Incre 5
    KK <- data1
    KK <- terra::app(c(data1, data2), fun = function(x) {
      a <- x[1]
      b <- x[2]
      if (is.na(a) || is.na(b)) return(NA)
      if (a < 2  && b < 2) return(0)
      if (a == 2 && b < 2) return(3)
      if (a == 2 && b == 2) return(4)
      if (a < 2 && b == 2) return(5)
    })
    KK_crs <-  terra::project(KK, dis_crs, method="near")
    points <- terra::as.data.frame(KK_crs, xy=T)
    colnames(points) <- c("x","y","value")
    All <- nrow(points)

    Decre <- nrow(points[which(points$value==3),])
    Keep <- nrow(points[which(points$value==4),])
    Incre <- nrow(points[which(points$value==5),])
    stat_data <- data.frame(Year = n,
                            PDA = scales::percent((Incre + Keep)/All,
                                                  accuracy=0.01),
                            Gain = scales::percent(Incre/All, accuracy=0.01),
                            Refuge = scales::percent(Keep/All, accuracy=0.01),
                            Loss = scales::percent(Decre/All, accuracy=0.01))
    stat_Data <- rbind(stat_Data, stat_data)
    terra::writeRaster(KK, sprintf("%s/%s_statistics.tif", out_path, n),
                overwrite = TRUE)
    n_pre <- n
  }
  utils::write.csv(stat_Data, sprintf("%s/ALL_statistics.csv", out_dir))
}
