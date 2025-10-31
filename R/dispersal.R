#' @title Using the dispersal ability of species to constrain future potential
#' distribution
#' @description This function adds the diffusion ability of species to the
#' potential distribution of binarization. Species can diffuse at different
#' frequencies and speeds, and combine them with potential suitable areas in
#' the future.
#' @usage dispersal(all_infor, occ, n = 8, out_dir=NULL, figure = TRUE)
#' @param all_infor data frame. it has seven columns ( "time", "path", "thr",
#' "max_disp","frequency_n", "exposure_threshold", "barrier") at least.
#' "time": numeric, a time of the future.
#' "path": logical, the path saved for the predicted raster (0-1) file at a
#' certain time in the future.
#' "thr": numeric, threshold for Binaring raster files.
#' "max_disp": numeric, The maximum diffusion distance of a species per unit
#' time.
#' "frequency_n": numeric, the quantity per unit time during a time quantum.
#' "exposure_threshold": How long can species survive in unsuitable
#' environments.
#' "barrier": logical, The path of raster data for the formation of obstacles
#' to species diffusion. Obstacle terrain includes altitude, water system, etc.
#' @param occ data.frame. Species occurrence data.
#' @param n numeric. The number of directions a cell diffusion. the default is
#' 4.
#' @param out_dir The complete path to the directory where the raster is written
#' as a. tif file. The diffusion results and exposure threshold results at each
#' time will be saved. if NULL, The result will be output in the working
#' directory.
#' @param figure logical. If TRUE, the intermediate results will plot. the
#' default is FALSE.
#' @param native logical. If TRUE, the species is native species,  the
#' default is TRUE.
#' @param k_value logical. Mean intensity of the driving factor effect. If
#' k_value is NULL, the default value set to max_disp.
#' @param buffer_dist logical. When native is FALSE, the current range of
#' species can be provided. If buffer_dist is NULL, the default range is 10
#' times the maximum distance from the occurrence points.
#' @param distribution Dispersal capacity distribution.
#' @return Return the SpatRaster of the corresponding time within the given
#' storage path.
#' @export
#'
#' @examples
#' #Change to actual data frame
#' ALL_infor <- readRDS(system.file("extdata", "ALL_infor.rda",
#'                                  package = "dispersalENM"))
#' occourence <- readRDS(system.file("extdata", "occourence.rda",
#'                                  package = "dispersalENM"))
#' #Result storage path
#' dirout <- getwd()
#'
#'#Copy the predicted raster files from the system
#' file.copy(system.file("extdata", "tif", package = "dispersalENM"), dirout,
#'                      recursive = TRUE)
#' file.copy(system.file("extdata", "barrier", package = "dispersalENM"), dirout,
#'                       recursive = TRUE)
#' #Dispersal
#' dispersal(all_infor = ALL_infor, occ = occourence, out_dir=dirout, n = 8,
#' figure = TRUE)
#'
dispersal <- function(all_infor, occ, n = 8, out_dir = NULL, figure = TRUE,
                      native = TRUE, buffer_dist = NULL, k_value  = NULL ,
                      distribution = "exponential") {
  dis_crs <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

  print(all_infor[1, 1])
  print(colnames(all_infor))
  if (is.null(out_dir)) {
    out_dir <- getwd()
  }else {
    out_dir <- out_dir
  }
  folder_path <- sprintf("%s/dispersal_tif", out_dir)
  dir.create(folder_path)
  folder_path_2 <- sprintf("%s/expose_tif", out_dir)
  dir.create(folder_path_2)

  cur_raster_1 <- terra::rast(sprintf("%s", all_infor[1, 2]))
  old_crs <- terra::crs(cur_raster_1)
  cur_raster <- terra::project(cur_raster_1, dis_crs)
  new_r <- terra::rast(terra::ext(cur_raster),
                       res = terra::res(cur_raster),
                       crs = terra::crs(cur_raster))

  if (figure) {
    terra::plot(cur_raster)
  }
  cur_binary <- cur_raster
  terra::values(cur_binary)[!is.na(terra::values(cur_binary))] <-
    ifelse(terra::values(cur_binary)[!is.na(terra::values(cur_binary))] >=
             all_infor$thr[1], 1, 0)
  c_binary <- cur_binary
  if ("barrier" %in% colnames(all_infor)) {
    if (!is.na(all_infor$barrier[1])) {
      barrier_10km <- terra::rast(sprintf("%s", all_infor[1, "barrier"]))
      barrier_10km <- terra::crop(barrier_10km, cur_raster_1, mask=T)
      barrier_10km <- terra::project(barrier_10km, dis_crs)
      barrier <- terra::resample(barrier_10km, new_r, method = "bilinear")
      terra::values(barrier)[terra::values(barrier) >= 1] <- 1
      barrier <- terra::match(barrier, 1)
    }
    cur_binary <- terra::mask(cur_binary, barrier)
  }else{
    cur_binary <- cur_binary
  }
  #Calculate the current potential distribution
  if(native == FALSE){
    max_disp <- all_infor$max_disp[1]
    if(is.null(buffer_dist)){
      buffer_dist = 10*max_disp
    }
    sf_points <- sf::st_as_sf(occ, coords = c("x", "y"), crs = old_crs)
    sf_points <- sf::st_transform(sf_points, crs = dis_crs)
    sf_points$value <- extract(cur_binary, vect(sf_points))[,2]
    sf_points <- subset(sf_points, value == 1)
    sf_points <- subset(sf_points, select = -value)
    tmp_sf_buff<-sf::st_buffer(sf_points, buffer_dist)
    pre_patch <- terra::mask(cur_binary, tmp_sf_buff)
    pre_patch <- terra::ifel(pre_patch == 1, 1, NA)
  }else{
    cur_patch <- terra::patches(cur_binary, zeroAsNA = TRUE, directions = 8)
    names(cur_patch) <- "patches"
    if (figure) {
      terra::plot(cur_patch)
    }
    sf_points <- sf::st_as_sf(occ, coords = c("x", "y"), crs = old_crs)
    sf_projected <- sf::st_transform(sf_points, crs = dis_crs)
    occ_new <- data.frame(sf::st_coordinates(sf_projected))
    occ_extract <- terra::extract(cur_patch, occ_new, ID = FALSE)
    truepatches <- unique(occ_extract$patches)
    truepatches <- truepatches[which(!is.na(truepatches))]
    pre_patch <- terra::match(cur_patch, truepatches)
  }
  if (figure) {
    terra::plot(pre_patch)
  }
  pre_patch <- terra::ifel(is.na(pre_patch),NA, 1)
  cur_dis <- pre_patch
  # expose
  expose <- terra::ifel(!is.na(c_binary), 999, c_binary)
  E_points <- terra::as.data.frame(expose, xy = TRUE)
  colnames(E_points) <- c("x","y","expose")
  val <- terra::extract(pre_patch, E_points[, c("x", "y")], ID=FALSE)
  E_points$new <- val[,1]
  E_points$expose <- ifelse(is.na(E_points$new), E_points$expose+1, 0)
  E_points$expose <- ifelse(E_points$expose > 999, 999, E_points$expose)
  E_points <- E_points[,c("x","y","expose")]
  exposeR <- terra::rast(E_points, type="xyz", crs=dis_crs)
  exposeR_old <- terra::project(exposeR, old_crs, method="near")
  exposeR_old <- terra::resample(exposeR_old, cur_raster_1,
                                 method="near")
  terra::writeRaster(exposeR_old,
                     sprintf("%s/%d_expose.tif", folder_path_2,
                             all_infor[1, 1]),
                     overwrite = TRUE)
  # mask
  mask_raster <- terra::project(cur_raster, dis_crs)
  mask_F <- mask_proj(mask_raster)
  names(mask_F) <-  "mask"

  dispersal_log <- list()
  j <- 2
  for (j in c(2:nrow(all_infor))) {
    print(all_infor[j, 1])
    pre_dis <- pre_patch
    max_disp <- all_infor$max_disp[j]
    diff_frequency <- all_infor$frequency_n[j]
    exposure_threshold <- all_infor$exposure_threshold[j]

    if ("barrier" %in% colnames(all_infor)) {
      if (!is.na(all_infor$barrier[j])) {
        barrier_10km <- terra::rast(sprintf("%s", all_infor[j, "barrier"]))
        barrier_10km <- terra::crop(barrier_10km, cur_raster_1, mask=T)
        barrier_10km <- terra::project(barrier_10km, dis_crs)
        barrier <- terra::resample(barrier_10km, new_r, method = "bilinear")
        terra::values(barrier)[terra::values(barrier) >= 1] <- 1
        barrier <- terra::match(barrier, 1)
        mask <- terra::mask(mask_F, terra::crop(barrier, terra::ext(mask_F)))
        mask_points <- terra::as.data.frame(mask, xy = TRUE)
        mask_points <- sf::st_as_sf(mask_points,
          coords = c("x", "y"),
          crs = sf::st_crs(mask)
        )
      }
    }else{
      mask <- mask_F
      mask_points <- terra::as.data.frame(mask, xy = TRUE)
      mask_points <- sf::st_as_sf(mask_points,
                                  coords = c("x", "y"),
                                  crs = sf::st_crs(mask)
      )
    }

    if ("resistance" %in% colnames(all_infor)) {
      if (!is.na(all_infor$resistance[j])) {
        resistance <- terra::rast(sprintf("%s", all_infor[j, "resistance"]))
        resistance <- terra::crop(resistance, cur_raster_1)
        resistance <- terra::project(resistance, dis_crs)
        resistance <- terra::resample(resistance, new_r, method = "bilinear")
      }
    }

    if (all(c("data_v", "data_u") %in% colnames(all_infor))) {
      if (!is.na(all_infor$data_v[j])) {
        data_v <- terra::rast(sprintf("%s", all_infor[j, "data_v"]))
        data_v <- terra::crop(data_v, cur_raster_1)
        data_v <- terra::project(data_v, dis_crs)
        data_v <- terra::resample(data_v, new_r, method = "bilinear")


        data_u <- terra::rast(sprintf("%s", all_infor[j, "data_u"]))
        data_u <- terra::crop(data_u, cur_raster_1)
        data_u <- terra::project(data_u, dis_crs)
        data_u <- terra::resample(data_u, new_r, method = "bilinear")
        if(is.null(k_value)){
          k_value <- max_disp
        }
        a_value <- afun(k_value, data_u, data_v)
      }
    }

    future_raster <- terra::rast(sprintf("%s", all_infor[j, "path"]))
    future_raster <- terra::project(future_raster, dis_crs)
    if (figure) {
      terra::plot(future_raster)
    }
    fu_binary <- future_raster
    terra::values(fu_binary)[!is.na(terra::values(fu_binary))] <-
      ifelse(terra::values(fu_binary)[!is.na(terra::values(fu_binary))]
             >= all_infor$thr[j], 1, 0)
    if ("barrier" %in% colnames(all_infor)){
      future_binary <- terra::mask(fu_binary,
                                   terra::crop(barrier, terra::ext(fu_binary)))
    }else{
      future_binary <- fu_binary
    }
    if (figure) {
      terra::plot(future_binary)
    }

    z <- 1
    pixel <- terra::res(future_binary)[1]
    for (z in c(1:diff_frequency)) {
      pre_lines <- sf::st_as_sf(terra::as.lines(terra::as.polygons(pre_patch)))
      pre_buffer <- sf::st_union(
        sf::st_buffer(pre_lines, dist = min(terra::res(pre_patch)) / 2))
      edge_buffer <- terra::mask(pre_patch, terra::vect(pre_buffer),
                                 touches = TRUE)
      shape <- sf::st_as_sf(terra::as.points(edge_buffer),
                            coords = c("x", "y"),
                            remove = FALSE,
                            crs = terra::crs(pre_patch))
      feather <- sf::st_sf(ID = seq_len(nrow(shape)), geom = shape$geometry)
      points <- sf::st_cast(feather$geom, "POINT")

      edge_polys <- NULL

      for (i in seq_along(points)) {
        point <- data.frame(sf::st_coordinates(points[i]))

        if(exists("resistance")){
          point$resistance <- terra::extract(resistance, point[,c("X","Y")],
                                             ID=FALSE)[, 1]
          point$resistance <- ifelse(is.na(point$resistance), 0, point$resistance)
        }else{
          point$resistance <- 0
        }

        if(exists("data_u") && exists("data_v")){
          point$ku <- terra::extract(data_u, point[,c("X","Y")], ID=FALSE)[,1]
          point$kv <- terra::extract(data_v, point[,c("X","Y")], ID=FALSE)[,1]
          a <- a_value
        }else{
          point$ku <- 0
          point$kv <- 0
          a <- 0
        }
        if(is.na(point$ku) | is.na(point$kv)){
          point$ku <- 0
          point$kv <- 0
          a <- 0
        }


        if (i == 1) {
          density <- distri(distribution)
        }
        buffer <- expand_buffer(point, n = n, max_disp = max_disp,
                                pixel = pixel, density = density,
                                crs = sf::st_crs(points), a = a)

        edge_polys <- if (is.null(edge_polys)) {
          buffer
        } else {
          rbind(edge_polys, buffer)
        }
      }
      edge_polys <- rbind(edge_polys, feather)
      if (figure) {
        terra::plot(edge_polys)
      }
      #### Internal diffusion
      pre_points <- terra::as.data.frame(pre_patch, xy = TRUE)[, c("x", "y")]
      #### edge diffusion
      edge_points <- terra::as.data.frame(edge_buffer, xy = TRUE)[, c("x", "y")]
      moveable_points <- dplyr::anti_join(pre_points, edge_points,
                                          by = c("x", "y"))
      moveable <- sf::st_as_sf(moveable_points, coords = c("x", "y"),
                               remove = FALSE, crs = terra::crs(pre_patch))
      moveable <- sf::st_sf(ID = seq_len(nrow(moveable)),
                            geom = moveable$geometry)
      move_buffer <- sf::st_buffer(moveable, max_disp)
      if (figure) {
        terra::plot(move_buffer)
      }

      all_polys <- rbind(edge_polys, move_buffer)
      geom <- sf::st_union(all_polys)
      if(any(st_geometry_type(geom) == "GEOMETRYCOLLECTION")){
        geom <- sf::st_collection_extract(geom, "POLYGON")
      }

      pts_buf_union <- sf::st_cast(geom, "POLYGON")
      if (figure) {
        terra::plot(pts_buf_union)
      }

      if (z != diff_frequency) {
        diff_union <- terra::vect(pts_buf_union)
        diff_raster <- terra::crop(future_binary, diff_union, mask=T)
        pre_patch <- terra::match(diff_raster, 1)
        if (figure) {
          terra::plot(pre_patch)
        }
      }
    }


    if (nrow(edge_points) != 0) {
      diff_union <- terra::vect(pts_buf_union)
      diff_raster <- terra::crop(future_binary, diff_union, mask=T)
      new_raster <- terra::ifel(diff_raster == 1, 2, diff_raster)
      new_patch <- terra::merge(new_raster, future_binary)

      terra::crs(new_patch) <- dis_crs
      new_patch_old <- terra::project(new_patch, old_crs, method="near")
      new_patch_old <- terra::resample(new_patch_old, cur_raster_1,
                                           method="near")
      terra::writeRaster(new_patch_old,
                         sprintf("%s/%d_dispersal.tif", folder_path,
                                 all_infor[j, 1]),
                         overwrite = TRUE)
    }
    pre2_patch <- terra::ifel(new_patch == 2, 1, NA)
    dispersal_log[[j]] <- new_patch
    #expose
    val <- terra::extract(pre2_patch, E_points[, c("x", "y")], ID=FALSE)
    E_points$new <- val[,1]
    E_points$expose <- ifelse(!is.na(E_points$new), 0, E_points$expose+1)
    E_points$expose <- ifelse(E_points$expose > exposure_threshold, 999,
                              E_points$expose)
    E_points <- E_points[,c("x","y","expose")]
    exposeR <- terra::rast(E_points, type="xyz", crs=dis_crs)
    exposeR_old <- terra::project(exposeR, old_crs, method="near")
    exposeR_old <- terra::resample(exposeR_old, cur_raster_1,
                                     method="near")
    terra::writeRaster(exposeR_old,
                       sprintf("%s/%d_expose.tif", folder_path_2,
                               all_infor[j, 1]),
                       overwrite = TRUE)
    pre_patch <- terra::ifel(exposeR <= exposure_threshold, 1, NA)
  }
  cur_dis <- terra::ifel(cur_dis  == 1, 2, NA)
  current_patch <- terra::merge(cur_dis, cur_binary, first = TRUE)
  current_patch_old <- terra::project(current_patch, old_crs, method="near")
  current_patch_old <- terra::resample(current_patch_old, cur_raster_1,
                                       method="near")
  terra::writeRaster(current_patch_old,
                     sprintf("%s/%d_dispersal.tif", folder_path,
                             all_infor[1, 1]), overwrite = TRUE)

  saveRDS(dispersal_log, sprintf("%s/dispersal.rda", out_dir))
}
