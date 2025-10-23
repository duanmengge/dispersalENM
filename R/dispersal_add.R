### Diffusion distance probability —— Exponential function

####distri####
#' distri function
#'
#' @title Diffusion distribution
#' @description Diffusion distribution
#' @param distribution Diffusion distribution
#' @return value
#' @export
#'
#'
distri <- function(distribution) {
  if (distribution == "exponential") {
    density <- stats::rexp(n = 10000, rate = 0.1)
    density <- density / max(density)
    density <- density[density <= 0.5]
    density <- density * 2
  }else if (distribution == "uniform") {
    density <- rep(1, 10000)
  }
  return(density)
}



####get_disp####
#' get_disp function
#'
#' @title Calculate the diffusion distance in each direction
#' @description Each cell has n diffusion directions, and the diffusion distance
#' of species in each direction is calculated based on the maximum diffusion
#' distance and exponential distribution.
#' @param n numeric. The number of directions a cell diffusion. the default is
#' 4.
#' @param max_disp numeric. The maximum diffusion distance of a species per unit
#' time.
#' @param pixel numeric. the resolution of environmental data.
#' @param R numeric. the resistance of environmental data.
#' @param density the species dispersal distance density distribution.
#'
#' @return vector
#' @export
#'
#'

get_disp <- function(n, max_disp, pixel, density, R) {
  dist_value <- density * max_disp *(1 - R) + pixel / 2 * sqrt(2)
  dist_value[sample(length(density), n, replace = TRUE)]
}

####afun####
#' afun function
#'
#' @title Estimate the coefficients of the driving factors
#' @description Calculate the proportional coefficients of the driving factors.
#' @param k_value Mean intensity of the driving factor effect, with the default
#' value set to max_disp.
#' @param kv The data of driving factors along the x-axis.
#' @param ku The data of driving factors along the y-axis.
#'
#' @return a
#' @export
#'
#'

afun <- function(k_value, kv, ku){
  t <- sqrt(ku^2 + kv^2)
  t_mean  <- terra::global(t, "mean", na.rm = TRUE)
  a <- k_value / t_mean$mean
  return(a)
}


####expand_buffer####
#' expand_buffer function
#'
#' @title Calculate the buffer of a cell after diffusion
#' @description After species diffuse at n angles on a certain cell, n points
#' are connected to form a buffer, which is the diffusion buffer of a certain
#' cell.
#' @usage expand_buffer(point, n = 8, distances = NA, max_disp, pixel, crs)
#' @param point A point containing x and y
#' @param n numeric. The number of directions a cell diffusion. the default is
#' 4.
#' @param distances numeric. The diffusion distance. if NA,  calculated based
#' on the exponential distribution and maximum diffusion distance.
#' @param max_disp numeric. the maximum diffusion distance.
#' @param pixel numeric. the resolution of environmental data.
#' @param crs the points’ projection.
#' @param density the species dispersal distance density distribution.
#' @param a the coefficients of the driving factors.
#'
#' @return polygon
#' @export
#'
#'
expand_buffer <-
  function(point, n = 8, distances = NA, max_disp, density, pixel, crs,
           a = NULL) {
    if (is.na(distances)) {
      distances <- get_disp(n, max_disp, pixel, density, R = point$resistance)
    }
    step <- 2 * pi / n
    new_points<-data.frame(X_center=point$X, Y_center=point$Y,
                           angle=seq(0, 2*pi-step, by=step),
                           kv=point$kv, ku=point$ku,
                           distance=distances,row.names = NULL)
    colnames(new_points)<-c("X_center","Y_center","angle","kv","ku","distance")
    new_points$X<-new_points$X_center +
      cos(new_points$angle) * new_points$distance + a * new_points$kv
    new_points$Y<-new_points$Y_center +
      sin(new_points$angle) * new_points$distance + a * new_points$ku



    xys <- sf::st_as_sf(new_points, coords = c("X", "Y"), crs = crs)
    polys <- sf::st_cast(sf::st_combine(xys), "POLYGON")
    polys_df <- sf::st_sf(ID = 1, geom = polys, crs = crs)
    return(polys_df)
  }
