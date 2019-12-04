#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @export
#' @examples
#' @importFrom magrittr "%>%"

#
get_soil_features <- function(plt){
  #get position of trees and the extent the extract soil data for
  centers <- sf::st_as_sf(plt, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

  bbclip <- sf::st_bbox(centers) %>% as.numeric
  names(bbclip) <- c("min_long", "min_lat", "max_long", "max_lat")
  bbclip = cbind.data.frame(bbclip[["min_long"]] - 0.05
                       , bbclip[["min_lat"]] - 0.05
                       , bbclip[["max_long"]] + 0.05
                       , bbclip[["max_lat"]] + 0.05)
  names(bbclip) <- c("min_long", "min_lat", "max_long", "max_lat")

  #bbclip <- bbclip[c(1,3,2,4)]
  #harmonize the SSURGO clip with the vegetation structure data
  soil_geometry_to_clip <- soilDB::mapunit_geom_by_ll_bbox(bbclip)
  soil_geometry_to_clip <- sf::st_as_sf(soil_geometry_to_clip)
  sf::st_crs(soil_geometry_to_clip)<- 4326
  #clip
  soil_features <- sf::st_join(centers, soil_geometry_to_clip, join = sf::st_intersects)
  #soil_features$plotID<- plt$plotID
  #colnames(soil_features)<- c("individualID", paste("soil", 1:(ncol(soil_features)-1), sep="."))
  return(soil_features)
}
