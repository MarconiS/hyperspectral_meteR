#' from coordinates, extract data for each of the data-products of interest
#'
#'
#' @return
#' @export
#' @examples
#' @import tidyverse
#' @importFrom magrittr "%>%"
download_point_daymet <- function(listSites = NULL,
                                  field_data = NULL,
                                  outpath = './tmp/climate/',
                                  startdate=1998,
                                  enddate= 2018){
  #field_data <- dplyr::select(field_data, individualID, latitude, longitude)
  if(is.null(field_data[["decimalLatitude"]]) || is.na(field_data[["decimalLongitude"]])){
    new_dat <- get_lat_long(field_data)
    daymet_coords <- cbind(as.character(field_data[["individualID"]]),
                           new_dat[["longitude"]],
                           new_dat[["latitude"]]) %>% unique
  }else{
    daymet_coords <- cbind(as.character(field_data[["individualID"]])
                           , field_data[["decimalLatitude"]]
                           , field_data[["decimalLongitude"]]) %>%
     unique
  }
  readr::write_csv(data.frame(daymet_coords), './tmp/daymet_metadata.csv')

  library(daymetr)
  download_daymet_batch(file_location = './tmp/daymet_metadata.csv',
                        start = startdate,
                        end = enddate,
                        internal = F,
                        path = outpath)
}
