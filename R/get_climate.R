#' from coordinates, extract data for each of the data-products of interest
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"

get_climate <- function(field_data = NULL, tmppath = './tmp/climate/'
                        , startdate = 1998
                        , enddate = 2018
                        , provider = "daymet"
                        ){
  if(provider=="daymet"){
    download_point_daymet(field_data = field_data,
                          outpath = tmppath,
                          startdate=startdate,
                          enddate = enddate)
    melt_daymet_trends(path = tmppath)
    climate  = readr::read_csv("./indir/climate_features.csv")
    return(climate)
  }
}
