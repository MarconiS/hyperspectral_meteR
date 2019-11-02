#' from coordinates, extract data for each of the data-products of interest
#'
#'
#' @return
#' @export
#' @import daymetr lubridate
#' @examples
#' @importFrom magrittr "%>%"

melt_daymet_trends <- function(path = "./tmp/climate/"
                               , outfile = "./indir/climate_features.csv"
                               ){

  dataset = data.frame(matrix(NA, ncol = 9, nrow = 0))
  colnames(dataset) <- c("month","ts_daylength", "ts_prec","ts_rad",
                         "ts_melt","ts_tmax","ts_tmin","ts_vp", "individualID")
  ls_dat = list.files(path, pattern = "csv")
  for(ii in ls_dat){
    id_clim <- read.csv(paste(path,ii, sep="/"), skip=7)
    colnames(id_clim) <- c("year", "month", "daylength", "prec", "srad", "snow_melt",
                           "tmax", "tmin", "vp")
    id_clim$month <- as.Date(id_clim$month-1, origin = "1995-01-01") %>%
      month
    point_features <- id_clim %>%
      dplyr::group_by(month, year) %>%
      dplyr::summarise(daylength = mean(daylength), prec = sum(prec),
                       srad = mean(srad), snow_melt = sum(snow_melt),
                       tmax = max(tmax), tmin = min(tmin), vp = mean(vp))

    point_features <- point_features %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(daylength = mean(daylength), prec = mean(prec),
                       rad = mean(srad), snow_melt = mean(snow_melt),
                       tmax = mean(tmax), tmin = mean(tmin), vp = mean(vp))
    point_features[["individualID"]] <- gsub('.{14}$', '', ii)
    dataset = rbind(dataset, point_features)
  }
  #write_csv(dataset, "./DMT_retriever/climate_features.csv")
  library(data.table) ## v >= 1.9.6
  clim_dat <- dcast(melt(dataset, id.vars=c("individualID", "month")),
                    individualID~variable+month)
  readr::write_csv(clim_dat, outfile)
  return(clim_dat)
}

