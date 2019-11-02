#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#'
get_epsg_from_utm <- function(utm){
  utm <- as.character(utm)
  utm <-  substr(utm,1,nchar(utm)-1)
  utm[as.numeric(utm)<10] <- paste('0', utm[as.numeric(utm)<10], sep="")
  epsg <- paste("326", utm, sep="")
  return(epsg)
}
