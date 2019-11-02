#' from coordinates, extract data for each of the data-products of interest
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"

get_dataset_from_point <- function(treeid = NULL,
                                   prdr = c("_CHM", "_DSM", "_DTM",
                                             "lope", "pect", "hsi"),
                                   inpath = "~/Documents/Data/Chapter3"){
  #get list of products that could be used to extract information from
  flrds = list.dirs(paste(inpath, "plots", sep="/"))
  flrds = grep(paste(prdr, collapse = "|"), flrds, value=T)
  for(ii in flrds){
    var_nm = substr(ii, nchar(ii)-2, nchar(ii))

    #get files path
    pt = list.files(ii, treeid[["plotID"]], full.names = T)
    if(length(pt)==0){
      if(var_nm %in% c("kld", "hsi")){
        warning("missing hyperspectral data for plot", treeid[["plotID"]])
        return(NULL)
      }else{
        pix_val = NA
      }
    }else{
      #read file in raster brick
      r = raster::brick(pt)

      #extract pixel of interest
      pix_val = raster::extract(r, treeid)
    }

    if(var_nm %in% c("kld", "hsi")){
      #add value to the sf
      var_nm = paste(var_nm, 1:length(pix_val), sep="_")
      for(bnd in 1:length(pix_val)){
        treeid[[var_nm[bnd]]] = as.numeric(pix_val)[bnd]
      }
    }else{
      #add value to the sf
      treeid[[var_nm]] = as.numeric(pix_val)
    }
  }
  return(treeid)
}
