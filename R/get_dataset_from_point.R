#' from coordinates, extract data for each of the data-products of interest
#'
#'
#' @return
#' @export
#' @examples
#' @importFrom magrittr "%>%"
#treeid = dat[9,]
get_dataset_from_point <- function(treeid = NULL,
                                   prdr = c("hsi"),
                                   itc = T,
                                   inpath = "//ufrc/ewhite/s.marconi/Chapter3/neonVegWrangleR/outdir/"){
  #get list of products that could be used to extract information from
  #treeid = sf::st_as_sf(treeid, coords=c("decimalLongitude", "decimalLatitude"), crs=4326)
  epsg= get_epsg_from_utm(treeid$utmZone)
  treeid = sf::st_as_sf(treeid, coords=c("easting", "northing"), crs=as.numeric(epsg))
  
  flrds = paste(paste(inpath, "plots", sep="/"), prdr, sep="/")
  #flrds = grep(paste(prdr, collapse = "|"), flrds, value=T)
  #flrds = "//ufrc/ewhite/s.marconi/Chapter3/neonVegWrangleR/outdir/plots/hsi"
  indvd_attr = list()
  for(ii in flrds){
    var_nm = substr(ii, nchar(ii)-2, nchar(ii))
    #get files path
    if(itc==T){
      pt = list.files(ii, paste(treeid[["individualID"]], ".tif", sep="")
                      , full.names = T)
    }else{
      pt = list.files(ii, treeid[["plotID"]], full.names = T)
    }
    if(length(pt)>1){
      pt = list.files(ii, treeid[["individualID"]], full.names = T)
    }
    #print(paste("path:", pt))
    
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
    }
    if(var_nm %in% c("kld", "hsi")){
      #add value to the sf
      # var_nm = paste(var_nm, 1:length(pix_val), sep="_")
      # for(bnd in 1:length(pix_val)){
      #   treeid[[var_nm[bnd]]] = as.numeric(pix_val)[bnd]
      # }
      
      #extract pixel of interest
      pix_val = raster::as.matrix(r)
      #dim(pix_val) = c(dim(pix_val)[1]*dim(pix_val)[2], dim(pix_val)[3])
      pix_val = data.frame(individualID = treeid[["individualID"]], pix_val)
      colnames(pix_val) = c("individualID", paste("band", 1:369, sep="_"))
      indvd_attr[[var_nm]] = pix_val
      
    }else{
      #add value to the sf
      pix_val = raster::as.matrix(r)
      pix_val = data.frame(individualID = treeid[["individualID"]], pix_val)
      indvd_attr[[var_nm]] = pix_val
    }
  }
  indvd_attr = do.call(cbind.data.frame, indvd_attr)
  todrop = indvd_attr[-1] %>% dplyr::select((contains("individualID"))) %>% colnames
  indvd_attr = indvd_attr  %>% dplyr::select(-one_of(todrop))
  #indvd_attr %>% reduce(inner_join, by = "individualID")
  return(indvd_attr)
}
