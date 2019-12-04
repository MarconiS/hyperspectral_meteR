# utility functions
#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @import rotl ape
#' @examples
#' @importFrom magrittr "%>%"
na.rm <- function(dat){
  na_rm <- rowSums(is.na(dat[grepl("band", names(dat))])) == 0
  dat=dat[na_rm, ]
}

#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @import rotl ape
#' @examples
#' @importFrom magrittr "%>%"
no.val.rm <- function(dat, no_value = 10000){
  library(dplyr)
  spectra <- dat[grepl("band", names(dat))] %>%
    dplyr::select(-one_of("band_site", "band_species"))
  dat <-  dat[!rowSums(spectra>no_value),]
  return(dat)
}

#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @import rotl ape
#' @examples
#' @importFrom magrittr "%>%"
optical.filter<- function(dat){
  ndvi <- (dat$band_90- dat$band_58)/(dat$band_58 + dat$band_90) <0.3
  nir860 <- (dat$band_96 + dat$band_97)/20000 < 0.1
  naval = as.logical(ndvi | nir860)
  dat[naval,] = NA
  return(dat)
}

mean.no.na <- function(x){mean(x, na.rm = T)}

#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @import rotl ape
#' @examples
#' @importFrom magrittr "%>%"
hiper_features <- function(dat, normalization = "norm2"){
    dat <- as.matrix(dat)
  if(normalization == "norm2"){
    normMat <- sqrt(apply(dat^2,FUN=sum,MAR=1, na.rm=TRUE)) 
    normMat <- matrix(data=rep(normMat,ncol(dat)),ncol=ncol(dat)) 
    normMat=dat/normMat
  }else if(normalization == "max_min"){
    min_x = apply(dat, 1, min)
    max_x <- apply(dat, 1, max)
    normMat = dat
    for(jj in 1: dim(normMat)[1]){
      #normMat[jj,] <- (normMat[jj,]-min_x[jj])/(normMat[jj,]-max_x[jj])
      normMat[jj,] <- (normMat[jj,]-min_x[jj])/(max_x[jj]-min_x[jj])
      
    }
    # normMat <- apply(dat,FUN=function(x, na.rm){(x-min_x)/(max_x-min_x)},MAR=1, na.rm=TRUE) 
    #normMat <- t(normMat)
  }else if(normalization == "max_min_bnd"){
    scaled <- as.data.frame(scale(dat, center = mins, scale = maxs - mins))
  }else{
    normMat = dat
  }
  return(normMat)
}

#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @import rotl ape
#' @examples
#' @importFrom magrittr "%>%"
get_mod_r2 <- function(pred, obs){
  #1 - sum((pred - obs)^2) / sum((obs - mean(obs, na.rm=T))^2)
  1 - sum((pred - obs)^2) / sum((obs - mean(obs, na.rm=T))^2)
} 
#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @import rotl ape
#' @examples
#' @importFrom magrittr "%>%"
category_to_dummy <- function(cat){
  #dummy variable example
  foo=list()
  foo$band_ <- cat
  mat <- model.matrix( ~ band_ - 1, data=foo)
  return(mat)
}

#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @import rotl ape
#' @examples
#' @importFrom magrittr "%>%"
#' 
plot_spectra<-function(dat){
  #plot reflectances
  plot_data <- dat %>% 
    dplyr::select(-one_of(c("site_ID", "species_ID",  "band_site","band_species", "flightpath"))) 
  plot_data <- plot_data[-1] %>%
    t %>%
    data.frame
  colnames(plot_data) = unlist(dat[1]) # the first row will be the header
  plot_data <- data.frame(bnd = 1:dim(plot_data)[1], plot_data)
  ggdat <- tidyr::gather(plot_data, treeID, Reflectance,-bnd)
  
  return(ggplot(ggdat, aes(x = bnd, y = Reflectance)) + 
           geom_line(aes(color = factor(treeID), alpha= 1), size = 0.2) +
           theme_bw()+
           theme(legend.position="none"))
  
}

