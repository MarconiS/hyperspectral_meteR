sf::st_crs(dat) = 4326
inpath="~/Documents/Data/Chapter3"
prdr = c("hsi")
for(x in 1:dim(dat)[1]){
  foo[[x]] <- get_dataset_from_point(dat[x,], prdr, inpath)
}
get_list_augmented_matrix = pbapply::pblapply(1:dim(dat)[1],
                                              function(x, prdr, inpath) get_dataset_from_point(dat[x,]),
                                              prdr = prdr, inpath= inpath)

neon_data = do.call(rbind.data.frame, get_list_augmented_matrix)

plts_needed <- dat$plotID %>% unique
plts_got <- list.files("~/Documents/Data/Chapter3/plots/hsi/")
plts_got <- substr(plts_got, 1, nchar(plts_got)-4)
plts_missing <- plts_needed[!plts_needed %in%plts_got]




#stack hiperspectral images
list_hsi <- list.files("~/Documents/Data/Chapter3/plots/hsi/", full.names = T)

good_pixs <-lapply(1:length(list_hsi), function(x) optical_filter(list_hsi[x]))
good_pixs <- do.call(rbind.data.frame, good_pixs)
#turn pixels non interesting out
optical_filter<- function(hsi_path){
  dat <- raster::brick(hsi_path)
  dat <- raster::as.data.frame(dat)
  colnames(dat) <- paste("band_", seq(1,ncol(dat)), sep="")
  ndvi <- (dat[[90]]- dat[[58]])/(dat[[58]] + dat[[90]]) <0.5
  nir860 <- (dat[[96]] + dat[97])/20000 < 0.3
  naval = as.logical(ndvi | nir860)
  dat[naval,] = NA
  return(dat)
}
#apply KLD and define bands to aggregate

write_csv(good_pixs, "./indir/all_hsi_for_kld.csv")

reticulate::source_python("../hsi_toolkit_py/dim_reduction/hdr.py")
aa<- getClusters(data.matrix(good_pics), 15)
write_csv(data.frame(aa), "./indir/bands_clusters_kld.csv")


