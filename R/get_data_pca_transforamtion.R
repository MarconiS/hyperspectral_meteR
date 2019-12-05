get_data_pca_transforamtion <- function(climate_trends) {
  #PCA by variable
  var <- c("daylength", "prec", "rad", "snow_melt", "tmax", "tmin", "vp")
  climate_pca = NULL
  climate_features = NULL
  for(ii in var){
    climate_pca[[ii]] <- climate_trends %>%
      select(contains(ii)) %>%
      FactoMineR::PCA(ncp=1)
    #get pca transformed data
    climate_features[[ii]] = climate_pca[[ii]]$ind$coord
  }
  climate_features = do.call(cbind.data.frame, climate_features)
  colnames(climate_features) <- var
  climate_features <- data.frame(climate_trends[["individualID"]], climate_features)
  colnames(climate_features)[1] <- "individualID"
  return(list(climate_pca=climate_pca, climate_features=climate_features))
}