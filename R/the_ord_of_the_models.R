the_one <- function(data){
  library(brms)
  library(tidyverse)
  traits_names <- c("nitrogenPercent"#,  "leafMassPerArea"#, "carbonPercent"
                    #"extractChlBConc",  "extractCarotConc" , "extractChlAConc"
                   #, "ligninPercent",  "cellulosePercent"
                   , "leafArea" #, "dryMassFraction"
                   )
  hc=5
  sc=5
  climate_features<-"tmin"
  train_data <- data %>% group_by(siteID, taxonID)%>%
    sample_frac(0.8)
  test_data <- data%>%filter(!individualID %in% train_data[["individualID"]])
  traits_formula <- paste(paste("mvbind("
        , paste(c(traits_names), collapse = " , "), ") ~ ", sep = "")
       # , paste("s(",colnames(climate_features)[-1],  ")", collapse = " + ", sep = "")
        #, " + "
        , paste("s(",paste("band_", 1:5, sep=""),  ")", collapse = " + ", sep = "")
        , " + "
        , paste("s(",paste("soil", 1:sc, sep=""),  ")", collapse = " + ", sep = "")
        , " + "
        , paste("s(",c("ope", "ect", "DTM"),")", collapse = " + ", sep = "") , " + "
        ,"s(CHM)"
        , "+"
        , "(1 | eco | nlcdClass)"
        , "+"
        , "(1 | scale | siteID:domainID)")
  species_formula <- paste(paste("mvbind("
                              , paste("taxonID", collapse = " , "), ") ~ ", sep = "")
                        # , paste("s(",colnames(climate_features)[-1],  ")", collapse = " + ", sep = "")
                        #, " + "
                        , paste("s(",paste("band_", 1:5, sep=""),  ")", collapse = " + ", sep = "")
                        , " + "
                        #, paste("s(",paste("soil", 1:sc, sep=""),  ")", collapse = " + ", sep = "")
                        , " + "
                        , paste("s(",c("ope", "ect", "DTM"),")", collapse = " + ", sep = "") , " + "
                        ,"s(CHM)"
                        , "+"
                        , "(1 | eco | nlcdClass)"
                        , "+"
                        , "(1 | scale | siteID:domainID)")
  health_formula <- paste(paste("mvbind("
                                , paste("plantStatus", collapse = " , "), ") ~ ", sep = "")
                          # , paste("s(",colnames(climate_features)[-1],  ")", collapse = " + ", sep = "")
                          #, " + "
                          , paste("s(",paste("band_", 1:5, sep=""),  ")", collapse = " + ", sep = "")
                          , " + "
                          #, paste("s(",paste("soil", 1:sc, sep=""),  ")", collapse = " + ", sep = "")
                          , " + "
                          , paste("s(",c("ope", "ect", "DTM"),")", collapse = " + ", sep = "") , " + "
                          ,"s(CHM)"
                          , "+"
                          , "(1 | eco | nlcdClass)"
                          , "+"
                          , "(1 | scale | siteID:domainID)")
  bf_traits <- bf(traits_formula, family = "gaussian")
  bf_health <- bf(health_formula, family = "binomial")
  bf_species <- bf(species_formula, family = "categorical")
  fit <- brm(bf_traits+bf_health+bf_species, data = train_data, cores =2, chains = 2, iter = 40,  seed = 1987)             #, cov_ranef = list(taxonID = as.matrix(taxaPD$mat))

}
