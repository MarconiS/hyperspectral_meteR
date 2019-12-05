#load libraries and functions first
library(tidyverse)
library(brms)
source("./R/get_phenetic_distance.R")
source("./R/utilities.R")
source("./R/soil_physics.R")
nComp = 20
refl_dr = "kld"
base_data <- readr::read_csv("./indir/ch3.csv")
reflectance <- colnames(base_data)[grep("hsi", colnames(base_data))]

set.seed(1987)
train_data <- base_data %>% 
  dplyr::select(individualID, taxonID, siteID) %>%
  group_by(taxonID, siteID) %>%
  sample_frac(0.8)

traits_names <- c("nitrogenPercent",  "carbonPercent", "extractChlBConc",  "extractCarotConc"
                  , "extractChlAConc",  "ligninPercent",  "cellulosePercent"  
                  , "leafArea" , "leafMassPerArea",  "dryMassFraction")
random_effects <- c("domainID",  "nlcdClass") #taxonID, soilID, plantStatus
local_environment <- c("ope", "ect", "DTM")
plant_status <- c("CHM", "stemDiameter", "plantStatus")
bands_clusters <- readr::read_csv("./indir/bands_20clusters_kld.csv")

# load climate trends
invisible(climate_trends <- readr::read_csv("./indir/climate_features.csv"))
colnames(bands_clusters) <- "KLDcluster"


leaf_traits_data <- base_data[colnames(base_data) %in% traits_names] %>%
  mutate_all(log)# %>%
#scale

leaf_traits_data <-  cbind.data.frame(base_data[["individualID"]], leaf_traits_data)
colnames(leaf_traits_data)[1] <- "individualID"

reflectance_data <- base_data[colnames(base_data) %in% reflectance] 

if(refl_dr == "kld"){
  reflectance_data[reflectance_data<0] <- NA
  reflectance_data = hiper_features(reflectance_data, normalization = "norm2")
  reflectance_data <- reflectance_data %>%
    t %>%
    cbind.data.frame(bands_clusters) %>%
    group_by(KLDcluster) %>%
    summarise_all(list(mean)) %>%
    t %>% data.frame
  
  reflectance_data <-  cbind.data.frame(base_data[["individualID"]], reflectance_data[-1,])
  reflectance_data <-reflectance_data[complete.cases(reflectance_data),]
}else if(refl_dr == "pls"){
  library(plsdepot)
  pls_transform <- base_data %>% filter(individualID %in% train_data[["individualID"]])
  pls_transform <- pls_transform[colnames(pls_transform) %in% reflectance] 
  pls_transform <- hiper_features(pls_transform, normalization = "norm2")
  pls_transform[pls_transform<0] <- NA
  pls_responses <- leaf_traits_data %>% filter(individualID %in% train_data[["individualID"]])
  pls_responses <- pls_responses[complete.cases(pls_transform), ]
  pls_transform <- pls_transform[complete.cases(pls_transform), ]
  
  #get center vars
  x_mean <- apply(pls_transform, 2, mean)
  x_sd <-  apply(pls_transform, 2, sd)
  
  x_weights <- plsreg2(pls_transform, pls_responses[-1],comps = nComp)
  x_scores <- x_weights$x.scores
  colnames(x_scores) <- paste("X_", 1:nComp, sep="")
  #x_scores %>% as.data.frame %>% dplyr::select(-one_of("individualID")) %>% plot_spectra
  
  #transform all data accordingly
  reflectance_data <- hiper_features(reflectance_data, normalization = "norm2") %>%
    scale(center = x_mean, scale = x_sd)
  
  x_rotation <- solve(t(x_weights$x.loads) %*% x_weights$mod.wgs)
  x_rotation <- x_weights$mod.wgs %*% x_rotation
  reflectance_data <- reflectance_data %*% x_rotation  %>%
    data.frame
  #reflectance_data %>% as.data.frame %>% dplyr::select(-one_of("individualID")) %>% plot_spectra
  reflectance_data <-  cbind.data.frame(base_data[["individualID"]], reflectance_data)
  
}
colnames(reflectance_data)[1] <- "individualID"
# check_outliers <- apply(reflectance_data[-1], 1, function(x)sum(abs(x)))
# check_outliers <- cbind.data.frame(reflectance_data[["individualID"]], check_outliers) 
# check_outliers <-  filter(check_outliers, check_outliers > 300)
# 
# climate_features <- climate_trends %>%
#   dplyr::filter(individualID %in% base_data[["individualID"]])
# prm <-  FactoMineR::PCA(climate_features[-1], ncp=10)
# climate_features <- data.frame(climate_features[["individualID"]], prm$ind$coord)
# colnames(climate_features)[1] <- "individualID"
climate <- get_data_pca_transforamtion(climate_trends)
climate_features = climate$climate_features[-1]
colnames(climate_features) <- c("individualID","daylength", "prec", "rad", "snow_melt", "tmax", "tmin", "vp")
# 
# soil_data = readRDS("./indir/soil_data2.rds")
# s_data <- do.call(rbind.data.frame, soil_data) %>%
#   dplyr::select(individualID, areasymbol, 
#                 musym, nationalmusym, mukey, geometry) %>% unique

# local_soil <- c("musym")
# soil_features <- s_data %>%dplyr::select(individualID, local_soil)
# soil_features <- unique(soil_features)
#soil_features <- get_soil_physics(s_data)
soil_features <- readr::read_csv("./indir/soil_features.csv")
species_list <- base_data %>%
  dplyr::select(taxonID, scientificName)
sp_list2 <- readr::read_csv("./indir/vegetation_structure.csv") %>%
  dplyr::select(taxonID, scientificName) 

species_list <- rbind.data.frame(species_list, sp_list2) %>%
  unique 
species_list <-  get_cophenetic_distance(species_list, sc_name = "scientificName")

species_features <- base_data %>% 
  dplyr::select(individualID, scientificName)
species_features[["scientificName"]] <- word(species_features[["scientificName"]], 1,2)
species_features[["scientificName"]] <- 
  str_replace(species_features[["scientificName"]], " spp.", "")
species_features[["scientificName"]] <- 
  str_replace(species_features[["scientificName"]], " sp.", "")
species_features[["scientificName"]] <- tolower(species_features[["scientificName"]])
species_features <-left_join(species_features, species_list$resolved_names, 
                             by = c("scientificName" = "search_string")) %>%
  select(individualID, taxonID)
colnames(species_features) <- c("individualID", "taxonID")


site_features <- base_data[colnames(base_data) %in% local_environment] 
site_features <-  cbind.data.frame(base_data[["individualID"]], site_features)
colnames(site_features)[1] <- "individualID"

status_features <- base_data[colnames(base_data) %in% plant_status] 
status_features <-  cbind.data.frame(base_data[["individualID"]], status_features)
status_features$plantStatus[status_features$plantStatus!="OK"]="damaged"
colnames(status_features)[1] <- "individualID"


random_effects_data <- base_data[colnames(base_data) %in% random_effects] 
random_effects_data <-  cbind.data.frame(base_data[["individualID"]], random_effects_data)
colnames(random_effects_data)[1] <- "individualID"




full_dataset <- base_data %>% select(individualID, siteID, taxonID) %>%
  ungroup() %>%
  dplyr::select(-one_of("taxonID")) %>%
  inner_join(leaf_traits_data) %>%
  inner_join(climate_features) %>%
  inner_join(site_features) %>%
  inner_join(random_effects_data) %>%
  inner_join(species_features) %>%
  inner_join(status_features) %>%
  inner_join(reflectance_data) %>% 
  inner_join(soil_features) %>%
  unique

full_dataset <- full_dataset %>%
  filter(!individualID %in% as.character(check_outliers[[1]]))
full_dataset<- full_dataset[!is.na(full_dataset[["taxonID"]]),]
full_dataset[!full_dataset[["taxonID"]] 
             %in% colnames(species_list$cov_taxa_eff),"taxonID"] %>% 
  unique %>% 
  print
full_dataset[full_dataset[["taxonID"]]=="ABIES","taxonID"] = "ABBA"
full_dataset[full_dataset[["taxonID"]]=="ACSA3","taxonID"] = "ACSA2"
full_dataset[full_dataset[["taxonID"]]=="BETUL","taxonID"] = "BELE"
full_dataset[full_dataset[["taxonID"]]=="BOURR","taxonID"] = "BOSU2"
full_dataset[full_dataset[["taxonID"]]=="FRAXI","taxonID"] = "FRAM2"
full_dataset[full_dataset[["taxonID"]]=="HALES","taxonID"] = "HACA3"
full_dataset[full_dataset[["taxonID"]]=="NYSY","taxonID"] = "NYBI"
full_dataset[full_dataset[["taxonID"]]=="OXYDE","taxonID"] = "OXAR"
full_dataset[full_dataset[["taxonID"]]=="PINUS","taxonID"] = "PIEC2"
full_dataset[full_dataset[["taxonID"]]=="SASSA","taxonID"] = "SAAL5"

full_dataset = full_dataset %>% filter(siteID != "JORN")


set.seed(1987)
train_data <- full_dataset %>% #dplyr::select(individualID, taxonID, siteID) %>%
  filter(individualID %in% train_data[["individualID"]])
# group_by(taxonID, siteID) %>%
# sample_frac(0.8)

test_data <- full_dataset %>% #dplyr::select(individualID, taxonID, siteID) %>%
  filter(!individualID %in% train_data[["individualID"]])


colnames(train_data)
train_data=train_data[complete.cases(train_data),]
#train_data=train_data %>% filter(!individualID %in% check_outliers[1])


#define formula and build the model
list_formulas <- paste(paste("mvbind("
                      , paste(colnames(leaf_traits_data)[-1], collapse = " , "), ") ~ ", sep = "")
                       #climate variables
                       , paste("s(",colnames(climate_features),  ")", collapse = " + ", sep = "")  
                       , " + "
                       , paste("s(",colnames(soil_features)[2:6],  ")", collapse = " + ", sep = "")  
                       , " + "
                       , paste(local_environment, collapse = " + ", sep = "") , " + "
                       , paste(c("CHM", "stemDiameter"), collapse = " + ", sep = "") , " + "
                       , "(1 | ind | plantStatus)" , "+"
                       , "(1 | eco | domainID:nlcdClass)" , "+"
                       #, "(1 | evo | taxonID"
                       #, ")"
                      , collapse = "+")

test_data <- test_data[complete.cases(test_data),]
#test_data <- filter(test_data, musym !="31D")

ind <- sapply(train_data, is.numeric)
train_data[ind] <- lapply(train_data[ind], scale)

taxaPD <- Matrix::nearPD( species_list$cov_taxa_eff)
fit <- brm(list_formulas, data = train_data, cores =2, chains = 2, iter = 2000,  seed = 1987,
           family=brmsfamily("gaussian"), cov_ranef = list(taxonID = as.matrix(taxaPD$mat)))
cls=1:ncol(test_data)
tst_scaled <- lapply(cls[as.vector(ind)], function(x){
  scale(test_data[x], center = attr(train_data[[x]],"scaled:center"),
        scale = attr(train_data[[x]],"scaled:scale"))})
tst_scaled <- do.call(cbind.data.frame, tst_scaled)
tst_scaled <- cbind.data.frame(test_data[!ind], tst_scaled)
bR2 <- bayes_R2(fit, newdata= tst_scaled, robust = T)

print(bR2)

parameters_summary <- VarCorr(fit)
obj = list(mod = fit, R2 = bR2, params =parameters_summary )
saveRDS(obj, "./models/ind_env_mod.rds")
# 
# me_regional <- predict(fit, effects = colnames(climate_features)[-1])
# me_local <- predict(fit, newdata = tst_no_cilm)
# 
# formula_reflectance <- paste(paste("mvbind("
#                                    , paste(colnames(leaf_traits_data)[-1], collapse = " , "), ") ~ ", sep = "")
#                              #climate variables
#                              , paste("s(",paste("X", 1:nComp, sep=""),  ")", collapse = " + ", sep = "")  
#                              , " + "
#                              , "(1 | site | siteID)", collapse = "+")
# 
# 
# fit_refl2 <- brm(formula_reflectance, data = train_data, cores =4, 
#                 chains = 4, iter = 4000,  seed = 1987,
#                 family=brmsfamily("gaussian"), prior = prior(horseshoe()))
# 
# test_data <- test_data[complete.cases(test_data),]
# cls=1:ncol(test_data)
# tst_scaled <- lapply(cls[as.vector(ind)], function(x){
#   scale(test_data[x], center = attr(train_data[[x]],"scaled:center"),
#         scale = attr(train_data[[x]],"scaled:scale"))})
# tst_scaled <- do.call(cbind.data.frame, tst_scaled)
# tst_scaled <- cbind.data.frame(test_data[!ind], tst_scaled)
# 
# bR2_refl2 <- bayes_R2(fit_refl, newdata= tst_scaled, robust = T)
# bR2_refl2
# objref = list(mod = fit_refl, R2 = bR2_refl)#, params =parameters_summary )
# saveRDS(objref, "./models/refl50.rds")
