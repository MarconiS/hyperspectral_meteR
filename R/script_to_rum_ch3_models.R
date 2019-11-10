library(tidyverse)
library(brms)
source("./R/get_phenetic_distance.R")
nComp = 10
invisible(base_data <- readr::read_csv("./indir/ch3.csv"))

traits_names <- c("nitrogenPercent",  "carbonPercent", "extractChlBConc",  "extractCarotConc"
                  , "extractChlAConc",  "ligninPercent",  "cellulosePercent"  
                  , "leafArea" , "leafMassPerArea",  "dryMassFraction")
random_effects <- c("domainID",  "nlcdClass") #taxonID, soilID, plantStatus
local_environment <- c("ope", "ect", "DTM")
plant_status <- c("CHM", "stemDiameter", "plantStatus")
reflectance <- colnames(base_data)[grep("hsi", colnames(base_data))]
bands_clusters <- readr::read_csv("./indir/bands_clusters_kld.csv")

# load climate trends
invisible(climate_trends <- readr::read_csv("./indir/climate_features.csv"))
colnames(bands_clusters) <- "KLDcluster"

reflectance_data <- base_data[colnames(base_data) %in% reflectance] %>%
  t %>%
  cbind.data.frame(bands_clusters) %>%
  group_by(KLDcluster) %>%
  summarise_all(list(mean)) %>%
  mutate_all(~ . /10000) %>%
  t %>% data.frame

reflectance_data <-  cbind.data.frame(base_data[["individualID"]], reflectance_data[-1,])
colnames(reflectance_data)[1] <- "individualID"

leaf_traits_data <- base_data[colnames(base_data) %in% traits_names] %>%
  mutate_all(log)# %>%
#scale

leaf_traits_data <-  cbind.data.frame(base_data[["individualID"]], leaf_traits_data)
colnames(leaf_traits_data)[1] <- "individualID"


climate_features <- climate_trends %>%
  dplyr::filter(individualID %in% base_data[["individualID"]])
prm <-  FactoMineR::PCA(climate_features[-1], ncp=10)
climate_features <- data.frame(climate_features[["individualID"]], prm$ind$coord)
colnames(climate_features)[1] <- "individualID"

soil_data = readRDS("./indir/soil_data.rds")
s_data <- do.call(rbind.data.frame, soil_data) %>%
  dplyr::select(individualID, areasymbol, spatialversion, 
                musym, nationalmusym, mukey, muareaacres, mupolygonkey, geometry)

local_soil <- c("musym")
soil_features <- s_data %>%dplyr::select(individualID, local_soil)
soil_features <- unique(soil_features)


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




full_dataset <- base_data %>% select(individualID, siteID) %>%
  ungroup() %>%
  dplyr::select(-one_of("taxonID")) %>%
  inner_join(leaf_traits_data) %>%
  inner_join(climate_features) %>%
  inner_join(site_features) %>%
  inner_join(random_effects_data) %>%
  inner_join(reflectance_data) %>% 
  inner_join(species_features) %>%
  inner_join(soil_features) %>%
  inner_join(status_features) %>%
  unique

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
  group_by(taxonID, siteID) %>%
  sample_frac(0.8)
test_data <- full_dataset %>% #dplyr::select(individualID, taxonID, siteID) %>%
  filter(!individualID %in% train_data[["individualID"]])


colnames(train_data)
train_data=train_data[complete.cases(train_data),]


#define formula and build the model
list_formulas <- paste(paste("mvbind("
                             , paste(colnames(leaf_traits_data)[-1], collapse = " , "), ") ~ ", sep = "")
                       #climate variables
                       , paste("s(",colnames(climate_features)[2:10],  ")", collapse = " + ", sep = "")  
                       , " + "
                       , paste(local_environment, collapse = " + ", sep = "") , " + "
                       , paste(c("CHM:plantStatus", "stemDiameter:plantStatus"), collapse = " + ", sep = "") , " + "
                       , "(1 | soil | musym)" , "+"
                       , "(1 | eco | domainID:nlcdClass)" , "+"
                       , "(1 | evo | taxonID"
                       , ")", collapse = "+")

test_data <- test_data[complete.cases(test_data),]
test_data <- filter(test_data, musym !="31D")

ind <- sapply(train_data, is.numeric)
train_data[ind] <- lapply(train_data[ind], scale)

taxaPD <- Matrix::nearPD( species_list$cov_taxa_eff)
fit <- brm(list_formulas, data = train_data, cores =2, chains = 2, iter = 2000,  
           family=brmsfamily("gaussian"), cov_ranef = list(taxonID = as.matrix(taxaPD$mat)))

bR2 <- bayes_R2(fit, newdata= test_data, robust = T)
print(bR2)

parameters_summary <- VarCorr(refl_fit)
obj = list(mod = refl_fit, R2 = bR2, params =parameters_summary )
saveRDS(obj, "./models/full_mod.rds")
