data <- readr::read_csv("~/Documents/Data/reflectance_vgstr.csv")
data <- readr::read_csv("~/Documents/GitHub/redSquirrel/Dat/JERC/cat_reget_ccbid.csv")


optical.filter <- function(dat){
  ndvi <- (dat[[91]]- dat[[59]])/(dat[[59]] + dat[[91]]) > 0.2
  nir860 <- (dat[[97]] + dat[[98]])/20000 > 0.1
  naval = as.logical(ndvi & nir860)
  dat = dat[naval,]
  return(dat)
}


data[-1] <- replace(data[-1], data[-1] < 0, NA)

readr::write_csv(dat, "~/Documents/GitHub/redSquirrel/Dat/JERC/cat_reget_ccbid.csv")

ancillary <- read_csv("~/Documents/GitHub/NeonTreeEvaluation/field_data.csv")
ancillary <- ancillary %>% filter(canopyPosition %in% c("Full sun",  "Partially shaded", "Open grown"))
ancillary <- ancillary %>% select(individualID, taxonID, scientificName, siteID)
dat = data %>% filter(individulaID %in% ancillary$individualID)
dat = optical.filter(data)
dat[-1] <- hiper_features(dat[-1], normalization = "norm2")
colnames(dat)[1] <- "individualID"
dat[-1] <- scale(dat[-1])
source_python("../hsi_toolkit_py/dim_reduction/hdr.py")
reflectance_data <- as.matrix(dat[-1])
bands = 10
KLD <- getClusters(reflectance_data, numBands = bands)
KLD

library(reticulate)

full_dataset <- left_join(cc3, ancillary, by="individualID")
full_dataset <- unique(full_dataset)
full_dataset$individualID <- factor(full_dataset$individualID)
full_dataset$individualID <- as.numeric(full_dataset$individualID)
dat$individualID %>% unique %>% length
write_csv(full_dataset, "../ccb-id/ccbid/support_files/withNA2019.csv")

