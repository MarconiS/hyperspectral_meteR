# structural equation modeling: path analysis for traits
library(brms)
get_paths_analysis <- function(train_dataset){

}

get_lavaan_sem <- function(train_dataset){
  library(lavaan)
  train_dataset <- train_dataset %>%
    left_join(leaf_traits_and_aop) %>% unique
  numerical_data <- train_dataset %>%
    select(-one_of(c("individualID", "domanID",
                     "plntStt", "nlcdCls", "taxonID", "musym"))) %>%
    scale
  fixed_effs <-  train_dataset %>%
    select((c( "nlcdCls", "plntStt", "taxonID", "musym")))
  fixed_effs$taxonID <-   substr(fixed_effs$taxonID,1,2) %>% factor
  mix_data <- cbind.data.frame(fixed_effs,  numerical_data)

  #min_dat <- table(mix_data$taxonID) > 34
  #levels(mix_data$taxonID)[min_dat]
  mix_data <- mix_data %>% filter(!nlcdCls=="mixedForest")
  taxaPD <- Matrix::nearPD()

  model <- '
    #latent variables definition
    climate =~ Dim.1 + Dim.2 + Dim.3 + Dim.4 + Dim.5 + Dim.6 + Dim.7 + Dim.8 + Dim.9 + Dim.10
    refl =~X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 +X9 + X10 + X11 + X12 + X13 + X14 + X15
    topo =~ DTM + ope + ect
    health =~ CHM + plntStt

    #responses
    pigms =~ extrCBC + extrcCC + extrCAC
    phys =~ lfMssPA + dryMssF + leafAre
    chem =~ ntrgnPr + crbnPrc
    strct =~ lgnnPrc + clllsPr

    #regressions
    phys ~ refl + climate + topo + health
    chem ~  refl + climate + topo + health
    strct ~ refl + climate + topo + health
    pigms ~ refl + climate + topo  + health
    refl =~ climate + topo + health
  '

  fit2 <- lavaan::sem(model
                     , data = mix_data
                     , group = c("nlcdCls")
                     , ordered=c("plntStt")
  )
  # Plot path diagram:
  semPaths(fit2, what = "est"
           , title = T
           , whatLabels = "name"
           , layout = "tree2"
           , curvePivot = TRUE)
}
