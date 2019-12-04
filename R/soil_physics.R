# in.statement <- soilDB::format_SQL_in_statement((s_data$mukey))
# 
# 
# q <- paste("SELECT
# component.mukey, taxclname,
#            taxorder, taxsuborder, taxgrtgroup, taxsubgrp
#            FROM legend
#            INNER JOIN mapunit ON mapunit.lkey = legend.lkey
# LEFT OUTER JOIN component ON component.mukey = mapunit.mukey
#            WHERE mapunit.mukey IN ", in.statement, sep="")
# # run query, process results, and return as data.frame object
# res <- soilDB::SDA_query(q)
# res <- res[complete.cases(res),]
# 

# query soil components by areasymbol and musym
#WHERE = "areasymbol = 'IN005' AND musym = 'MnpB2'"
get_soil_physics <- function(ii, s_data){
  soil_features= NULL
  for(ii in 1:nrow(s_data)){
    profile_data = soilDB::fetchSDA_component(WHERE = paste("areasymbol = '", 
                                                            s_data[ii,"areasymbol"] ,
                                                            "' AND musym = '", s_data[ii,"musym"],
                                                            "' AND mukey = '", s_data[ii,"mukey"],
                                                            "' AND nationalmusym = '", s_data[ii,"nationalmusym"],
                                                            "'", sep=""))
    soil_features[[ii]] = cbind(s_data$individualID[ii], profile_data@horizons)
  }                         
  soil_physics= do.call(rbind.data.frame, soil_features)
  colnames(soil_physics)[1]<- "individualID"
  soil_physics <- soil_physics %>% group_by(individualID) %>%
    summarize_if(is.numeric,  mean, na.rm = TRUE)
  
  soil_physics <- soil_physics[, colSums(is.na(soil_physics)) == 0]
  soil_physics %>% select(-one_of(c("hzID", "chkey", "hzname", "cokey")))
  prmsoil <-  FactoMineR::PCA(soil_physics[-1], ncp=10)
  plot(prmsoil)
  soil_features <- data.frame(soil_physics[["individualID"]], prmsoil$ind$coord)
  colnames(soil_features) <- c("individualID", paste("soil", 1:10, sep="."))
  return(soil_features)
}