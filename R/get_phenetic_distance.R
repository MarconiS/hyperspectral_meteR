#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @import rotl ape
#' @examples
#' @importFrom magrittr "%>%"
get_cophenetic_distance <- function(species_list, sc_name = "scientificName"){
  

  # Get scientific names
  species_list <- unique(species_list)
  species_list <- species_list[!is.na(species_list[[sc_name]]),]
  species_list[[sc_name]] <- word(species_list[[sc_name]], 1,2)
  species_list[[sc_name]] <- str_replace(species_list[[sc_name]], " spp.", "")
  species_list[[sc_name]] <- str_replace(species_list[[sc_name]], " sp.", "")
  taxa <- unique(species_list[[sc_name]])
  taxa <- taxa[!is.na(taxa)]
  taxa[taxa == "Rubus hawaiensis"] <- "Rubus"
  taxa[taxa == "Ceanothus cuneatus"] <- "Ceanothus"
  taxa[taxa == "Rubus armeniacus"] <- "Rubus"
  taxa[taxa == "Colubrina elliptica"] <- "Colubrina"
  taxa[taxa == "Colubrina arborescens"] <- "Colubrina"
  taxa[taxa == "Rubus oklahomus"] <- "Rubus"
  taxa[taxa == "Rubus laciniatus"] <- "Rubus"
  taxa[taxa == "Rubus leucodermis"] <- "Rubus"
  taxa[taxa == "Serenoa repens"] <- NA
  taxa[taxa == "Quercus ×heterophylla"] <-" Quercus phellos"
  taxa <- unique(taxa)
  taxa <- taxa[!is.na(taxa)]

  if(length(taxa)/250<=1){
    resolved_names <- rotl::tnrs_match_names(taxa[1:length(taxa)])
  }else{
    resolved_names <- rotl::tnrs_match_names(taxa[1:250])
    for(n in 1:floor(length(taxa)/250)){
      tmp <- rotl::tnrs_match_names(taxa[(1+n*250):min((n+1)*250, length(taxa))])
      resolved_names = rbind.data.frame(resolved_names, tmp)
    }
  }
  species_list[[sc_name]] <- tolower(species_list[[sc_name]])
  resolved_names = resolved_names[!is.na(resolved_names$approximate_match),]
  # format names to be linked to correct taxonIDs
  resolved_names <- dplyr::left_join(resolved_names, species_list,
                                     by = c("search_string" = sc_name))
  
  resolved_names <- dplyr::group_by(resolved_names, ott_id) %>% dplyr::slice(1)
  resolved_names = resolved_names[!is.na(resolved_names$ott_id),]


  normalize <- function(x) {
    return ((max(x) - x) / max(x))
  }
  resolved_names <- resolved_names[!is.na(resolved_names$ott_id),]
  #get tree
  my_tree <-  rotl::tol_induced_subtree(ott_ids = resolved_names$ott_id)
  #my_tree <- ape::compute.brlen(my_tree, 1)
  
  patristic_distMatrix<-ape::cophenetic.phylo(my_tree)
  #get cophenetic distance
  cross_species_phylogenetic_correlation = normalize(patristic_distMatrix)
  cross_species_phylogenetic_correlation =
    cross_species_phylogenetic_correlation[, !colnames(cross_species_phylogenetic_correlation) 
                                           %in% c("mrcaott31448ott33231", "Cladrastis_clade_ott7055956")]
  cross_species_phylogenetic_correlation =
    cross_species_phylogenetic_correlation[!rownames(cross_species_phylogenetic_correlation) 
                                           %in% c("mrcaott31448ott33231", "Cladrastis_clade_ott7055956"), ]
  #get taxonID from ott_ids
  get_ott_id = colnames(cross_species_phylogenetic_correlation)
  get_ott_id = sapply(strsplit(get_ott_id, "_ott"), "[", 2) %>%
    as.integer %>%
    data.frame
  taxa_id = left_join(get_ott_id, resolved_names,
                      by = c("." = "ott_id")) %>%select(taxonID)
  
  taxa_id <- taxa_id[!is.na(taxa_id)]
  # set maxtrix col and row names following taxonID
  colnames(cross_species_phylogenetic_correlation) =
    rownames(cross_species_phylogenetic_correlation) = taxa_id

  return(list(phylotree = my_tree,
              cov_taxa_eff = cross_species_phylogenetic_correlation,
              resolved_names = resolved_names))

}
