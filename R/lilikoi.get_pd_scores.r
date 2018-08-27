lilikoi.get_pd_scores <- function(metaboliteMeasurements, metabolitePathwayTable) {
  #' A lilikoi.get_pd_scores Function
  #'
  #' This function allows you to compute Pathway Deregulation Score deriving
  #' make sure that you have the below database for the metabolites and pathway list:
  #' meta_path.RData
  #'
  #' @param metabolitePathwayTable This is the metabolitePathwayTable from lilikoi.metab_to_pathway function. This table includes the metabolite ids and their corresponding HMDB ids
  #' @keywords PDS
  #' @export
  #' @import pathifier dplyr
  #' @examples PDSmatrix= lilikoi.get_pd_scores(metaboliteMeasurements, metabolitePathwayTable)
  #' lilikoi.get_pd_scores(metaboliteMeasurements, metabolitePathwayTable)
  
  metaboliteIds <- metabolitePathwayTable %>% filter(pathway != "NA") %>% select(Query, HMDB)
  
  # Keep just the columns which we have matched
  metaboliteLevelsWithMatchingQuery <- metaboliteMeasurements[, t(metaboliteIds["Query"])]
  
  # Now swap those names with the HMDB id
  colnames(metaboliteLevelsWithMatchingQuery) <- t(metaboliteIds["HMDB"])
  
  # Finally transpose that set
  metaboliteLevelsWithHmdbIds <- t(metaboliteLevelsWithMatchingQuery)
  
  # uses mapping dbs
  isHealthySample <- metaboliteMeasurements$Label %>% as.numeric %>% -1 %>% as.logical
  
  ## although pathifier is built for gene inputs, turns out it works for metabolites as well
  pds <- pathifier::quantify_pathways_deregulation(as.matrix(metaboliteLevelsWithHmdbIds), row.names(metaboliteLevelsWithHmdbIds), 
                                                   lilikoi.metabolites, lilikoi.pathways, normals = isHealthySample, attempts = 5, min_exp = 0, 
                                                   min_std = 0)
  
  pdsMatrix <- matrix(as.data.frame(pds$scores), nrow = length(names(pds$scores)), byrow = TRUE)
  
  colnames(pdsMatrix) <- colnames(metaboliteLevelsWithHmdbIds)
  rownames(pdsMatrix) <- names(pds$scores)
  mode(pdsMatrix) <- "numeric"
  
  pdsMatrix
}

