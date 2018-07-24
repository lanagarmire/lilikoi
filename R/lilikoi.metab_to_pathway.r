lilikoi.metab_to_pathway <- function(metaboliteNames, searchType) {
  #' A lilikoi.metab_to_pathway Function
  #'
  #' This function allows you to convert your metabolites id such as names, kegg ids, pubchem ids.
  #' into pathways. Metabolites which have no pathways will be excluded from any downstream analysis.
  #' @param metaboliteNames Your metabolite data
  #' @param searchType The type of the metabolites id such as 'name', 'kegg', 'hmdb','pubchem'
  #' @examples
  #' matches = lilikoi.metab_to_pathway(c("Asparagine", "Hypotaurine", "5-Oxoproline"), "name")
  #'
  #' matches = lilikoi.metab_to_pathway(c("C00152", "C00519", "C01879"), "kegg")
  #'
  #' @export

  if (searchType == "name") {
    matchIndexes <- .lilikoi.name_mapping_exact(metaboliteNames)
  } else {
    matchIndexes <- .lilikoi.id_mapping(metaboliteNames, searchType)
  }

  matchResults <- list(matchIndexes = matchIndexes, matchValues = lilikoi::data.metaboAnalyst$name[matchIndexes],
    matchStates = as.integer(!is.na(matchIndexes)))

  matchMatrix <- .lilikoi.get_match_matrix(matchResults, metaboliteNames)
  .lilikoi.print_results_summary(matchResults)
  .lilikoi.path_mapping(matchMatrix)
}

.lilikoi.print_results_summary <- function(results) {
  messages <- lilikoi._get_results_summary(results)
  print(messages$matched)
  print(messages$unmatched)
}


lilikoi._get_results_summary <- function(results) {
  nameHitLength <- length(results$matchIndexes)
  todoLength <- length(which(is.na(results$matchIndexes)))

  if (todoLength/nameHitLength > 0.5) {
    print("Over half of the compound IDs could not be matched to our database. Please make sure that correct compound IDs or common compound names are used.")
  } else if (todoLength > 100) {
    print("There are >100 compounds without matches. You can either proceed or if necessary, update these compound IDs and upload again.")
  }

  matched <- paste0((nameHitLength - todoLength), " ", "out of ", nameHitLength, " ", "Matched metabolites",
    " ", round((nameHitLength/(nameHitLength + todoLength)), 1) * 100, " ", "%")

  unmatched <- paste0(todoLength, " ", "out of ", nameHitLength, " ", "UnMatched metabolites", " ",
    round((todoLength/(nameHitLength + todoLength)), 1) * 100, " ", "%")
  paste0(matched, "\n", unmatched)

  list(matched = matched, unmatched = unmatched)
}

.lilikoi.name_mapping_exact <- function(metaboliteNames) {
  #' This function matches your metabolites names with databases
  #' with more than 20k metabolites. it allows you also to search among metabolite synonyms.
  #' This function was modified version of the name.match function in the below link:
  #' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R

  matchedNameIndices <- match(tolower(metaboliteNames), tolower(lilikoi::data.metaboAnalyst$name))
  matchedSynIndices <- match(tolower(metaboliteNames), tolower(lilikoi::data.metaboAnalyst$synonym))
  matchIndexes <- ifelse(is.na(matchedNameIndices), matchedSynIndices, matchedNameIndices)

  matchIndexes
}

.lilikoi.id_mapping <- function(metaboliteNames, sourceDatabase) {
  #' A .lilikoi.id_mapping Function
  #'
  #' This function allows you to exact matching of your metabolites hmdb, pubchem and kegg ids with data base
  #' consists of more than 20k metabolites.
  metaboliteNamesLower <- tolower(metaboliteNames)

  sourceDb <- lilikoi::data.metaboAnalyst[[paste0(sourceDatabase, "_id")]]
  if (is.null(sourceDb)) {
    matchIndexes <- match(metaboliteNamesLower, tolower(lilikoi::data.metaboAnalyst$hmdb))
    missHmdbIndexes <- is.na(matchIndexes)
    keggMatchIndexes <- match(metaboliteNamesLower, tolower(lilikoi::data.metaboAnalyst$kegg))
    # Give priority to HMDB matches
    matchIndexes[missHmdbIndexes] <- keggMatchIndexes[missHmdbIndexes]
  } else {
    matchIndexes <- match(metaboliteNamesLower, tolower(sourceDb))
  }

  matchIndexes
}

.lilikoi.get_match_matrix <- function(matchResults, metaboliteNames) {
  if (is.null(metaboliteNames)) {
    return(c())
  }
  metabolitesLength <- length(metaboliteNames)

  columns <- c("Query", "Match", "HMDB", "PubChem", "ChEBI", "KEGG", "METLIN", "FoundMatch")
  matchMatrix <- matrix("", nrow = metabolitesLength, ncol = length(columns))
  colnames(matchMatrix) <- columns

  for (i in 1:metabolitesLength) {
    hit <- lilikoi::data.metaboAnalyst[matchResults$matchIndexes[i], , drop = F]
    matchState <- matchResults$matchStates[i]
    matchMatrix[i, ] <- c(metaboliteNames[i], ifelse(matchState == 0, NA, matchResults$matchValues[i]),
      ifelse(matchState == 0, NA, hit$hmdb_id), ifelse(matchState == 0, NA, hit$pubchem_id), ifelse(matchState ==
        0, NA, hit$chebi_id), ifelse(matchState == 0, NA, hit$kegg_id), ifelse(matchState ==
        0, NA, hit$metlin_id), matchState)
  }

  matchMatrix
}

.lilikoi.path_mapping <- function(matchMatrix) {
  #' This function allows you to map user metabolites to the its corrsponding pathways
  #' @import dplyr
  hmdbPathways <- lilikoi::data.hmdb[match(matchMatrix[, "HMDB"], lilikoi::data.hmdb$Accession), ,
    drop = F]

  # We use arrange_ for SE vs NSE so R CMD Check passes
  table_result <- data.frame(matchMatrix) %>% mutate(pathway = hmdbPathways$Pathway_Name) %>% arrange_("pathway")

  table_result
}
