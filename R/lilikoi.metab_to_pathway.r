lilikoi.metab_to_pathway <- function(metaboliteNames, searchType, hmdb = T, pubchem = T, chebi = F, kegg = T, 
                                     metlin = F) {
  #' A lilikoi.metab_to_pathway Function
  #'
  #' This function allows you to convert your metabolites id such as names, kegg ids, pubchem ids.
  #' into pathways. Metabolites which have no pathways will be excluded from any downstream analysis.
  #' This function was modified version of the name.match function in the below link:
  #' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R
  #' @param metaboliteNames Your metabolite data
  #' @param searchType The type of the metabolites id such as 'name', 'kegg', 'hmdb','pubchem'
  #' @keywords Match
  #' @export
  #' @examples metabolitePathwayTable=lilikoi.metab_to_pathway('name')
  #' lilikoi.metab_to_pathway()
  #'
  #'
  
  searchOptionsHash <- hash()
  # record the filter for 8 major databases
  searchOptionsHash$returnColumns <- c(hmdb, pubchem, chebi, kegg , metlin)
  
  # distribute job
  searchOptionsHash$searchType <- searchType
  print(searchOptionsHash)
  if (searchType == "name") {
    results <- .lilikoi.name_mapping_exact(metaboliteNames, searchOptionsHash)
  } else {
    results <- .lilikoi.id_mapping(metaboliteNames, searchOptionsHash)
  }
  
  .lilikoi.print_results_summary(results$stats$matchIndexes)
  
  results$paths
}

.lilikoi.print_results_summary <- function(matchIndexes) {
  #' do some sanity check
  #' A .lilikoi.print_results_summary Function
  #'
  #' This function prints how many metabolites matched and unmatched from Metatopathway function
  #' @param name.map is the dataframe of metabolites and its corresponding ids.
  #' @param hit.inx is the status of the metabolite hits (number if the metabolites matched the database or NA otherwise.)
  #' @keywords santiy search
  #' @export
  #' @examples .lilikoi.print_results_summary(name.map,hit.inx)
  #' .lilikoi.print_results_summary
  #'
  #'
  
  messages = lilikoi._get_results_summary(length(matchIndexes), length(which(is.na(matchIndexes))))
  print(messages$matched)
  print(messages$unmatched)
}


lilikoi._get_results_summary = function (matchCount, missCount) {
  
  if (missCount / matchCount > 0.5) {
    return("Over half of the compound IDs could not be matched to our database. Please make sure that correct compound IDs or common compound names are used.")
  } else if (missCount > 100) {
    return("There are >100 compounds without matches. You can either proceed or if necessary, update these compound IDs and upload again.")
  }
  
  matched = paste0((matchCount - missCount), " ", "out of ", matchCount, " ", "Matched metabolites", " ",
                   round((matchCount / (matchCount + missCount)), 1) * 100, " ", "%")
  
  unmatched = paste0( missCount, " ", "out of ", matchCount, " ", "UnMatched metabolites", " ",
                      round((missCount / (matchCount + missCount)), 1) * 100, " ", "%")
  paste0(matched, "\n", unmatched)
  
  list(
    matched = matched,
    unmatched = unmatched
  )
}

.lilikoi.name_mapping_exact <- function(metaboliteNames, searchOptionsHash) {
  #' A .lilikoi.name_mapping_exact Function
  #'
  #' This function matches your metabolites names with databases
  #' with more than 20k metabolites. it allows you also to search among metabolite synonyms.
  #' Make sure that you have the below two databases for the exact matching:
  #' cmpd_db.rda, and syn_nms_db.rda. These should be loaded from lilikoi/R/sysdata.rda
  #' This function was modified version of the name.match function in the below link:
  #' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R
  #' @param metaboliteNames Dataset
  #' @keywords Exact Matching
  #' @export
  #' @examples .lilikoi.name_mapping_exact()
  #' .lilikoi.name_mapping_exact()
  #'
  #'
  
  metabolitesLength <- length(metaboliteNames)
  
  # first find exact match to the common compound names
  matchIndexes <- match(tolower(metaboliteNames), tolower(lilikoi.cmpd$name))
  matchValues <- lilikoi.cmpd$name[matchIndexes]
  matchStates <- as.integer(!is.na(matchIndexes))
  
  syns.list <- lilikoi.syn$syns.list
  
  indexesWithNoMatch <- which(is.na(matchIndexes))
  if (length(indexesWithNoMatch) > 0) {
    for (i in 1:length(syns.list)) {
      syns <- syns.list[[i]]
      hitInx <- match(tolower(metaboliteNames[indexesWithNoMatch]), tolower(syns))
      
      hitPos <- which(!is.na(hitInx))
      if (length(hitPos) > 0) {
        orig.inx <- indexesWithNoMatch[hitPos]
        matchIndexes[orig.inx] <- i
        matchValues[orig.inx] <- lilikoi.cmpd$name[i]  # show common name
        matchStates[orig.inx] <- 1
        
        # update unmatched list
        indexesWithNoMatch <- indexesWithNoMatch[is.na(hitInx)]
      }
      if (length(indexesWithNoMatch) == 0) 
        break
    }
  }
  
  matchResults <- list(matchIndexes = matchIndexes, matchValues = matchValues, matchStates = matchStates)
  
  csv_result <- .lilikoi.get_mapping_result_table(matchResults, metaboliteNames, searchOptionsHash)
  .lilikoi.store_result(csv_result)
  
  list(paths = .lilikoi.path_mapping(csv_result), stats = matchResults)
}

.lilikoi.id_mapping <- function(metaboliteNames, searchOptionsHash) {
  #' A .lilikoi.id_mapping Function
  #'
  #' This function allows you to exact matching of your metabolites hmdb, pubchem and kegg ids with data base
  #' consists of more than 20k metabolites.
  #' make sure that you have the below database for the exact matching:
  #' cmpd_db.rda
  #' This function was modified version of the name.match function in the below link:
  #' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R
  #' @param sourceDatabase which is one of kegg id, pubchem is, hmdb id
  #' @keywords Exact Matching
  #' @export
  #' @examples .lilikoi.id_mapping(sourceDatabase)
  #' .lilikoi.id_mapping()
  #'
  
  sourceDatabase <- searchOptionsHash$searchType
  
  metaboliteNamesLower <- tolower(metaboliteNames)
  
  # for now
  sourceDb <- lilikoi.cmpd[[paste0(sourceDatabase, "_id")]]
  if (!is.null(sourceDb)) {
    matchIndexes <- match(metaboliteNamesLower, tolower(sourceDb))
  } else {
    matchIndexes <- match(metaboliteNamesLower, tolower(lilikoi.cmpd$hmdb))
    missHmdbIndexes <- is.na(matchIndexes)
    keggMatchIndexes <- match(metaboliteNamesLower, tolower(lilikoi.cmpd$kegg))
    matchIndexes[missHmdbIndexes] <- keggMatchIndexes[missHmdbIndexes]
  }
  
  matchResults <- list(matchIndexes = matchIndexes, matchValues = lilikoi.cmpd$name[matchIndexes], 
                       matchStates = as.integer(!is.na(matchIndexes)))
  
  csv_result <- .lilikoi.get_mapping_result_table(matchResults, metaboliteNames, searchOptionsHash)
  .lilikoi.store_result(csv_result)
  
  list(paths = .lilikoi.path_mapping(csv_result), stats = matchResults)
}

.lilikoi.store_result <- function(csv_result) {
  if (length(csv_result)) {
    # store the value for report
    write.csv(csv_result, file = "name_map.csv", row.names = F)
  }
}

.lilikoi.get_mapping_result_table <- function(matchResults, metaboliteNames, searchOptionsHash) {
  #' A .lilikoi.get_mapping_result_table Function
  #'
  #' This function allows you to generate csv file with user metabolites matched names, kegg, pubchem
  #' and hmdb ids
  #' This function was modified version of the name.match function in the below link:
  #' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R
  #' @param nonething
  #' @keywords data.map
  #' @export
  #' @examples .lilikoi.get_mapping_result_table()
  #' .lilikoi.get_mapping_result_table()
  #'
  
  if (is.null(metaboliteNames)) {
    return(c())
  }
  metabolitesLength <- length(metaboliteNames)
  
  # style for highlighted background for unmatched names
  pre.style <- NULL
  post.style <- NULL
  
  # style for no matches
  if (searchOptionsHash$searchType == "name") {
    no.prestyle <- "<strong style=\"background-color:yellow; font-size=125%; color=\"black\">"
    no.poststyle <- "</strong>"
  } else {
    no.prestyle <- "<strong style=\"background-color:red; font-size=125%; color=\"black\">"
    no.poststyle <- "</strong>"
  }
  
  matchIndexes <- matchResults$matchIndexes
  matchValues <- matchResults$matchValues
  matchStates <- matchResults$matchStates
  
  # contruct the result table with cells wrapped in html tags the unmatched will be highlighted in
  # different background
  html.res <- matrix("", nrow = metabolitesLength, ncol = 8)
  csv.res <- matrix("", nrow = metabolitesLength, ncol = 8)
  colnames(csv.res) <- c("Query", "Match", "HMDB", "PubChem", "ChEBI", "KEGG", "METLIN", "Comment")
  
  for (i in 1:metabolitesLength) {
    if (matchStates[i] == 1) {
      pre.style <- ""
      post.style <- ""
    } else {
      # no matches
      pre.style <- no.prestyle
      post.style <- no.poststyle
    }
    hit <- lilikoi.cmpd[matchIndexes[i], , drop = F]
    html.res[i, ] <- c(paste(pre.style, metaboliteNames[i], post.style, sep = ""), paste(ifelse(matchStates[i] == 
                                                                                                  0, "", matchValues[i]), sep = ""), paste(ifelse(matchStates[i] == 0 || is.na(hit$hmdb_id) || 
                                                                                                                                                    hit$hmdb_id == "" || hit$hmdb_id == "NA", "-", paste("<a href=http://www.hmdb.ca/metabolites/", 
                                                                                                                                                                                                         hit$hmdb_id, " target='_blank'>", hit$hmdb_id, "</a>", sep = "")), sep = ""), paste(ifelse(matchStates[i] == 
                                                                                                                                                                                                                                                                                                      0 || is.na(hit$pubchem_id) || hit$pubchem_id == "" || hit$pubchem_id == "NA", "-", paste("<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=", 
                                                                                                                                                                                                                                                                                                                                                                                               hit$pubchem_id, " target='_blank'>", hit$pubchem_id, "</a>", sep = "")), sep = ""), paste(ifelse(matchStates[i] == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0 || is.na(hit$chebi_id) || hit$chebi_id == "" || hit$chebi_id == "NA", "-", paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     hit$chebi_id, " target='_blank'>", hit$chebi_id, "</a>", sep = "")), sep = ""), paste(ifelse(matchStates[i] == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    0 || is.na(hit$kegg_id) || hit$kegg_id == "" || hit$kegg_id == "NA", "-", paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    hit$kegg_id, " target='_blank'>", hit$kegg_id, "</a>", sep = "")), sep = ""), paste(ifelse(matchStates[i] == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 0 || is.na(hit$metlin_id) || hit$metlin_id == "" || hit$metlin_id == "NA", "-", paste("<a href=http://metlin.scripps.edu/metabo_info.php?molid=", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       hit$metlin_id, " target='_blank'>", hit$metlin_id, "</a>", sep = "")), sep = ""), ifelse(matchStates[i] != 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  1, "View", ""))
    csv.res[i, ] <- c(metaboliteNames[i], ifelse(matchStates[i] == 0, "NA", matchValues[i]), ifelse(matchStates[i] == 
                                                                                                      0, "NA", hit$hmdb_id), ifelse(matchStates[i] == 0, "NA", hit$pubchem_id), ifelse(matchStates[i] == 
                                                                                                                                                                                         0, "NA", hit$chebi_id), ifelse(matchStates[i] == 0, "NA", hit$kegg_id), ifelse(matchStates[i] == 
                                                                                                                                                                                                                                                                          0, "NA", hit$metlin_id), matchStates[i])
  }
  # return only columns user selected
  
  # add query and match columns at the the beginning, and 'Detail' at the end
  returnColumns <- c(TRUE, TRUE, searchOptionsHash$returnColumns, TRUE)
  # html.res <- html.res[, returnColumns, drop = F]
  csv.res <- csv.res[, returnColumns, drop = F]
  csv.res
}

.lilikoi.path_mapping <- function(qvec) {
  #' A .lilikoi.path_mapping Function
  #'
  #' This function allows you to map user metabolites to the its corrsponding pathways
  #' make sure that you have the below database for the exact matching:
  #' Sijia_pathway.rda
  #' @param the list of user metabolites and its corrsponding keeg, pubchem and hmdb ids
  #' @import dplyr
  #' @keywords pathway mapping
  #' @export
  #' @examples .lilikoi.path_mapping(qvec)
  #' .lilikoi.path_mapping(qvec)
  
  
  # if(!exists('sijia_pathway')){
  # load('/home/fadl2/lilikoi/lilikoi_Fadhl/lilikoi/lilikoi/data/Sijia_pathway.rda',.GlobalEnv);
  # load(paste0(getwd(),'/','data/','Sijia_pathway.rda'), .GlobalEnv); } could we add the pathway for
  # these metabolites from Sijia pathway database
  met_W_pathway <- lilikoi.hmdb[match(qvec[, "HMDB"], lilikoi.hmdb$Accession), , drop = F]
  table_result <- data.frame(qvec) %>% mutate(pathway = met_W_pathway$Pathway_Name) %>% arrange(pathway)
  
  .lilikoi.store_result(table_result)
  table_result
}

