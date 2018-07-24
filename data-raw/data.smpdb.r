library(tidyverse)
system("wget http://smpdb.ca/downloads/smpdb_metabolites.csv.zip -O temp/smpdb_metabolites.csv.zip")
system("cd temp; unzip smpdb_metabolites.csv.zip")
smpdb = read_csv("temp/metabolites.csv")

# This will be our list of pathways to use
pathways = smpdb$`Pathway Name`

smpdbpathways = data.frame(list(pathway = smpdb$`Pathway Name`, hmdb_id = smpdb$`HMDB ID`))

# smpdb has old shorter version of hmdbids
smpdbpathways$hmdb_id = gsub("HMDB", "HMDB00", smpdbpathways$hmdb_id)

data.smpdb = smpdbpathways

use_data(data.smpdb, overwrite = TRUE)
