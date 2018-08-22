library(tidyverse)
print("Notice: you may need to manually update the MetaboAnalyst link from: http://www.metaboanalyst.ca/faces/docs/Resources.xhtml")
system("wget https://www.dropbox.com/s/adkps7jq810dyl4/MetaboAnalyst-4.09.war -O temp/meta.war")
system("cd temp; jar xvf meta.war")

# Strip non-ascii characters because R sucks at utf8
system("perl -i.bak -pe 's/[^[:ascii:]]//g' temp/resources/libs/cmpd_name.csv")

data.metaboAnalyst = read.csv("temp/resources/libs/cmpd_name.csv") %>% separate_rows("synonym", sep = "; *")

data.metaboAnalyst$lipid = NULL
use_data(data.metaboAnalyst, overwrite = TRUE)

# Add SMPDB ids to our SMPDB table. Todo: do this directly from SMPDB (which is currently down for indexing).

system("perl -i.bak -pe 's/[^[:ascii:]]//g' temp/resources/libs/smp_path.csv")
data.metaboAnalystSmpPath = read.csv("temp/resources/libs/smp_path.csv")
data.metaboAnalystSmpPath$V2 = NULL
data.metaboAnalystSmpPath$id = substr(data.metaboAnalystSmpPath$V3, 26, 40)
data.metaboAnalystSmpPath$V3 = NULL


findCell = function (dataframe, searchString, idColumnName, returnColumnName) {
  matches = dataframe[as.character(dataframe[[idColumnName]]) == searchString,][[returnColumnName]]
  if (length(matches) == 0)
    return("-1")
  matches[[1]]
}

data.smpdb$ids = lapply(data.smpdb$pathway, function (pathway) findCell(data.metaboAnalystSmpPath, pathway, "V1", "id"))
data.smpdb$ids = unlist(data.smpdb$ids)
use_data(data.smpdb, overwrite = TRUE)

