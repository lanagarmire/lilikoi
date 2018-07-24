library(tidyverse)
system("wget https://www.dropbox.com/s/adkps7jq810dyl4/MetaboAnalyst-4.09.war -O temp/meta.war")
system("cd temp; jar xvf meta.war")

# Strip non-ascii characters because R sucks at utf8
system("sed -i 's/[\d128-\d255]//g' temp/resources/libs/cmpd_name.csv")

data.metaboAnalyst = read.csv("temp/resources/libs/cmpd_name.csv") %>% separate_rows("synonym", sep = "; *")
data.metaboAnalyst$lipid = NULL
use_data(data.metaboAnalyst, overwrite = TRUE)
