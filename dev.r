list.of.packages <- c("ggplot2", "caret", "devtools", "dplyr", "hash", "RWeka", "infotheo", "pROC", "reshape2", 
  "corrplot", "stringr", "Hmisc", "Matrix", "randomForest", "glmnet", "gbm", "e1071", "pamr")

# Load packages
lapply(list.of.packages, library, character.only = TRUE)
library(pathifier)

# Load package where ~/ is path to your local lilikoi repo
library(devtools)
load_all("~/lilikoi")
setwd("~/lilikoi")

# 
install_github("lanagarmire/lilikoi")

# Saving datasets
use_data(lilikoi.cmpd, lilikoi.hmdb, lilikoi.syn, lilikoi.pathways, lilikoi.metabolites, overwrite = TRUE)

#Fadhl_modifications
#lilikoi.metabolites=lapply(lilikoi.metabolites,function (x)  {gsub("HMDB", "HMDB00", x)})
#use_data(lilikoi.metabolites, overwrite = TRUE)

# Install packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)
source("https://bioconductor.org/biocLite.R")
biocLite("pathifier")

# using older version of prinCurve package
devtools::install_version("princurve", version = "1.1-12", repos = "http://cran.us.r-project.org")

# Adding Bioconductor to R repos
options(repos = BiocInstaller::biocinstallRepos())
getOption("repos")
