list.of.packages <- c("ggplot2", "caret", "devtools", "dplyr", "hash", "RWeka", "infotheo", "Hmisc", "pROC", "reshape2",
                      "corrplot", "stringr", "Matrix", "randomForest", "glmnet", "gbm", "e1071", "pamr")

lilikoi.develop = function () {
  library(devtools)
  library(pathifier)
  lapply(list.of.packages, library, character.only = TRUE)

  load_all("~/lilikoi")
  setwd("~/lilikoi")
}

lilikoi.develop()

#If you have a problem installing Rweka as it requires Java, you can reconfigure R from the command line by running the below line:
# R CMD javareconf

lilikoi.install = function () {
  # Install packages
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages)) install.packages(new.packages)
  source("https://bioconductor.org/biocLite.R")
  biocLite("pathifier")
}

lilikoi.other = function () {
  # Check package
  devtools::check()

  # Update datasets

  # Saving datasets
  use_data(data.hmdb, data.smpdb, data.metaboAnalyst, overwrite = TRUE)

  install.packages("roxygen2")
  install.packages("testthat")
  install.packages("usethis")
  install.packages("covr")
  install.packages("DT")

  usethis::use_test("get_pd_scores")
  devtools::test()
  lilikoi.codeCoverage()

  devtools::build_win()


  devtools::release()

  # devtools::revdep_check()

  lilikoi.makeDoc(data.hmdb)

  install.packages("lilikoi_0.1.0.tar.gz", repos = NULL, type ="source")

  # installing Rweka on ubuntu may require
  # sudo apt-get install openjdk-6-jdk
  # sudo R CMD javareconf
}

lilikoi.makeDoc = function (dataFrame, title = substitute(dataFrame)) {
  output = c(paste("#'", title), "#' @format data.frame", gsub("^","#'",capture.output(str(dataFrame))), dQuote(title))
  cat(output, sep="\n")
}

lilikoi.codeCoverage = function () {
  library(covr)
  report()
}

# R CMD build .
# R CMD check --as-cran lilikoi_0.1.0.tar.gz
