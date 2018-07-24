lilikoi.update_database <- function(dataRawPath = paste0(getwd(), "/data-raw/")) {
  #' A lilikoi.update_database Function
  #'
  #' This function allows you to update the database sources used
  #' @param dataRawPath String path to the location of the data-raw folder.
  #' @export

  wd = getwd()
  setwd(dataRawPath)
  source("data.metaboAnalyst.r")
  source("data.smpdb.r")
  source("data.hmdb.r")
  setwd(wd)
}
