lilikoi.select_features <- function(PDSmatrix, metaboliteMeasurements, threshold = 0.5, method = "info") {
  #' A lilikoi.select_features Function
  #'
  #' This function allows you to reduce the pathway diemsion using xxxx
  #' @param PDSmatrix from lilikoi.get_pd_scores function
  #' @keywords features selection
  #' @import caret RWeka infotheo ggplot2
  #' @export
  #' @examples selected_Pathways_Weka= lilikoi.select_features(PDSmatrix, metaboliteMeasurements)
  #' lilikoi.select_features(PDSmatrix, metaboliteMeasurements)
  #'
  pds_matrix <- (as.data.frame(cbind(t(PDSmatrix), Label = metaboliteMeasurements$Label)))
  # head(pds_matrix)
  set.seed(2000)
  training_ID <- createDataPartition(pds_matrix$Label, p = 0.8, list = FALSE, times = 1)
  training_diagnosis <- pds_matrix[training_ID, ]
  # head(training_diagnosis)
  if (method == "info") {
    infogainfeatures <- InfoGainAttributeEval(as.logical(training_diagnosis$Label - 1) ~ ., data = training_diagnosis)
    selected_pathways <- names(infogainfeatures[infogainfeatures > threshold])
  } else {
    infogainfeatures <- GainRatioAttributeEval(as.logical(training_diagnosis$Label - 1) ~ ., data = training_diagnosis)
    selected_pathways <- names(infogainfeatures[infogainfeatures > threshold])
  }
  
  
  info.paireddiagnosis.R <- discretize(training_diagnosis[, selected_pathways])
  info.paireddiagnosis.R <- cbind(info.paireddiagnosis.R, as.numeric(as.matrix(training_diagnosis[, 
    ncol(training_diagnosis)])))
  I.R <- mutinformation(info.paireddiagnosis.R, method = "emp")
  I.R.paireddiagnosis <- I.R[, ncol(I.R)]
  # I.R.paireddiagnosis
  theTable <- within(as.data.frame(I.R.paireddiagnosis), I.R.paireddiagnosis <- as.numeric(I.R.paireddiagnosis, 
    levels = names(sort(table(I.R.paireddiagnosis), decreasing = TRUE))))
  # theTable
  theTable <- cbind(row.names(theTable), theTable)
  theTable <- theTable[-ncol(I.R), ]
  colnames(theTable)[1] <- c("name")
  theTable <- transform(theTable, name = reorder(name, order(I.R.paireddiagnosis, decreasing = TRUE)))
  # theTable
  
  p <- ggplot(theTable, aes(name, I.R.paireddiagnosis)) + geom_col() + xlab(NULL) + ylab(NULL)
  
  p + theme(axis.text.x = element_text(angle = 90))
  p + coord_flip()
  q <- p + aes(stringr::str_wrap(name, 20), I.R.paireddiagnosis) + ylab("Mutual information") + xlab("Pathways")
  plot(q + coord_flip())
  
  selected_pathways
}
