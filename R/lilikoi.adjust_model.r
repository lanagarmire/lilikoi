lilikoi.adjust_model <- function(mlResults, PDSmatrix, selected_Pathways_Weka, metaboliteMeasurements, clinical_factors_data,
  factors) {
  #' Model Adjustment method using clinical factors
  #'
  #' This function can adjust the best performing model using additional clinical factors specified by the user
  #' It plots ROC for three models: a model building using only the selected pathways, a model built using the clinical factors
  #' and a model built using both the selected pathways and the clinical factors.
  #'
  #' @param mlResults The models and other results generated from machine_learning method
  #' @param PDSmatrix the PDS matrix generated using lilikoi.get_pd_scores function
  #' @param selected_Pathways_Weka selected pathway using WEKA algorithm generated from lilikoi.select_pathways function
  #' @param metaboliteMeasurements the metaboliteMeasurements for the samples
  #' @param clinical_factors_data Metadata
  #' @param factors which the users want to add to the model
  #' @import corrplot pROC
  #' @importFrom Hmisc rcorr
  #' @export
  prostate_df <- data.frame(t((PDSmatrix[selected_Pathways_Weka, ])), Label = metaboliteMeasurements$Label,
    check.names = T)
  colnames(prostate_df)[which(names(prostate_df) == "Label")] <- "subtype"
  trainDf <- prostate_df[mlResults$train_inx, ]
  testDf <- prostate_df[-mlResults$train_inx, ]

  best_model <- mlResults$models[which.max(mlResults$performance[1, ])]  # the best model has
  # the high AUC
  method <- (unlist(best_model)[[1]])
  # Todo: remove the below statement. We need to refactor the plotting code b/c right now it doesn't
  # work for lda or svmRadial
  method <- "gbm"

  cartClasses <- stats::predict(best_model, newdata = testDf, type = "prob")
  cartClasses1 <- stats::predict(best_model, newdata = testDf)
  cartConfusion <- confusionMatrix(data = unlist(cartClasses1), testDf$subtype)
  # ROC_pathway <-
  # roc(predictor=cartClasses[[1]]$Normal,response=testDf$subtype,levels=rev(levels(testDf$subtype)))
  ROC_pathway <- roc(predictor = as.numeric(unlist(cartClasses[[1]][1])), response = testDf$subtype,
    levels = rev(levels(testDf$subtype)))
  # graphics::plot(smooth(ROC,method='fitdistr'),print.auc=TRUE,col='green')
  smooth_method <- "binormal"

  # pdf('factors.pdf',width=10,height=10)
  graphics::plot(pROC::smooth(ROC_pathway, method = smooth_method), col = "black", cex.lab = 1.5)
  # graphics::plot(ROC_pathway,col='black')
  graphics::par(new = TRUE)
  train_index <- mlResults$train_inx

  factor_data <- cbind(Label = clinical_factors_data[, 2], clinical_factors_data[, factors])
  colnames(factor_data)[which(names(factor_data) == "Label")] <- "subtype"

  ROC_factor <- .lilikoi.create_the_model(factor_data, train_index, method)

  graphics::plot(pROC::smooth(ROC_factor$ROC, method = smooth_method), col = "red", cex.lab = 1.5)
  # graphics::plot(ROC_factor,col='red',print.auc=T)
  graphics::par(new = TRUE)

  pathway_factors_data <- cbind(factor_data[-1], prostate_df)
  colnames(pathway_factors_data)[which(names(pathway_factors_data) == "Label")] <- "subtype"
  # print(head(pathway_factors_data))
  ROC_pathway_factors <- .lilikoi.create_the_model(pathway_factors_data, train_index, method)
  graphics::plot(ROC_pathway_factors$ROC, col = "blue", cex.lab = 1.5)

  # graphics::plot(ROC_pathway_factors,col='blue')

  graphics::legend(0.5, 0.4, legend = c("Selected Pathways", "Clinical factors", "Selected pathways + clinical factors"),
    col = c("black", "red", "blue"), lty = 1:2, cex = 1.2)

  # dev.off()



  # plot the correlation between the pathways and the clinical factors

  df <- select(pathway_factors_data, -"subtype")
  df <- sapply(df, as.numeric)
  all_corr <- Hmisc::rcorr(as.matrix(df), type = "pearson")
  # pdf('corr.pdf',width=20,height=15)
  list(corrPlot = .lilikoi.corr_matrix_plot(all_corr), importancePlotFactors = graphics::plot(graphics::plot(varImp(ROC_factor$model,
    scale = FALSE, top = 20), main = method)), importancePlotPathwayFactors = graphics::plot(graphics::plot(varImp(ROC_pathway_factors$model,
    scale = FALSE, top = 20), main = method)))

}

.lilikoi.corr_matrix_plot <- function(all_corr) {
  corrplot(all_corr$r, method = "circle", type = "upper", order = "hclust", tl.col = "black", tl.srt = 45,
    tl.cex = 1, p.mat = all_corr$P, sig.level = 0.05, insig = "blank", number.cex = 0.3, addrect = 4)
}

.lilikoi.create_the_model <- function(data, train_index, method) {
  #' create the model using caret package
  #'
  #' This function train the model using 80% of the data and retrive the ROC using the 20% of the data
  #' @param data which you want to build the model for
  #' @param train_index theis index was created using lilikoi.machine_learning function and it is required to make sure
  #' that we train and test the data on the same set
  #' @param method which method you want to train the model such as 'lda' and 'gbm'
  res <- list()
  control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
  trainDf <- data[train_index, ]
  testDf <- data[-train_index, ]
  set.seed(7)
  garbage <- suppressWarnings(utils::capture.output(fit <- train(subtype ~ ., data = trainDf, method = method,
    trControl = control, metric = "ROC")))
  cartClasses <- stats::predict(fit, newdata = testDf, type = "prob")
  cartClasses1 <- stats::predict(fit, newdata = testDf)
  cartConfusion <- confusionMatrix(data = cartClasses1, testDf$subtype)
  # graphics::plot(graphics::plot(varImp(fit, scale = FALSE,top=20),main=method))

  ROC <- roc(predictor = as.numeric(unlist(cartClasses[1])), response = testDf$subtype, levels = rev(levels(testDf$subtype)),
    smooth = FALSE)
  res$model <- fit
  res$ROC <- ROC
  res
}
