lilikoi.adjust_model <- function(result, PDSmatrix, selected_Pathways_Weka, metaboliteMeasurements, clinical_factors_data, 
  factors) {
  #' Model Adjustemnt function using clinical factors
  #'
  #' This function for adjusted the best performed model using the clinical factors inserted by the user
  #' It plots ROC for three models: model1 build using only seelcted pathways, model2 build using clinical factors
  #' and model3 build using selected pathways and the clinical factors
  #' @param PDSmatrix the PDS matrix geneerated using lilikoi.get_pd_scores function
  #' @param selected_Pathways_Weka selected pathway using WEKA algorithm geneerated from lilikoi.select_features function
  #' @param clinical_factors_data the metaboliteMeasurements for the samples
  #' @param factors which the users want to add to the model
  #' @keywords adjustment
  #' @import Hmisc corrplot
  #' @export
  #' @examples lilikoi.adjust_model(result,PDSmatrix,selected_Pathways_Weka,clinical_factors_data,factors=c('Age','Race'))
  #' lilikoi.adjust_model(result,PDSmatrix,selected_Pathways_Weka,metaboliteMeasurements, clinical_factors_data,factors=c('Age','Race'))
  # 
 
  prostate_df <- data.frame(t((PDSmatrix[selected_Pathways_Weka, ])), Label = metaboliteMeasurements$Label, 
    check.names = T)
  colnames(prostate_df)[which(names(prostate_df) == "Label")] <- "subtype"
  trainDf <- prostate_df[result$train_inx, ]
  testDf <- prostate_df[-result$train_inx, ]
  
  best_model <- result$models[which.max(result$performance[1, ])]  # the best model has
  # the high AUC
  method <- (unlist(best_model)[[1]])
  
  
  cartClasses <- predict(best_model, newdata = testDf, type = "prob")
  cartClasses1 <- predict(best_model, newdata = testDf)
  cartConfusion <- confusionMatrix(data = unlist(cartClasses1), testDf$subtype)
  # ROC_pathway <-
  # roc(predictor=cartClasses[[1]]$Normal,response=testDf$subtype,levels=rev(levels(testDf$subtype)))
  ROC_pathway <- roc(predictor = as.numeric(unlist(cartClasses[[1]][1])), response = testDf$subtype, 
    levels = rev(levels(testDf$subtype)))
  # plot(smooth(ROC,method='fitdistr'),print.auc=TRUE,col='green')
  smooth_method <- "binormal"
  
  # pdf('factors.pdf',width=10,height=10)
  plot((ROC_pathway), col = "black", cex.lab = 1.5)
  # plot(ROC_pathway,col='black')
  par(new = TRUE)
  train_index <- result$train_inx
  
  factor_data <- cbind(Label = clinical_factors_data[, 2], clinical_factors_data[, factors])
  colnames(factor_data)[which(names(factor_data) == "Label")] <- "subtype"
  
  ROC_factor <- .lilikoi.create_the_model(factor_data, train_index, method)
  
  plot((ROC_factor$ROC), col = "red", cex.lab = 1.5)
  # plot(ROC_factor,col='red',print.auc=T)
  par(new = TRUE)
  
  pathway_factors_data <- cbind(factor_data[-1], prostate_df)
  colnames(pathway_factors_data)[which(names(pathway_factors_data) == "Label")] <- "subtype"
  # print(head(pathway_factors_data))
  ROC_pathway_factors <- .lilikoi.create_the_model(pathway_factors_data, train_index, method)
  plot(ROC_pathway_factors$ROC, col = "blue", cex.lab = 1.5)
  
  # plot(ROC_pathway_factors,col='blue')
  
  legend(0.5, 0.4, legend = c("Selected Pathways", "Clinical factors", "Selected pathways + clinical factors"), 
    col = c("black", "red", "blue"), lty = 1:2, cex = 1.2)
  
  # dev.off()
  
  plot(plot(varImp(ROC_factor$model, scale = FALSE, top = 20), main = method))
  plot(plot(varImp(ROC_pathway_factors$model, scale = FALSE, top = 20), main = method))
  # plot the correlation between the pathways and the clinical factors
  
  df <- select(pathway_factors_data, -subtype)
  df <- sapply(df, as.numeric)
  all_corr <- rcorr(as.matrix(df), type = "pearson")
  #pdf('corr.pdf',width=20,height=15)
  corrplot(all_corr$r, method = "circle", type = "upper", order = "hclust", tl.col = "black", tl.srt = 45, 
    tl.cex = 1, p.mat = all_corr$P, sig.level = 0.05, insig = "blank", number.cex = 0.3, addrect = 4)
  #dev.off()
  
  modelplot <- list()
  
  modelplot$ROCpath <- ROC_pathway
  modelplot$ROCfac <- ROC_factor
  modelplot$ROCpathfac <- ROC_pathway_factors
  modelplot$all_corr <- all_corr
  
  modelplot

}

.lilikoi.create_the_model <- function(data, train_index, method) {
  #' create the model using caret package
  #'
  #' This function train the model using 80% of the data and retrive the ROC using the 20% of the data
  #' @param data which you want to build the model for
  #' @param train_index theis index was created using lilikoi.machine_learning function and it is required to make sure
  #' that we train and test the data on the same set
  #' @param method which method you want to train the model such as 'lda' and 'gbm'
  #' @keywords creat the model
  #' @export
  #' @examples .lilikoi.create_the_model(data,train_index,method)
  #' .lilikoi.create_the_model(data,train_index,method)
  res <- list()
  control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
  trainDf <- data[train_index, ]
  testDf <- data[-train_index, ]
  set.seed(7)
  garbage <- suppressWarnings(capture.output(fit <- train(subtype ~ ., data = trainDf, method = method, 
    trControl = control, metric = "ROC")))
  cartClasses <- predict(fit, newdata = testDf, type = "prob")
  cartClasses1 <- predict(fit, newdata = testDf)
  cartConfusion <- confusionMatrix(data = cartClasses1, testDf$subtype)
  # plot(plot(varImp(fit, scale = FALSE,top=20),main=method))
  
  ROC <- roc(predictor = as.numeric(unlist(cartClasses[1])), response = testDf$subtype, levels = rev(levels(testDf$subtype)), 
    smooth = FALSE)
  res$model <- fit
  res$ROC <- ROC
  
  return(res)

  
}
