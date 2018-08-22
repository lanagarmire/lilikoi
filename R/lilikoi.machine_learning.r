lilikoi.machine_learning <- function(PDSmatrix, measurementLabels, significantPathways, dependentVariableName = "Normal") {
  #' A machine learning Function
  #'
  #' This function for classification using 7 different machine learning algorithms
  #' and it plots the ROC curves and the AUC, SEN, and specificity
  #' @param PDSmatrix from lilikoi.get_pd_scores
  #' @param measurementLabels Cancer or not
  #' @param significantPathways Pathways with high significance from lilikoi.select_pathways function
  #' @import caret pROC ggplot2 reshape2 gbm
  #' @export
  cancer_df <- data.frame(t((PDSmatrix[significantPathways, ])), Label = measurementLabels, check.names = T)
  colnames(cancer_df)[which(names(cancer_df) == "Label")] <- "subtype"

  performance_training <- matrix(rep(0, len = 21), nrow = 3)  # AUC SENS SPECF
  performance_testing <- matrix(rep(0, len = 56), nrow = 8)  # ROC SENS SPEC

  performance <- matrix(rep(0, len = 56), nrow = 8)  # ROC SENS SPEC

  models <- list()

  ############### Shuffle stat first rand <- sample(nrow(cancer_df)) cancer_df=cancer_df[rand, ]

  ############### Randomly Split the data in to training and testing
  set.seed(2000)
  trainIndex <- createDataPartition(cancer_df$subtype, p = 0.8, list = FALSE, times = 1)
  trainDf <- cancer_df[trainIndex, ]
  testDf <- cancer_df[-trainIndex, ]
  # trainDf$subtype=as.factor(paste0('X',trainDf$subtype))
  # testDf$subtype=as.factor(paste0('X',testDf$subtype)) Training and tunning parameters prepare
  # training scheme
  control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)

  # 1- RPART ALGORITHM
  set.seed(7)  # This ensures that the same resampling sets are used,which will
  # come in handy when we compare the resampling profiles between models.

  # supress the warning messgae options(warn=-1) options(warn=0) ?suppressWarnings()

  garbage <- utils::capture.output(fit.cart <- train(subtype ~ ., data = trainDf, method = "rpart",
    trControl = control, metric = "ROC"))
  # fit.cart <- train(subtype~., data=trainDf, method = 'rpart', trControl=control,metric='ROC')
  # #loclda
  models[[1]] <- fit.cart
  performance_training[1, 1] <- max(fit.cart$results$ROC)  # AUC
  performance_training[2, 1] <- fit.cart$results$Sens[which.max(fit.cart$results$ROC)]  # sen
  performance_training[3, 1] <- fit.cart$results$Spec[which.max(fit.cart$results$ROC)]  # spec

  # Model Testing
  cartClasses <- stats::predict(fit.cart, newdata = testDf, type = "prob")
  cartClasses1 <- stats::predict(fit.cart, newdata = testDf)
  cartConfusion <- confusionMatrix(data = cartClasses1, testDf$subtype)

  cart.ROC <- roc(predictor = cartClasses[[dependentVariableName]], response = testDf$subtype, levels = rev(levels(testDf$subtype)))

  performance_testing[1, 1] <- as.numeric(cart.ROC$auc)  # AUC
  performance_testing[2, 1] <- cartConfusion$byClass[1]  # SENS
  performance_testing[3, 1] <- cartConfusion$byClass[2]  # SPEC
  performance_testing[4, 1] <- cartConfusion$overall[1]  # accuracy
  performance_testing[5, 1] <- cartConfusion$byClass[5]  # precision
  performance_testing[6, 1] <- cartConfusion$byClass[6]  # recall = sens
  performance_testing[7, 1] <- cartConfusion$byClass[7]  # F1
  performance_testing[8, 1] <- cartConfusion$byClass[11]  # BALANCED ACCURACY

  # 2-LDA ALGORITHM
  set.seed(7)
  # assign(paste0('fit.lda',k),train(subtype~., data=trainDf, method='pls',
  # trControl=control,metric='ROC'))
  garbage <- suppressWarnings(utils::capture.output(fit.lda <- train(subtype ~ ., data = trainDf, method = "lda",
    trControl = control, metric = "ROC", trace = F)))  # loclda)
  # fit.lda <- train(subtype~., data=trainDf, method = 'lda', trControl=control,metric='ROC') #loclda
  models[[2]] <- fit.lda
  performance_training[1, 2] <- max(fit.lda$results$ROC)  # AUC
  performance_training[2, 2] <- fit.lda$results$Sens[which.max(fit.lda$results$ROC)]  # sen
  performance_training[3, 2] <- fit.lda$results$Spec[which.max(fit.lda$results$ROC)]  # spec

  # model Testing
  ldaClasses <- stats::predict(fit.lda, newdata = testDf, type = "prob")
  ldaClasses1 <- stats::predict(fit.lda, newdata = testDf)
  ldaConfusion <- confusionMatrix(data = ldaClasses1, testDf$subtype)

  lda.ROC <- roc(predictor = ldaClasses[[dependentVariableName]], response = testDf$subtype, levels = rev(levels(testDf$subtype)))

  performance_testing[1, 2] <- as.numeric(lda.ROC$auc)  # AUC
  performance_testing[2, 2] <- ldaConfusion$byClass[1]  # SENS
  performance_testing[3, 2] <- ldaConfusion$byClass[2]  # SPEC
  performance_testing[4, 2] <- ldaConfusion$overall[1]  # accuracy
  performance_testing[5, 2] <- ldaConfusion$byClass[5]  # precision
  performance_testing[6, 2] <- ldaConfusion$byClass[6]  # recall = sens
  performance_testing[7, 2] <- ldaConfusion$byClass[7]  # F1
  performance_testing[8, 2] <- ldaConfusion$byClass[11]  # BALANCED ACCURACY

  # 3- SVM ALGORITHM
  set.seed(7)
  garbage <- utils::capture.output(fit.svm <- train(subtype ~ ., data = trainDf, method = "svmRadial",
    trControl = control, metric = "ROC"))
  # fit.svm <- train(subtype~., data=trainDf, method='svmRadial', trControl=control,metric='ROC')
  # assign(paste0('fit.svm',k),train(subtype~., data=trainDf, method='svmRadical',
  # trControl=control,metric='ROC'))
  models[[3]] <- fit.svm
  performance_training[1, 3] <- max(fit.svm$results$ROC)  # AUC
  performance_training[2, 3] <- fit.svm$results$Sens[which.max(fit.svm$results$ROC)]  # sen
  performance_training[3, 3] <- fit.svm$results$Spec[which.max(fit.svm$results$ROC)]  # spec

  # Model Testing
  svmClasses <- stats::predict(fit.svm, newdata = testDf, type = "prob")
  svmClasses1 <- stats::predict(fit.svm, newdata = testDf)
  svmConfusion <- confusionMatrix(data = svmClasses1, testDf$subtype)

  svm.ROC <- roc(predictor = svmClasses[[dependentVariableName]], response = testDf$subtype, levels = rev(levels(testDf$subtype)))

  performance_testing[1, 3] <- as.numeric(svm.ROC$auc)  # AUC
  performance_testing[2, 3] <- svmConfusion$byClass[1]  # SENS
  performance_testing[3, 3] <- svmConfusion$byClass[2]  # SPEC
  performance_testing[4, 3] <- svmConfusion$overall[1]  # accuracy
  performance_testing[5, 3] <- svmConfusion$byClass[5]  # precision
  performance_testing[6, 3] <- svmConfusion$byClass[6]  # recall = sens
  performance_testing[7, 3] <- svmConfusion$byClass[7]  # F1
  performance_testing[8, 3] <- svmConfusion$byClass[11]  # BALANCED ACCURACY

  # 4-RF ALGORITHM
  set.seed(7)
  garbage <- utils::capture.output(fit.rf <- train(subtype ~ ., data = trainDf, method = "rf", trControl = control,
    metric = "ROC"))
  # fit.rf <- train(subtype~., data=trainDf, method='rf', trControl=control,metric='ROC')

  models[[4]] <- fit.rf
  performance_training[1, 4] <- max(fit.rf$results$ROC)  # AUC
  performance_training[2, 4] <- fit.rf$results$Sens[which.max(fit.rf$results$ROC)]  # sen
  performance_training[3, 4] <- fit.rf$results$Spec[which.max(fit.rf$results$ROC)]  # spec

  # Model Testing
  rfClasses <- stats::predict(fit.rf, newdata = testDf, type = "prob")
  rfClasses1 <- stats::predict(fit.rf, newdata = testDf)
  rfConfusion <- confusionMatrix(data = rfClasses1, testDf$subtype)


  rf.ROC <- roc(predictor = rfClasses[[dependentVariableName]], response = testDf$subtype, levels = rev(levels(testDf$subtype)))

  performance_testing[1, 4] <- as.numeric(rf.ROC$auc)  # AUC
  performance_testing[2, 4] <- rfConfusion$byClass[1]  # SENS
  performance_testing[3, 4] <- rfConfusion$byClass[2]  # SPEC
  performance_testing[4, 4] <- rfConfusion$overall[1]  # accuracy
  performance_testing[5, 4] <- rfConfusion$byClass[5]  # precision
  performance_testing[6, 4] <- rfConfusion$byClass[6]  # recall = sens
  performance_testing[7, 4] <- rfConfusion$byClass[7]  # F1
  performance_testing[8, 4] <- rfConfusion$byClass[11]  # BALANCED ACCURACY

  # 5- GBM ALGORITHM
  set.seed(7)
  garbage <- suppressWarnings(utils::capture.output(fit.gbm <- train(subtype ~ ., data = trainDf, method = "gbm",
    trControl = control, metric = "ROC")))
  # fit.gbm <- train(subtype~., data=trainDf, method='gbm', trControl=control,metric='ROC')
  models[[5]] <- fit.gbm
  performance_training[1, 5] <- max(fit.gbm$results$ROC)  # AUC
  performance_training[2, 5] <- fit.gbm$results$Sens[which.max(fit.gbm$results$ROC)]  # sen
  performance_training[3, 5] <- fit.gbm$results$Spec[which.max(fit.gbm$results$ROC)]  # spec

  # Model Testing
  gbmClasses <- stats::predict(fit.gbm, newdata = testDf, type = "prob")
  gbmClasses1 <- stats::predict(fit.gbm, newdata = testDf)
  gbmConfusion <- confusionMatrix(data = gbmClasses1, testDf$subtype)


  gbm.ROC <- roc(predictor = gbmClasses[[dependentVariableName]], response = testDf$subtype, levels = rev(levels(testDf$subtype)))

  performance_testing[1, 5] <- as.numeric(gbm.ROC$auc)  # AUC
  performance_testing[2, 5] <- gbmConfusion$byClass[1]  # SENS
  performance_testing[3, 5] <- gbmConfusion$byClass[2]  # SPEC
  performance_testing[4, 5] <- gbmConfusion$overall[1]  # accuracy
  performance_testing[5, 5] <- gbmConfusion$byClass[5]  # precision
  performance_testing[6, 5] <- gbmConfusion$byClass[6]  # recall = sens
  performance_testing[7, 5] <- gbmConfusion$byClass[7]  # F1
  performance_testing[8, 5] <- gbmConfusion$byClass[11]  # BALANCED ACCURACY

  # 6- PAM ALGORITHM
  set.seed(7)
  garbage <- utils::capture.output(fit.pam <- train(subtype ~ ., data = trainDf, method = "pam", trControl = control,
    metric = "ROC"))  # plr) #loclda)
  # fit.pam <- train(subtype~., data=trainDf, method='pam', trControl=control,metric='ROC')#plr
  models[[6]] <- fit.pam
  performance_training[1, 6] <- max(fit.pam$results$ROC)  # AUC
  performance_training[2, 6] <- fit.pam$results$Sens[which.max(fit.pam$results$ROC)]  # sen
  performance_training[3, 6] <- fit.pam$results$Spec[which.max(fit.pam$results$ROC)]  # spec

  # model Testing
  pamClasses <- stats::predict(fit.pam, newdata = testDf, type = "prob")
  pamClasses1 <- stats::predict(fit.pam, newdata = testDf)
  pamConfusion <- confusionMatrix(data = pamClasses1, testDf$subtype)


  pam.ROC <- roc(predictor = pamClasses[[dependentVariableName]], response = testDf$subtype, levels = rev(levels(testDf$subtype)))

  performance_testing[1, 6] <- as.numeric(pam.ROC$auc)  # AUC
  performance_testing[2, 6] <- pamConfusion$byClass[1]  # SENS
  performance_testing[3, 6] <- pamConfusion$byClass[2]  # SPEC
  performance_testing[4, 6] <- pamConfusion$overall[1]  # accuracy
  performance_testing[5, 6] <- pamConfusion$byClass[5]  # precision
  performance_testing[6, 6] <- pamConfusion$byClass[6]  # recall = sens
  performance_testing[7, 6] <- pamConfusion$byClass[7]  # F1
  performance_testing[8, 6] <- pamConfusion$byClass[11]  # BALANCED ACCURACY

  # 7- logistic regression

  set.seed(7)
  garbage <- suppressWarnings(utils::capture.output(fit.log <- train(subtype ~ ., data = trainDf, method = "glmnet",
    trControl = control, metric = "ROC")))
  # fit.log <- train(subtype~., data=trainDf, method='glm', trControl=control,metric='ROC')#
  models[[7]] <- fit.log
  performance_training[1, 7] <- max(fit.log$results$ROC)  # AUC
  performance_training[2, 7] <- fit.log$results$Sens[which.max(fit.log$results$ROC)]  # sen
  performance_training[3, 7] <- fit.log$results$Spec[which.max(fit.log$results$ROC)]  # spec

  # Model Testing
  logClasses <- stats::predict(fit.log, newdata = testDf, type = "prob")
  logClasses1 <- stats::predict(fit.log, newdata = testDf)
  logConfusion <- confusionMatrix(data = logClasses1, testDf$subtype)
  log.ROC <- roc(predictor = logClasses[[dependentVariableName]], response = testDf$subtype, levels = rev(levels(testDf$subtype)))

  performance_testing[1, 7] <- as.numeric(log.ROC$auc)  # AUC
  performance_testing[2, 7] <- logConfusion$byClass[1]  # SENS
  performance_testing[3, 7] <- logConfusion$byClass[2]  # SPEC
  performance_testing[4, 7] <- logConfusion$overall[1]  # accuracy
  performance_testing[5, 7] <- logConfusion$byClass[5]  # precision
  performance_testing[6, 7] <- logConfusion$byClass[6]  # recall = sens
  performance_testing[7, 7] <- logConfusion$byClass[7]  # F1
  performance_testing[8, 7] <- logConfusion$byClass[11]  # BALANCED ACCURACY

  # performance_training=matrix( rep( 0, len=21), nrow = 3) #AUC SENS SPECF performance_testing=matrix(
  # rep( 0, len=56), nrow = 8) # ROC SENS SPEC

  ##### plot the variable importance graphics::par(mfrow=c(7,1))
  p5 <- graphics::plot(graphics::plot(varImp(fit.cart, scale = FALSE, top = 20), main = "RPART"))
  p6 <- graphics::plot(graphics::plot(varImp(fit.lda, scale = FALSE, top = 20), main = "LDA"))
  p7 <- graphics::plot(graphics::plot(varImp(fit.svm, scale = FALSE, top = 20), main = "SVM"))
  p8 <- graphics::plot(graphics::plot(varImp(fit.rf, scale = FALSE, top = 20), main = "RF"))
  p9 <- graphics::plot(graphics::plot(varImp(fit.gbm, scale = FALSE, top = 20), main = "GBM"))
  p10 <- graphics::plot(graphics::plot(varImp(fit.pam, scale = FALSE, top = 20), main = "PAM"))
  p11 <- graphics::plot(graphics::plot(varImp(fit.log, scale = FALSE, top = 20), main = "LOG"))

  ############# plot ROC
  smooth_method <- "binormal"  # 'density'

  graphics::plot(cart.ROC, col = "red", cex.lab = 1.5)
  graphics::par(new = TRUE)
  graphics::plot(lda.ROC, col = "green", cex.lab = 1.5)
  graphics::par(new = TRUE)
  graphics::plot(svm.ROC, col = "black", cex.lab = 1.5)
  graphics::par(new = TRUE)
  graphics::plot(rf.ROC, col = "orange", cex.lab = 1.5)
  graphics::par(new = TRUE)
  graphics::plot(gbm.ROC, col = "blue", cex.lab = 1.5)
  graphics::par(new = TRUE)
  graphics::plot(pam.ROC, col = "hotpink", cex.lab = 1.5)
  graphics::par(new = TRUE)

  graphics::plot(log.ROC, col = "lightgoldenrod2", main = "Testing ROC", cex.lab = 1.5)

  graphics::legend(0.2, 0.4, legend = c("RPART", "LDA", "SVM", "RF", "GBM", "PAM", "LOG"), col = c("red",
    "green", "black", "orange", "blue", "hotpink", "lightgoldenrod2"), lty = 1:2, cex = 1)
  # dev.off()

  ###################### performance plotting
  performance_training_list <- list()
  performance_testing_list <- list()
  performance_testing_list[[1]] <- performance_testing
  performance_training_list[[1]] <- performance_training
  list_test <- performance_testing_list
  list_train <- performance_training_list

  AUC_train <- lapply(list_train, function(x) x[1, ])
  AUC_test <- lapply(list_test, function(x) x[1, ])

  SENS_train <- lapply(list_train, function(x) x[2, ])
  SENS_test <- lapply(list_test, function(x) x[2, ])

  SPEC_train <- lapply(list_train, function(x) x[3, ])
  SPEC_test <- lapply(list_test, function(x) x[3, ])

  F1_test <- lapply(list_test, function(x) x[7, ])
  Balanced_accuracy_test <- lapply(list_test, function(x) x[8, ])


  output1 <- do.call(rbind, lapply(AUC_train, matrix, ncol = 7, byrow = TRUE))
  output2 <- do.call(rbind, lapply(AUC_test, matrix, ncol = 7, byrow = TRUE))

  output3 <- do.call(rbind, lapply(SENS_train, matrix, ncol = 7, byrow = TRUE))
  output4 <- do.call(rbind, lapply(SENS_test, matrix, ncol = 7, byrow = TRUE))

  output5 <- do.call(rbind, lapply(SPEC_train, matrix, ncol = 7, byrow = TRUE))
  output6 <- do.call(rbind, lapply(SPEC_test, matrix, ncol = 7, byrow = TRUE))

  output7 <- do.call(rbind, lapply(F1_test, matrix, ncol = 7, byrow = TRUE))
  output8 <- do.call(rbind, lapply(Balanced_accuracy_test, matrix, ncol = 7, byrow = TRUE))

  AUC_train_mean <- apply(output1, 2, mean)
  AUC_test_mean <- apply(output2, 2, mean)
  AUC <- data.frame(AUC = t(cbind(t(AUC_train_mean), t(AUC_test_mean))))


  SENS_train_mean <- apply(output3, 2, mean)
  SENS_test_mean <- apply(output4, 2, mean)
  SENS <- data.frame(SENS = t(cbind(t(SENS_train_mean), t(SENS_test_mean))))

  SPEC_train_mean <- apply(output5, 2, mean)
  SPEC_test_mean <- apply(output6, 2, mean)
  SPEC <- data.frame(SPEC = t(cbind(t(SPEC_train_mean), t(SPEC_test_mean))))

  F1_test_mean <- apply(output7, 2, mean)
  F1 <- data.frame(F1 = t(t(F1_test_mean)))

  Balanced_accuracy_test_mean <- apply(output8, 2, mean)
  Balanced_accuracy <- data.frame(Balanced_accuracy = t(t(Balanced_accuracy_test_mean)))


  trainingORtesting <- t(cbind(t(rep("training", 7)), t(rep("testing", 7))))
  testing_only <- t(t(rep("testing", 7)))

  performance_data <- data.frame(AUC = AUC, SENS = SENS, SPEC = SPEC, trainingORtesting, Algorithm = (rep(t(c("RPART",
    "LDA", "SVM", "RF", "GBM", "PAM", "LOG")), 2)))

  performance_data_test <- data.frame(AUC = data.frame(AUC = t((t(AUC_test_mean)))), SENS = data.frame(SENS = t((t(SENS_test_mean)))),
    SPEC = data.frame(SPEC = t((t(SPEC_test_mean)))), F1 = F1, Balanced_accuracy = Balanced_accuracy,
    testing_only, Algorithm = (rep(t(c("RPART", "LDA", "SVM", "RF", "GBM", "PAM", "LOG")), 1)))

  # print(performance_data_test)

  # performance_data
  melted_performance_data <- suppressMessages(melt(performance_data))
  melted_performance_data_test <- suppressMessages(melt(performance_data_test))
  # melted_performance_data


  # We do this purely to fix the R CMD check message, due to us using name in NSE below.
  Algorithm <- NULL
  value <- NULL
  variable <- NULL

  textLabels <- geom_text(aes(x = Algorithm, label = round(value, 2), fill = variable), position = position_dodge(width = 1),
    vjust = -0.5, size = 2)


  # pdf('pdf1.pdf',width=10,height=10)
  p1 <- ggplot(data = melted_performance_data[trainingORtesting == "training", ], aes(x = Algorithm,
    y = value, fill = variable)) + geom_bar(stat = "identity", position = position_dodge()) + xlab("") +
    ylab("") + ggtitle("Training") + theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 15,
    face = "bold"), axis.title = element_text(size = 14, face = "bold")) + labs(fill = "") + textLabels
  print(p1)

  # dev.off()

  # pdf('pdf2.pdf',width=10,height=10)
  p2 <- ggplot(data = melted_performance_data[trainingORtesting == "testing", ], aes(x = Algorithm,
    y = value, fill = variable)) + geom_bar(stat = "identity", position = position_dodge()) + xlab("") +
    ylab("") + ggtitle("Testing") + theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 15,
    face = "bold"), axis.title = element_text(size = 14, face = "bold")) + labs(fill = "") + textLabels
  print(p2)
  # dev.off()

  # pdf('pdf3.pdf',width=10,height=10)
  p3 <- ggplot(data = melted_performance_data_test, aes(x = Algorithm, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) + xlab("") + ylab("") + ggtitle("Testing") +
    theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 14, face = "bold")) + labs(fill = "") + textLabels
  print(p3)
  # dev.off()

  # Which algorithm performs better based on the its AUC on testing

  # the best algorithms
  best_model <- models[which.max(performance_testing[1, ])]  # the best model has
  # the high AUC
  method <- (unlist(best_model)[[1]])

  if (method == "glmnet") {
    method <- "log"
  }
  if (method == "svmRadial") {
    method <- "svm"
  }

  dd <- filter(melted_performance_data_test, Algorithm == toupper(method))
  p4 <- ggplot(data = dd, aes(x = Algorithm, y = value, fill = variable)) + geom_bar(stat = "identity",
    position = position_dodge()) + xlab("") + ylab("") + ggtitle("Testing") + theme(plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 14, face = "bold")) +
    labs(fill = "")
  print(p4)

  mlResults <- list()
  mlResults$models <- models
  mlResults$performance <- performance_testing
  mlResults$train_inx <- trainIndex

  plots <- list()
  plots$mlResults <- mlResults
  plots$p1 <- p1
  plots$p2 <- p2
  plots$p3 <- p3
  plots$p4 <- p4
  plots$p5 <- p5
  plots$p6 <- p6
  plots$p7 <- p7
  plots$p8 <- p8
  plots$p9 <- p9
  plots$p10 <- p10
  plots$p11 <- p11

  plots
}
