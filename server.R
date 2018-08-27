server <- function(input,output,session){
  
  
  
  output$table <- renderDataTable({
    validate(
      need(input$Click1, "")
    )
    if (!is.null(input$Data$datapath)){
      metaboliteMeasurements <<- read.csv(input$Data$datapath, check.names = F, row.names = 1)
      metaboliteNames <<- colnames(metaboliteMeasurements)[-1]
      Samples <<- rownames(metaboliteMeasurements)
      DT::datatable(head(metaboliteMeasurements,20),
                    options = list(scrollX=TRUE))
    }})
  
  output$MTP <- renderTable({
    metabolitePathwayTable <<- lilikoi.metab_to_pathway(metaboliteNames, input$identifier )
    head(metabolitePathwayTable)
  })
  
  output$matchedBox <- renderInfoBox({
    validate(
      need(input$Click1, ""))
    
    match <- length(which(metabolitePathwayTable$Comment==1))
    unmatched <- length(which(!metabolitePathwayTable$Comment==1))
    percentage <- round((match / (match + unmatched)), 1) * 100
    infoBox(
      "Percent Matched", percentage,
      color = if (percentage>=80){"green"}
      else {"yellow"},
      fill = TRUE)
  })
  output$unmatchedBox <- renderInfoBox({
    validate(
      need(input$Click1, ""))
    match <- length(which(metabolitePathwayTable$Comment==1))
    unmatched <- length(which(!metabolitePathwayTable$Comment==1))
    percentage <- round((unmatched / (unmatched + match)), 1) * 100
    infoBox(
      "Percent Unmatched", percentage,
      color = if (percentage<=20){"yellow"}
      else {"yellow"},
      fill=TRUE)
  })
  output$PDS <- renderDataTable({
    validate(
      need(input$Click2, "")
    )
    PDSmatrix <<- lilikoi.get_pd_scores(metaboliteMeasurements, metabolitePathwayTable)
    DT::datatable(t(PDSmatrix),options = list(scrollX=TRUE))
  })
  output$SL <- renderPlot({
    validate(
      need(input$Click3, "")
    )
    selectedPathwaysWeka <<- lilikoi.select_features(PDSmatrix, metaboliteMeasurements, threshold = input$thresh,
                                                     method = "gain")
    ggsave("FeatureSelection.pdf")
    
  })
  output$CL <- renderPlot({
    validate(
      need(input$Click4, "")
    )
    mlplots <<-  lilikoi.machine_learning(PDSmatrix, metaboliteMeasurements, selectedPathwaysWeka)
    result <<- mlplots$res
    ggsave("ClassBest.pdf")
  })
  output$CLPL1 <- renderPlot({
    validate(
      need(input$Click4, "")
    )
    mlplots <<-  lilikoi.machine_learning(PDSmatrix, metaboliteMeasurements, selectedPathwaysWeka)
    return(mlplots$p1)
    ggsave("ClassTrain1.pdf")
  })
  output$CLPL2 <- renderPlot({
    validate(
      need(input$Click4, "")
    )
    return(mlplots$p2)
    ggsave("ClassTrain2.pdf")
  })
  # output$CLPL3 <- renderPlot({
  #   validate(
  #     need(input$Click4, "")
  #   )
  #   return(mlplots$p3)
  #   ggsave("ClassTest.pdf")
  # })
  output$CLPL4 <- renderPlot({
    validate(
      need(input$Click4, "")
    )
    smooth_method <- "binormal"
    plot((mlplots$cart.ROC), col = "red", cex.lab = 1.5)
    par(new = TRUE)
    plot((mlplots$lda.ROC), col = "green", cex.lab = 1.5)
    par(new = TRUE)
    plot((mlplots$svm.ROC), col = "black", cex.lab = 1.5)
    par(new = TRUE)
    plot((mlplots$rf.ROC), col = "orange", cex.lab = 1.5)
    par(new = TRUE)
    plot((mlplots$gbm.ROC), col = "blue", cex.lab = 1.5)
    par(new = TRUE)
    plot((mlplots$pam.ROC), col = "hotpink", cex.lab = 1.5)
    par(new = TRUE)
    plot((mlplots$log.ROC), col = "lightgoldenrod2", main = "Testing ROC", cex.lab = 1.5)
    legend(0.2, 0.4, legend = c("RPART",
                                "LDA", "SVM", "RF", "GBM", "PAM", "LOG"),
           col = c("red", "green","black", "orange", "blue", "hotpink", "lightgoldenrod2"),
           lty = 1:2, cex = 1)
    
    
    ggsave("CLROC.pdf")
  })
  output$CLVI <- renderPlot({
    validate(
      need(input$Click4, "")
    )
    return(plot(plot(varImp(mlplots$res$models[[as.numeric(input$variable)]], scale = FALSE, top = 20))))
    ggsave("VariableImportance.pdf")
  })
  
  output$CF <- renderTable({
    validate(
      need(input$Click5, "")
    )
    if(!is.null(input$Data2$datapath)){
      
      clinicalFactorsData <<- read.csv(file = input$Data2$datapath)
      
      
      CFoptions <<- colnames(clinicalFactorsData)[3:ncol(clinicalFactorsData)]
      
      head(clinicalFactorsData)
    }})
  output$MAPL1 <- renderPlot({
    validate(
      need(input$Click6, "")
    )
    model <<- lilikoi.adjust_model(result, PDSmatrix, selectedPathwaysWeka, metaboliteMeasurements, clinicalFactorsData,
                                   factors = as.character(CFoptions))
    ggsave("Corr.pdf")
  })
  output$MA <- renderPlot({
    validate(
      need(input$Click6, "")
    )
    plot((model$ROCpath), col = "black", cex.lab = 1.5)
    par(new=TRUE)
    plot((model$ROCfac$ROC), col = "red", cex.lab = 1.5)
    par(new=TRUE)
    plot(model$ROCpathfac$ROC, col = "blue", cex.lab = 1.5)
    legend(0.5, 0.4, legend = c("Selected Pathways", "Clinical factors", "Selected pathways + clinical factors"),
           col = c("black", "red", "blue"), lty = 1:2, cex = 1.2)
    ggsave("MAROC.pdf")
  })
  output$MAPL2 <- renderPlot({
    validate(
      need(input$Click6, "")
    )
    model <<- lilikoi.adjust_model(result, PDSmatrix, selectedPathwaysWeka, metaboliteMeasurements, clinicalFactorsData,
                                   factors = as.character(CFoptions))
    return(plot(plot(varImp(model$ROCfac$model, scale = FALSE, top = 20))))
    ggsave("CorrVI.pdf")
  })
  output$MAPL3 <- renderPlot({
    validate(
      need(input$Click6, "")
    )
    return(plot(plot(varImp(model$ROCpathfac$model, scale = FALSE, top = 20))))
    ggsave("CorrVI2.pdf")
  })
  
  
  
  
  output$DownloadCorr <- downloadHandler(
    filename = function() {
      "Corr.pdf"
    },
    content = function(file) {
      file.copy("Corr.pdf", file, overwrite=TRUE)
    }
  )
  output$DownloadMAVI1 <- downloadHandler(
    filename = function() {
      "CorrVI1.pdf"
    },
    content = function(file) {
      file.copy("CorrVI1.pdf", file, overwrite=TRUE)
    }
  )
  output$DownloadMAVI2 <- downloadHandler(
    filename = function() {
      "CorrVI2.pdf"
    },
    content = function(file) {
      file.copy("CorrVI2.pdf", file, overwrite=TRUE)
    }
  )
  output$DownloadCorrROC <- downloadHandler(
    filename = function() {
      "MAROC.pdf"
    },
    content = function(file) {
      file.copy("MAROC.pdf", file, overwrite=TRUE)
    }
  )
  output$DownloadFeatureSelection <- downloadHandler(
    filename = function() {
      "FeatureSelection.pdf"
    },
    content = function(file) {
      file.copy("FeatureSelection.pdf", file, overwrite=TRUE)
    }
  )
  output$DownloadCLBest <- downloadHandler(
    filename = function() {
      "ClassBest.pdf"
    },
    content = function(file) {
      file.copy("ClassBest.pdf", file, overwrite=TRUE)
    }
  )
  output$DownloadCLTrain1 <- downloadHandler(
    filename = function() {
      'ClassTrain1.pdf'
    },
    content = function(file) {
      file.copy("ClassTrain1.pdf", file, overwrite=TRUE)
    }
  )
  output$DownloadCLTrain2 <- downloadHandler(
    filename = function() {
      'ClassTrain2.pdf'
    },
    content = function(file) {
      file.copy("ClassTrain2.pdf", file, overwrite=TRUE)
    }
  )
  output$DownloadCLTest <- downloadHandler(
    filename = function() {
      'ClassTest.pdf'
    },
    content = function(file) {
      file.copy("ClassTest.pdf", file, overwrite=TRUE)
    }
  )
  output$DownloadCLROC <- downloadHandler(
    filename = function() {
      'CLROC.pdf.pdf'
    },
    content = function(file) {
      file.copy("CLROC.pdf.pdf", file, overwrite=TRUE)
    }
  )
  output$DownloadCLVI <- downloadHandler(
    filename = function() {
      'VariableImportance.pdf'
    },
    content = function(file) {
      file.copy("VariableImportance.pdf", file, overwrite=TRUE)
    }
  )
  
  
  url1 <- a("https://github.com/lanagarmire/lilikoi", href="https://github.com/lanagarmire/lilikoi")
  output$tab1 <- renderUI({
    tagList(url1)
  })
  url2 <- a("https://mybinder.org/v2/gh/FADHLyemen/lilikoi_Fadhl/master", href="https://mybinder.org/v2/gh/FADHLyemen/lilikoi_Fadhl/master")
  output$tab2 <- renderUI({
    tagList(url2)
  })
  url3 <- a("https://www.biorxiv.org/content/early/2018/03/16/283408", href="https://www.biorxiv.org/content/early/2018/03/16/283408")
  output$tab3 <- renderUI({
    tagList(url3)
  })

output$Sampledata <- downloadHandler(
    filename = 'Sample.csv',
    content = function(con) {
        data <- read_csv('inst/extdata/plasma_breast_cancer.csv')
      write_csv(data, con)
    }
  )
output$SampledatapubChem <- downloadHandler(
  filename = 'Samplepub.csv',
  content = function(con) {
    data <- read_csv('inst/extdata/plasma_breast_cancer_data_pubChem.csv')
    write_csv(data, con)
  }
)
output$SampledataHMDB <- downloadHandler(
  filename = 'SampleHMDB.csv',
  content = function(con) {
    data <- read_csv('inst/extdata/plasma_breast_cancer_data_HMDB.csv')
    write_csv(data, con)
  }
)
output$SampledataChEBI <- downloadHandler(
  filename = 'SampleChEBI.csv',
  content = function(con) {
    data <- read_csv('inst/extdata/plasma_breast_cancer_data_chEBI.csv')
    write_csv(data, con)
  }
)
output$SampledataKEGG <- downloadHandler(
  filename = 'SampleKegg.csv',
  content = function(con) {
    data <- read_csv('inst/extdata/plasma_breast_cancer_data_KEGG (1).csv')
    write_csv(data, con)
  }
) 

output$SampleMeta<- downloadHandler(
  filename = 'SampleMeta.csv',
  content = function(con) {
    data <- read_csv('inst/extdata/plasma_breast_cancer_Meta.csv')
    write_csv(data, con)
  }
) 

  observeEvent(input$Shinybutton,{
    updateTabItems(session, "tabs", "Logo")})
  observeEvent (input$Click1,{
    updateTabItems(session, "tabs", "Upload")})
  observeEvent (input$Click2,{
    updateTabItems(session, "tabs", "Preview")})
  observeEvent (input$Click3,{
    updateTabItems(session, "tabs", "Mapping")})
  observeEvent (input$Click4,{
    updateTabItems(session, "tabs", "Transformation")})
  observeEvent(input$Click5,{
    updateTabItems(session, "tabs", "Selection")})
  observeEvent(input$Click6,{
    updateTabItems(session, "tabs", "Classification")})
  observeEvent(input$Click7,{
    updateTabItems(session, "tabs", "Clinical")})
  observeEvent(input$Click8,{
    updateTabItems(session, "tabs", "Model")})
  

}