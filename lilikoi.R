library(shiny)
library(dplyr)
library(ggplot2)
library(reshape2)
library(DT)
library(RWeka)
library(caret)
library(infotheo)
library(shinydashboard)
library(devtools)
library(pROC)
library(gbm)
library(e1071)
library(lilikoi)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(hash)
library(colourpicker)
library(Hmisc)
library(corrplot)
library(rsconnect)


ui <- 
  dashboardPage(skin = "green",
    dashboardHeader(title = ""),
    dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        menuItem("Begin",tabName = "Logo", icon = icon('play')),
        menuItem("Resources", tabName = "Database", icon = icon("database")),
        menuItem("Data Upload", tabName = "Upload", icon = icon("upload")),
        menuItem("Data Preview", tabName = "Preview"),
        menuItem("Feature Mapping", tabName = "Mapping", icon = icon("table")),
        menuItem("Dimension Tranformation", tabName = "Transformation", icon = icon("random")),
        menuItem("Feature Selection", tabName = "Selection", icon = icon("check-circle")),
        menuItem("Classification and Prediction", tabName = "Classification", icon = icon("barcode")),
        menuItem("Clinical Factors", tabName = "Clinical",icon = icon("industry")),
        menuItem("Model Adjustment", tabName = "Model", icon = icon("adjust")),
        menuItem("Download Report", tabName = "Download", icon = icon("arrow-down")))),
  dashboardBody(
    tabItems(
      tabItem(tabName = "Logo",
              fluidPage(
                fluidRow(column(12, align="center",
                                div(style="display: inline-block;",
                                    img(src="LogoScopic.png", height=700, width=700,
                                        tags$h3("A personalized pathway analysis of metabolomic data", style="color: black"),
                                        actionBttn('Click00','Begin', style= "minimal", color = "success"))))))),
      tabItem(tabName = "Database",
              fluidPage(column(12, align="left",
                       div(style="display:inline-block;",
                           tags$h3(tags$b("Resources:")),
                           br(),
                           tags$p(tags$b("Update DataBase"),"link Here"))),
                br(),
                column(12, align="left",
                       div(style="display:inline-block;",
                           tags$p(tags$b("Github:"),"https://github.com/lanagarmire/lilikoi"))),
                br(),
                column(12, align="left",
                       div(style="display:inline-block;",
                           tags$p(tags$b("Mybinder:"),"https://mybinder.org/v2/gh/FADHLyemen/lilikoi_Fadhl/master"))),
                br(),
                column(12, align="left",
                       div(style="display:inline-block;",
                           tags$p(tags$b("Docker"),
                                  "docker pull fadhlyemen/lilikoi", "docker run -d --rm -ti -p 
                            5001:8888 fadhlyemen/lilikoi start-notebook.sh --NotebookApp.token='‘"))),
                br(),
                column(12, align="left",
                       div(style="display:inline-block;",
                           tags$p(tags$b("Citation"),"F. Alakwaa, S. Huang, and L. Garmire, “Lilikoi: 
                            an R package for personalized pathway-based classification modeling using 
                            metabolomics data,” bioRxiv, 2018."))),
                br(),
                column(12, align="left",
                       div(style="display:inline-block;",
                           tags$p(tags$b("Link"),"https://www.biorxiv.org/content/early/2018/03/16/283408"))),
                column(2, offset = 10,
                       actionBttn("reset", "Refresh", style="minimal",color = "success", size="sm"),
                       actionBttn("Click0", label = "Submit",style="minimal",color = "success", size="sm")))),
      tabItem(tabName = "Upload",
              fluidPage(
                  column(12,
                       box(height = '280px', 
                       fileInput("Data","Upload Your Data"),
                       checkboxInput("Ourdata", tags$b("Or Use Our Data:")))),
                  column(2, offset = 10,
                       actionBttn("reset", "Refresh", style="minimal",color = "success", size="sm"),
                       actionBttn("Click1", label = "Submit",style="minimal",color = "success", size="sm")))),
      tabItem(tabName = "Preview",
              fluidPage(
                column(12,
                       box(width='700px',height = '500px', withSpinner(DT::dataTableOutput("table")))),
                column(6,
                       infoBoxOutput("tablestat1")),
                column(6,
                       infoBoxOutput("tablestat2")),
                column(2, offset = 10,
                       actionBttn("reset", "Refresh", style="minimal",color = "success", size="sm"),
                       actionBttn("Click11", label = "Submit",style="minimal",color = "success", size="sm")))),
      
      tabItem(tabName = "Mapping",
              fluidPage(
                
                column(12,
                       tags$h3("Transforms the metabolite names to the matching IDs using Lilikoi MetaTOpathway function")),
              verticalLayout(  
                box(height = '350px',
                       withSpinner(tableOutput("MTP"))),
                box(height = '100px',
                  radioButtons(
                    inputId = "identifier", label = "Please Choose Identifier:", 
                    choices = c("name","pubchem","chebi","kegg","metlin","hmdb"),
                    selected = "name",
                    inline = TRUE)),
                column(12,
                  infoBoxOutput("matchedBox")),
                column(12,
                  infoBoxOutput("unmatchedBox"))),
                column(2,offset = 10,
                       actionBttn("refresh", label = "Refresh", style = "minimal",color = "success", size = "sm"),
                       actionBttn("Click2", label = "Submit",style="minimal", color = "success", size="sm")))),
      
      tabItem(tabName = "Transformation",
              fluidPage(
                tags$h2(("Transforms metabolites into pathway using Pathifier algorithm")),
                tags$h6(tags$i("Pss.. This will take a couple of minutes.")),
                verticalLayout(
                  column(12,
                         box(withSpinner(DT::dataTableOutput("PDS")))),
                  column(2, offset = 10,
                          actionBttn("refresh", label = "Refresh", style = "minimal",color = "success", size = "sm"),
                           actionBttn("Click3", label = "Submit", style="minimal", color = "success", size="sm"))))),
      
      tabItem(tabName = "Selection",
              fluidPage(
                tags$h2(("Selects the most signficant pathway related to phenotype.")),
                tags$h6(tags$i("If you did not get any selected pathwys, you can lower the threshold")),
                verticalLayout(
                  withSpinner(plotOutput("SL")),
                  column(12,
                         sliderInput("thresh", "Threshold",min =0.5 ,max = 0.6,value = 0.54)),
                  column(12,
                         downloadButton("DownloadFeatureSelection")),
                  column(2,offset = 10,
                         fixedRow(
                           actionBttn("refresh", label = "Refresh", style = "minimal",color = "success", size = "sm"),
                           actionBttn("Click4", label = "Submit",style="minimal", color = "success", size="sm")))))),
      
      tabItem(tabName = "Classification",
              fluidPage(
                tags$h5("This function will randomly separate the PDS score matrix with only
                        the selected pathways into training and testing sets. It will use seven
                        widely used machine learning algorithms to build the classification model
                        from the training set. It plots the pathway importance from each model and
                        its accuracy (AUC, sensitivity, specificity)."),
                tabsetPanel(
                  tabPanel("Best Model",
                           box(withSpinner(plotOutput("CL"))),
                           downloadButton("DownloadCLBest")),
                  tabPanel("Training and Testing",
                           box(withSpinner(plotOutput("CLPL1"))),
                           box(withSpinner(plotOutput("CLPL2"))),
                           box(withSpinner(plotOutput("CLPL3"))),
                           tags$ol(
                             tags$p("Download First Plot"),
                             downloadButton("DownloadCLTrain1"),
                             tags$p("Download Second Plot"),
                             downloadButton("DownloadCLTrain2"),
                             tags$p("Download Third Plot"),
                             downloadButton("DownloadCLTest"))),
                           
                  tabPanel("ROC",
                           box(withSpinner(plotOutput("CLPL4"))),
                           downloadButton("DownloadCLROC")),
                  
                  tabPanel("Variable Importance",
                           selectInput(inputId = "variable",
                                       label = "Select Algorithm Importance:",
                                       choices = c("cart"="1",
                                                   "lda"="2",
                                                   "svm"="3",
                                                   "rf"="4",
                                                   "gmb"="5",
                                                   "pam"="6",
                                                   "log"="7"),
                                       selected = NULL), 
                           box(withSpinner(plotOutput("CLVI"))),
                           downloadButton("DownloadCLVI")),
                  column(2, offset = 10,
                         fixedRow(
                           actionBttn("refresh", label = "Refresh", style = "minimal",color = "success", size = "sm"),
                           actionBttn("Click5", label = "Submit",style="minimal",color = "success", size="sm")))))),
      
      tabItem(tabName = "Clinical",
              fluidPage(
                tags$h3("This function include the clinical factors and compute model performance using these factors."),
                tags$h6(tags$i("Model will be created using only the best performed model from the last step.")),
                  
                verticalLayout(
                  column(12,offset = 1,
                         withSpinner(tableOutput("CF"))),
                  column(4,offset = 1,
                         fileInput("Data2","Upload Metadata (Optional)")),
                  column(2, offset = 10,
                         fluidRow(
                           actionBttn("refresh", label = "Refresh", style = "minimal",color = "success", size = "sm"),
                           actionBttn("Click6", label = "Submit",style="minimal", color = "success", size="sm")))))),
      tabItem(tabName = "Model",
              fluidPage(
                tags$h3("Model Adjustemnt using clinical factors such as Age and ethnicity"),
                tabsetPanel(
                  tabPanel("Correlation",
                           withSpinner(plotOutput("MAPL1", width = 800, height = 800))),
                  tabPanel("Variable Importance",
                           box(withSpinner(plotOutput("MAPL2"))),
                           box(withSpinner(plotOutput("MAPL3"))),
                           downloadButton("DownloadCorr")),
                  tabPanel("ROC",
                           withSpinner(plotOutput("MA")))),
                column(2,offset = 10,
                       actionBttn("refresh", label = "Refresh", style = "minimal", size = "sm")))))))
     
      
      
    
        

server <- function(input,output,session){
observeEvent(input$Ourdata,{
    validate(
      need(input$Ourdata,"")
    )
    ReturnData <<- read.csv('inst/extdata/plasma_breast_cancer_data.csv')
    ReturnData
  })
  
output$table <- renderDataTable({
  validate(
    need(input$Click1, "")
    
  )
    if (!is.null(input$Data$datapath)){
      metaboliteMeasurements <<- read.csv(input$Data$datapath, check.names = F, row.names = 1)
      metaboliteNames <<- colnames(metaboliteMeasurements)[-1]
      Samples <<- rownames(metaboliteMeasurements)
      DT::datatable(head(metaboliteMeasurements),
                    options = list(scrollX=TRUE))
}})
  
output$tablestat1 <- renderInfoBox({
      validate(
        need(input$Data$datapath, "")
      )
       infoBox("Metabolites",as.numeric(NROW(metaboliteNames))) 
  })
output$tablestat2 <- renderInfoBox({
      validate(
        need(input$Data$datapath, "")
      )
      infoBox("Samples", as.numeric(NROW(Samples)))
})


output$MTP <- renderTable({
    metabolitePathwayTable <<- lilikoi.metab_to_pathway(metaboliteNames, "name")
    head(metabolitePathwayTable)
})
output$matchedBox <- renderInfoBox({
  validate(
    need(input$Click1, "")
  )
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
    need(input$Click1, "")
  )
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
    PDSmatrix #<<- lilikoi.get_pd_scores(metaboliteMeasurements, metabolitePathwayTable)
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
output$CLPL3 <- renderPlot({
  validate(
    need(input$Click4, "")
  )
    return(mlplots$p3)
    ggsave("ClassTest.pdf")
})
output$CLPL4 <- renderPlot({
  validate(
    need(input$Click4, "")
  )
    plot(smooth(mlplots$cart.ROC, method = "binormal"), col = "red", cex.lab = 1.5)
    par(new=TRUE)
    plot(smooth(mlplots$lda.ROC, method = "binormal"), col = "green", cex.lab = 1.5)
    par(new=TRUE) 
    plot(smooth(mlplots$svm.ROC, method = "binormal"), col = "black", cex.lab = 1.5)
    par(new=TRUE)
    plot(smooth(mlplots$rf.ROC , method = "binormal"), col = "orange", cex.lab = 1.5)
    par(new=TRUE)
    plot(smooth(mlplots$gbm.ROC, method = "binormal"), col = "blue", cex.lab = 1.5)
    par(new=TRUE)
    plot(smooth(mlplots$pam.ROC, method = "binormal"), col = "hotpink", cex.lab = 1.5)
    par(new=TRUE)
    plot(smooth(mlplots$log.ROC, method = "binormal"), col = "lightgoldenrod2", main = "Testing ROC", cex.lab = 1.5)
    legend(0.2, 0.4, legend = c("RPART", "LDA", "SVM", "RF", "GBM", "PAM", "LOG"),
           col = c("red", "green", "black", "orange", "blue", "hotpink", "lightgoldenrod2"), lty = 1:2, cex = 1)
    
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
      head(clinicalFactorsData)
}})
output$MAPL1 <- renderPlot({
  validate(
    need(input$Click6, "")
  )
    model <<- lilikoi.adjust_model(result, PDSmatrix, selectedPathwaysWeka, metaboliteMeasurements, clinicalFactorsData,
                         factors = c("Age", "Race"))
    ggsave("Corr.pdf")
})
output$MA <- renderPlot({
  validate(
    need(input$Click6, "")
  )
    plot(pROC::smooth(model$ROCpath, method = "binormal"), col = "black", cex.lab = 1.5)
    par(new=TRUE)
    plot(pROC::smooth(model$ROCfac$ROC, method = "binormal"), col = "red", cex.lab = 1.5)
    par(new=TRUE)
    plot(model$ROCpathfac$ROC, col = "blue", cex.lab = 1.5)
    legend(0.5, 0.4, legend = c("Selected Pathways", "Clinical factors", "Selected pathways + clinical factors"),
           col = c("black", "red", "blue"), lty = 1:2, cex = 1.2)
})
output$MAPL2 <- renderPlot({
  validate(
    need(input$Click6, "")
  )
  return(plot(plot(varImp(model$ROCfac$model, scale = FALSE, top = 20))))
  ggsave("MA1")
  })
output$MAPL3 <- renderPlot({
  validate(
    need(input$Click6, "")
  )
  return(plot(plot(varImp(model$ROCpathfac$model, scale = FALSE, top = 20))))
  ggsave("MA2")
})




output$DownloadCorr <- downloadHandler(
  filename = function() {
    "Corr.pdf"
  },
  content = function(file) {
    file.copy("Corr.pdf", file, overwrite=TRUE)
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
    'ClassROC.pdf'
  },
  content = function(file) {
    file.copy("ClassROC.pdf", file, overwrite=TRUE)
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



  observeEvent (input$Click00,{
    updateTabItems(session, "tabs", "Database")})
  observeEvent (input$Click0,{
    updateTabItems(session, "tabs", "Upload")})
  observeEvent (input$Click1,{
    updateTabItems(session, "tabs", "Preview")})
  observeEvent (input$Click11,{
     updateTabItems(session, "tabs", "Mapping")})
  observeEvent(input$Click2,{
    updateTabItems(session, "tabs", "Transformation")})
  observeEvent(input$Click3,{
    updateTabItems(session, "tabs", "Selection")})
  observeEvent(input$Click4,{
    updateTabItems(session, "tabs", "Classification")})
  observeEvent(input$Click5,{
    updateTabItems(session, "tabs", "Clinical")})
  observeEvent(input$Click6,{
    updateTabItems(session, "tabs", "Model")})
  
}  

shinyApp(ui=ui,server = server)

