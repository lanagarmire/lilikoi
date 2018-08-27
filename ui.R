library(shiny)
library(dplyr)
library(ggplot2)
library(reshape2)
library(DT)
library(RWeka)
library(caret)
library(infotheo)
library(shinydashboard)
library(pROC)
library(gbm)
library(e1071)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(hash)
library(colourpicker)
library(Hmisc)
library(corrplot)
library(rsconnect)
library(devtools)
library(tidyverse)
library(readr)
load_all(getwd())

library(pathifier)

ui <- dashboardPage(skin = "green",
                    dashboardHeader(title = ""),
                    dashboardSidebar(
                      sidebarMenu(
                        id = "tabs",
                        menuItem("Resources", tabName = "Database", icon = icon("database")),
                        menuItem("Begin",tabName = "Logo", icon = icon('play')),
                        menuItem("Data Upload", tabName = "Upload", icon = icon("upload")),
                        menuItem("Data Preview", tabName = "Preview", icon = icon("search")),
                        menuItem("Feature Mapping", tabName = "Mapping", icon = icon("table")),
                        menuItem("Dimension Tranformation", tabName = "Transformation", icon = icon("random")),
                        menuItem("Feature Selection", tabName = "Selection", icon = icon("check-circle")),
                        menuItem("Classification and Prediction", tabName = "Classification", icon = icon("barcode")),
                        menuItem("Clinical Factors", tabName = "Clinical",icon = icon("industry")),
                        menuItem("Model Adjustment", tabName = "Model", icon = icon("adjust")))),
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "Database",
                                fluidPage(column(12, align="left",
                                                 div(style="display:inline-block;",
                                                     tags$h3(tags$b("Introduction:")),
                                                     br(),
                                                     #tags$p(tags$b("Steps to update metabolite database (optional):")),
                                                     tags$ol(
                                                       #tags$li("Open Rstudio"),
                                                       #tags$li("Install and load lilikoi"),
                                                       #tags$li("Run lilikoi.update_database")
                                                       ))),
                                          br(),
                                          column(12, align="left",
                                                 div(style="display:inline-block;",
                                                     tags$p(tags$b("Github:")),
                                                     tags$p(uiOutput("tab1")))),
                                          br(),
                                          column(12, align="left",
                                                 div(style="display:inline-block;",
                                                     tags$p(tags$b("Mybinder:")),
                                                     tags$p(uiOutput("tab2")))),
                                          br(),
                                          column(12, align="left",
                                                 div(style="display:inline-block;",
                                                     tags$p(tags$b("Docker:"),
                                                            tags$p("docker pull fadhlyemen/lilikoi"),
                                                            tags$p("docker run -d --rm -ti -p 5001:8888 fadhlyemen/lilikoi start-notebook.sh --NotebookApp.token='‘")))),
                                          br(),
                                          column(12, align="left",
                                                 div(style="display:inline-block;",
                                                     tags$p(tags$b("Citation:")),
                                                     tags$p(tags$i("F. Alakwaa, S. Huang, and L. Garmire, “Lilikoi: 
                                                            an R package for personalized pathway-based classification modeling using 
                                                            metabolomics data,” bioRxiv, 2018.")))),
                                          br(),
                                          column(12, align="left",
                                                 div(style="display:inline-block;",
                                                     tags$p(tags$b("Link:")),
                                                     tags$p(uiOutput("tab3")))),
                                          column(12, align = "left",
                                                 div(style="display:inline-block;",
                                                     actionBttn("Shinybutton", "Proceed to Shiny APP", style = "material-flat", color = "success"))))),

                        tabItem(tabName = "Logo",
                                fluidPage(
                                  fluidRow(
                                    column(12, 
                                           align="center",
                                              div(style="display: inline-block;",
                                                  img(src="LogoScopic.png", height=700, width=700,
                                                      tags$h3("A personalized pathway analysis of metabolomic data", style="color: black"),
                                                      actionBttn('Click1','Begin', style= "minimal", color = "success"))))))),
                        
                        tabItem(tabName = "Upload",
                                fluidPage(
                                  column(12,
                                         box(height = '300px', 
                                             fileInput("Data","Upload Data"),
                                             tags$p("Or Use Our Data:", downloadLink('Sampledata')),
                                             tags$p("Kegg Sample:", downloadLink("SampledataKEGG")),
                                             tags$p("ChEBI Sample:", downloadLink("SampledataChEBI")),
                                             tags$p("pubChem Sample:", downloadLink("SampledatapubChem")),
                                             tags$p("HMDB Sample:", downloadLink("SampledataHMDB")),
                                             tags$p("META Sample:", downloadLink("SampleMeta")))),
                                  column(2, offset = 10,
                                         actionBttn("Click2", label = "Submit",style="minimal",color = "success", size="sm")))),
                        tabItem(tabName = "Preview",
                            verticalLayout(    
                                fluidPage(column(12,
                                         box(width = "700", withSpinner(DT::dataTableOutput("table")))),
                                  column(6,
                                         infoBoxOutput("tablestat1")),
                                  column(6,
                                         infoBoxOutput("tablestat2"))),
                                fluidPage(column(2, offset = 10,
                                         actionBttn("Click3", label = "Submit",style="minimal",color = "success", size="sm"))))),
                        
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
                                         actionBttn("Click4", label = "Submit",style="minimal", color = "success", size="sm")))),
                        
                        tabItem(tabName = "Transformation",
                                fluidPage(
                                  tags$h2(("Transforms metabolites into pathway using Pathifier algorithm")),
                                  tags$h6(tags$i("Please wait... This will take a couple of minutes.")),
                                  verticalLayout(
                                    column(12,
                                           box(withSpinner(DT::dataTableOutput("PDS")))),
                                    column(12, tags$p(tags$i("Drier Y, Sheffer M and Domany E. Pathway-based personalized analysis of cancer. Proceedings of the National Academy of Sciences. 2013;110 16:6388-93."))),
                                    column(2, offset = 10,
                                           actionBttn("Click5", label = "Submit", style="minimal", color = "success", size="sm"))))),
                        
                        tabItem(tabName = "Selection",
                                fluidPage(
                                  tags$h2(("To Selects the most significant pathways related to phenotype. ")),
                                  tags$h6(tags$i("If you did not get any selected pathwys, you can lower the threshold")),
                                  verticalLayout(
                                    withSpinner(plotOutput("SL")),
                                    column(12,
                                           sliderInput("thresh", "Threshold",min =0.1 ,max = 0.6,value = 0.54)),
                                    column(12,
                                           downloadButton("DownloadFeatureSelection")),
                                    column(2,offset = 10,
                                           fixedRow(
                                             actionBttn("Click6", label = "Submit",style="minimal", color = "success", size="sm")))))),
                        
                        tabItem(tabName = "Classification",
                                fluidPage(
                                  tags$h5("This function will randomly separate the PDS score matrix with only
                                          the selected pathways into training and testing sets. It will use seven
                                          widely used machine learning algorithms to build the classification model
                                          from the training set. It plots the pathway importance from each model and
                                          its accuracy (AUC, sensitivity, specificity)."),
                                  tabsetPanel(
                                    tabPanel("Training",
                                             box(withSpinner(plotOutput("CLPL1")))),
                                    tabPanel("Testing",
                                             box(withSpinner(plotOutput("CLPL2")))),
                                    tabPanel("ROC",
                                             box(withSpinner(plotOutput("CLPL4")))),
                                            
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
                                    tabPanel("Best Model",
                                             box(withSpinner(plotOutput("CL")))),
                                    column(2, offset = 10,
                                           fixedRow(
                                             actionBttn("Click7", label = "Submit",style="minimal",color = "success", size="sm")))))),
                        
                        tabItem(tabName = "Clinical",
                                fluidPage(
                                  tags$h3("This function adjusts the model using the imported clinical factors."),
                                  tags$h6(tags$i("Model will be created using only the best performed model from the last step.")),
                                  
                                  verticalLayout(
                                    column(12,offset = 1,
                                           withSpinner(tableOutput("CF"))),
                                    column(4,offset = 1,
                                           fileInput("Data2","Upload Metadata (Optional)")),
                                    column(2, offset = 10,
                                           fluidRow(
                                             actionBttn("Click8", label = "Submit",style="minimal", color = "success", size="sm")))))),
                        tabItem(tabName = "Model",
                                fluidPage(
                                  tags$h3("Model Adjustemnt using clinical factors such as Age and ethnicity"),
                                  tabsetPanel(
                                    tabPanel("Variable Importance",
                                             box(withSpinner(plotOutput("MAPL2"))),
                                             box(withSpinner(plotOutput("MAPL3")))),
                                    tabPanel("ROC",
                                             withSpinner(plotOutput("MA"))),
                                    tabPanel("Correlation",
                                             withSpinner(plotOutput("MAPL1", width = 800, height = 800)))))))))
                        
