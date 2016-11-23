library(shiny)
library(tiff)
#
library("parallel")
library("foreach")
library("doParallel")
library("iterators")
#

ui = fluidPage(
  
  includeCSS("www/bootstrap.css"),
  includeCSS("www/Progress.css"),
  titlePanel("ColocalizR"),
  
  tabsetPanel(
    tabPanel("Plate information",
             fluidRow(
               column(2,align='left',
                      h3('Plate information'),
                      wellPanel(
                        radioButtons("MulPlates", label = "Multiple plates ?",inline = T,
                                     choices = c("YES","NO"),selected = "NO"),
                        strong('Image location :'),
                        h4(''),
                        actionButton("folder", label = 'Select folder'),
                        h4(''),
                        uiOutput("PlateIDs")
                      ),
                      
                      wellPanel(
                        radioButtons("PLATEMAP", label = "Plate map :",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        conditionalPanel(
                          condition = "input.PLATEMAP == 'YES'",  
                          fileInput(inputId = 'platemaps', 
                                    label = 'Select plate map (CSV)',
                                    multiple = TRUE,
                                    accept=c(".csv"))
                        )
                      ),
                      br(),
                      (actionButton("ImpImg", label='IMPORT PLATE INFO',style='padding:10px; font-size:150%;align:center',width=290))
               ),
               
               column(2,align='left', offset = 0,
                      
                      h3('Colocalization Settings'),
                      wellPanel(
                        radioButtons("TimeCourse", label = "Time course experiment ?",inline = T,
                                     choices = c("YES","NO"),selected = "NO"),
                        radioButtons("CellIm", label = "Cell-by-cell analysis ?",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        conditionalPanel(
                          condition = "input.CellIm == 'YES'",
                          radioButtons("ExportFCS", label = "Export cell data as a fcs file ?",inline = T,
                                       choices = c("YES","NO"),selected = "YES")
                        ),
                        radioButtons("ExportResults", label = "Export results as :",inline = T,
                                     choices = c("RData", "CSV"),selected = "RData"),
                        radioButtons("ExpSeg", label = "Export segmentation ?",inline = T,
                                     choices = c("YES","NO"),selected = "NO")
                      ),
                      
                      wellPanel(
                        radioButtons("BlueChannel", label = "Nucleus :",inline = T,
                                     choices = c(1,2,3),selected = 1),
                        radioButtons("GreenChannel", label = "Compartment 1 :",inline = T,
                                     choices = c(1,2,3),selected = 2),
                        radioButtons("RedChannel", label = "Compartment 2 :",inline = T,
                                     choices = c(1,2,3),selected = 3),
                        radioButtons("Cyto", label = "Cytoplasm segmentation on :",inline = T,
                                     choices = c('Compartment 1','Compartment 2','Both'),selected = 'Compartment 1'),
                        style = 'padding-left:3px;padding-right:3px;')
               ),

               column(2,align='left', offset = 0,
                      h3('Test Settings'),
                      wellPanel(   
                        uiOutput("SampPlate"),
                        uiOutput("SampTime"),
                        uiOutput("SampWell"),
                        uiOutput("SampSite")
                      ),
                      br(),
                      h3('Export Results'),
                      wellPanel(
                        strong('Save results in :'),
                        h4(''),
                        actionButton("savefolder", label = 'Select folder')
                      )
               ),
               
               conditionalPanel(
                 condition = "input.PLATEMAP == 'YES'",
                 column(5,align='center', offset = 0,
                        h3('Plate Map'),
                        dataTableOutput('PlateMap')
                 )
               )
             )
    ),
    
    tabPanel("Nucleus",
             fluidRow(
               column(3,align='center',
                      br(),
                      wellPanel(
                        radioButtons("Nucrm", label = "Remove nucleus from mask?",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        conditionalPanel(
                          condition = "input.Nucrm == 'YES'",
                          numericInput('w1OFF',label='Offset for Nucleus',value = 0.1, step=0.01)),
                        numericInput('adj.step1',label='Extrema smoothing',value = 2, step=1,min=1,max=8),
                        sliderInput("Rm1", "Adjust image:",
                                    min=0, max=1, value=c(0,1.00)),
                        (actionButton('Test1',label='TEST SETTINGS'))
                      )
                      
               ),
               
               column(6, align='center',
                      br(),
                      imageOutput("LookUp1",height=600, width=600, inline=F)
               ))
    ),
    
    tabPanel("Compartment 1",
             
             fluidRow(
               column(3,align='center',
                      br(),
                      wellPanel(
                        radioButtons("auto2", label = "Automated segmentation :",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        conditionalPanel(
                          condition = "input.auto2 == 'NO'",
                          numericInput('w2OFF',label='Offset for Comparment1',value = 0.1, step=0.01)),
                        numericInput('TopSize2',label='Top Hat filter size :',value = 29, step=2),
                        numericInput('adj.step2',label='Extrema smoothing',value = 2, step=1,min=1,max=8),
                        sliderInput("Rm2", "Adjust image:", 
                                    min=0, max=1, value=c(0,1.00)),
                        (actionButton('Test2',label='TEST SETTINGS'))
                      )),
               column(6, align='center',
                      br(),
                      imageOutput("LookUp2",height=600, width=600, inline=F)
               ))
    ),
    
    tabPanel("Compartment 2",
             
             fluidRow(
               column(3,align='center',
                      br(),
                      wellPanel(
                        radioButtons("auto3", label = "Automated segmentation :",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        conditionalPanel(
                          condition = "input.auto3 == 'NO'",
                          numericInput('w3OFF',label='Offset for Comparment2',value = 0.1, step=0.01)),
                        numericInput('TopSize3',label='Top Hat filter size :',value = 29, step=2),
                        numericInput('adj.step3',label='Extrema smoothing',value = 2, step=1,min=1,max=8),
                        sliderInput("Rm3", "Adjust image:", 
                                    min=0, max=1, value=c(0,1.00)),
                        (actionButton('Test3',label='TEST SETTINGS'))
                      )),
               column(6, align='center',
                      br(),
                      imageOutput("LookUp3",height=600, width=600, inline=F)
               ))       
    ),
    
    tabPanel("Merge",align='center',
             fluidRow(
               column(3,align='center',
                      br(style='font-size:400%'),
                      strong('Pearson coefficient', style = 'font-size:150%'),
                      uiOutput("SampPCC", style = 'font-size:150%'),
                      br(style = 'font-size:800%'),
                      actionButton("launcher", label='LAUNCH ANALYSIS',style='padding:10px; font-size:150%;align:center',width=290),
                      h3(''),
                      actionButton("stop", label='QUIT ANALYSIS',style='padding:10px; font-size:150%;align:center',width=290)
                      
               ),
               column(6, align = 'center',
                      br(),
                      imageOutput("LookUp4",height=600, width=600, inline=F)
               )
             )
    )
  )
)

