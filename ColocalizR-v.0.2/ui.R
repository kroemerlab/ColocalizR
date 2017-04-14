library(shiny)
library(tiff)
library(MiXR)
library(pbapply)
library(gtools)
library(rChoiceDialogs)
#
library("doParallel")
library("foreach")
#

ui = fluidPage(
  
  includeCSS("www/bootstrap.css"),
  includeCSS("www/Progress.css"),
  titlePanel("ColocalizR"),
  
  tabsetPanel(
    tabPanel("Plate information",
             fluidRow(
               column(2,
                      fluidRow(align='left',
                      h3('Plate information'),
                      wellPanel(
                        radioButtons("UseSQL", label = "USE SQL ?",inline = T,
                                     choices = c("YES","NO"),selected = "NO"),
                        radioButtons("MulPlates", label = "Multiple plates ?",inline = T,
                                     choices = c("YES","NO"),selected = "NO"),
                        radioButtons("TimeCourse", label = "Time course experiment ?",inline = T,
                                     choices = c("YES","NO"),selected = "NO"),
                        conditionalPanel(
                          condition = "input.UseSQL == 'YES'",
                          h4(''),
                          uiOutput("Server")
                        ),
                        conditionalPanel(
                          condition = "input.UseSQL == 'NO'",
                          strong('Image location :'),
                          h4(''),
                          actionButton("folder", label = 'Select folder')
                        ),
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
                      )),
                      fluidRow(align='center',
                      br(),
                      (actionButton("ImpImg", label='IMPORT PLATE INFO',style='padding:10px; font-size:125%;align:center',width=200)),
                      br(),
                      br(),
                      br(),
                      htmlOutput("SetStatus", inline=F)
                      )
               ),
               
               column(2,align='left', offset = 0,
                      
                      h3('Colocalization Settings'),
                      wellPanel(
                        radioButtons("CellIm", label = "Cell-by-cell analysis ?",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        conditionalPanel(
                          condition = "input.CellIm == 'YES'",
                          radioButtons("ExportFCS", label = "Export cell data as a fcs file ?",inline = T,
                                       choices = c("YES","NO"),selected = "NO")
                        ),
                        radioButtons("ExportResults", label = "Export results as ?",inline = T,
                                     choices = c("RData", "CSV"),selected = "RData"),
                        radioButtons("ExpSeg", label = "Export segmentation ?",inline = T,
                                     choices = c("YES","NO"),selected = "NO"),
                        conditionalPanel(
                          condition = "input.CellIm == 'YES'",
                          radioButtons("ExpFea", label = "Export cell features ?",inline = T,
                                       choices = c("YES","NO"),selected = "NO")
                        )
                      ),
                      
                      wellPanel(
                        fluidRow(
                        column(width=4,
                        radioButtons("BlueChannel", label = "Nucleus :",inline = F,
                                     choices = c(1,2,3),selected = 1)
                        ),
                        column(width=4,
                        radioButtons("GreenChannel", label = "Compt 1 :",inline = F,
                                     choices = c(1,2,3),selected = 2)
                        ),
                        column(width=4,
                        radioButtons("RedChannel", label = "Compt 2 :",inline = F,
                                     choices = c(1,2,3),selected = 3)
                        )),
                        radioButtons("Cyto", label = "Cytoplasm segmentation on :",inline = T,
                                     choices = c('Compt 1','Compt 2','Both'),selected = 'Compt 1'),
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
                 column(5,align='left', offset = 0,
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
                          radioButtons("Seg", label = "Segmentation method",inline = T,
                                     choices = c("Fast","Robust"),selected = "Fast"),
                          radioButtons("Denoising", label = "Denoise image ?",inline = T,
                                       choices = c("YES","NO"),selected = "NO"),
                          conditionalPanel(
                            condition = "input.Denoising == 'YES'",
                            numericInput('RO.size',label='Denoising filter size', value = 25, step = 2)
                          ),
                          numericInput("w1OFF",label='Offset for Nucleus',value = 0.1, step=0.01)),
                        numericInput('adj.step1',label='Extrema smoothing',value = 2, step=1,min=1,max=8),
                        radioButtons("AutoAd1", label = "Auto-adjust",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        uiOutput("hRm1"),
                        (actionButton('Test1',label='TEST SETTINGS'))
                      )
                      
               ),
               
               column(6, align='center',
                      br(),
                      imageOutput("LookUp1",height=600, width=600, inline=F),
                      br(),
                      sliderInput("zoom1", "Zoom", min=1, max=5, step=0.5, value=1)
               ))
    ),
    
    tabPanel("Compartment 1",
             
             fluidRow(
               column(3,align='center',
                      br(),
                      wellPanel(
                        conditionalPanel(
                          condition = "input.Cyto == 'Compt 1' | input.Cyto == 'Both'",
                          numericInput('adj1',label='Adjustment for cytoplasm segmentation',value = 1, step = 0.01, min = -2, max = 2)),
                        radioButtons("auto2", label = "Automated segmentation",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        conditionalPanel(
                          condition = "input.auto2 == 'NO'",
                          numericInput('w2OFF',label='Offset for Comparment1',value = 0.1, step=0.01)),
                        numericInput('TopSize2',label='Top Hat filter size',value = 29, step=2),
                        numericInput('adj.step2',label='Extrema smoothing',value = 2, step=1,min=1,max=8),
                        radioButtons("AutoAd2", label = "Auto-adjust",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        uiOutput("hRm2"),
                        (actionButton('Test2',label='TEST SETTINGS'))
                      )),
               column(6, align='center',
                      br(),
                      imageOutput("LookUp2",height=600, width=600, inline=F),
                      br(),
                      sliderInput("zoom2", "Zoom", min=1, max=5, step=0.5, value=1)
               ))
    ),
    
    tabPanel("Compartment 2",
             
             fluidRow(
               column(3,align='center',
                      br(),
                      wellPanel(
                        conditionalPanel(
                          condition = "input.Cyto == 'Compt 2' | input.Cyto == 'Both'",
                          numericInput('adj2',label='Adjustment for cytoplasm segmentation', value = 1, step = 0.01, min = -2, max = 2)),
                        radioButtons("auto3", label = "Automated segmentation",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        conditionalPanel(
                          condition = "input.auto3 == 'NO'",
                          numericInput('w3OFF',label='Offset for Comparment2',value = 0.1, step=0.01)),
                        numericInput('TopSize3',label='Top Hat filter size',value = 29, step=2),
                        numericInput('adj.step3',label='Extrema smoothing',value = 2, step=1,min=1,max=8),
                        radioButtons("AutoAd3", label = "Auto-adjust",inline = T,
                                     choices = c("YES","NO"),selected = "YES"),
                        uiOutput("hRm3"),
                        (actionButton('Test3',label='TEST SETTINGS'))
                      )),
               column(6, align='center',
                      br(),
                      imageOutput("LookUp3",height=600, width=600, inline=F),
                      br(),
                      sliderInput("zoom3", "Zoom", min=1, max=5, step=0.5, value=1)
               ))       
    ),
    
    tabPanel("Merge",align='center',
             fluidRow(
               column(3,align='center',
                      br(style='font-size:400%'),
                      strong('Correlation', style = 'font-size:150%'),
                      uiOutput("SampPCC", style = 'font-size:150%'),
                      br(style = 'font-size:800%'),
                      actionButton("launcher", label='LAUNCH ANALYSIS',style='padding:10px; font-size:150%;align:center',width=290),
                      h3(''),
                      actionButton("stop", label='QUIT ANALYSIS',style='padding:10px; font-size:150%;align:center',width=290)
                      
               ),
               column(6, align = 'center',
                      br(),
                      imageOutput("LookUp4",height=600, width=600, inline=F),
                      br(),
                      sliderInput("zoom4", "Zoom",  min=1, max=5, step=0.5, value=1)
               )
             )
    )
  )
)

