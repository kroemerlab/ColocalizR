jsdis = "
shinyjs.disableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
}

shinyjs.enableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.unbind('click.tab');
  tab.removeClass('disabled');
}
"

cssdis = "
.nav li a.disabled {
  background-color: #aaa !important;
  color: #333 !important;
  cursor: not-allowed !important;
  border-color: #aaa !important;
}"

ui = fluidPage(
  shinyjs::useShinyjs(),
  shinyjs::extendShinyjs(text = jsdis, functions=c('disableTab','enableTab')),
  shinyjs::inlineCSS(cssdis),
  includeCSS("www/bootstrap.css"),
  titlePanel("ColocalizR"),
  
  tabsetPanel(id='tabs',
              tabPanel(title="Plate information",
                       value="Plate_information",
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
                                             uiOutput("ODBC")
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
                                                        choices = c("YES","NO"),selected = "NO"),
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
                                         (actionButton("ImpImg", label='IMPORT PLATE INFO',style='padding:10px; font-size:125%;align:center',width="100%")),
                                         br(),br(),br(),
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
                                  radioButtons("ExpPDF", label = "Export Pixel Profiling ?",inline = T,
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
                                                        choiceValues=1:3,choiceNames=paste0('w',1:3),selected = 1)
                                    ),
                                    column(width=4,
                                           radioButtons("GreenChannel", label = "Compt 1 :",inline = F,
                                                        choiceValues=1:3,choiceNames=paste0('w',1:3),selected = 2)
                                    ),
                                    column(width=4,
                                           radioButtons("RedChannel", label = "Compt 2 :",inline = F,
                                                        choiceValues=1:3,choiceNames=paste0('w',1:3),selected = 3)
                                    )))
                                ,id='ColSet'),
                         
                         column(2,align='center', offset = 0,
                                h3('Output Settings'),
                                wellPanel(
                                  strong('Save results in :'),
                                  h4(''),
                                  actionButton("savefolder", label = 'Select folder'),
                                  textInput("savefolder.str", label ='', value = paste0(getwd(),'/Results'))
                                  ,id='OutSet'),
                                br(),br(),br(),
                                h3('Hardware Settings'),
                                wellPanel(
                                  numericInput("UsedCores","CPU threads:",value=parallel::detectCores()),
                                  actionButton("CPU.OP", label='Optimize',style='padding:10px; font-size:125%;align:center',width="100%")
                                )
                         ),

                         column(2,align='left', offset = 0,
                                h3(id='Testid','Test Settings'),
                                wellPanel(   
                                  uiOutput("SampPlate"),
                                  uiOutput("SampTime"),
                                  uiOutput("SampWell"),
                                  uiOutput("SampSite")
                                  ,id='TestSet')
                                ),
                         
                         conditionalPanel(
                           condition = "input.PLATEMAP == 'YES'",
                           column(3,align='left', offset = 0,
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
                                  radioButtons("Nucrm", label = "Remove nucleus from mask ?",inline = T,
                                               choices = c("YES","NO"),selected = "YES"),
                                  
                                  conditionalPanel(
                                    condition = "input.Nucrm == 'YES' | input.CellIm == 'YES'",                        
                                    radioButtons("Seg", label = "Segmentation method :",inline = T,
                                                 choices = c("Fast","Robust"),selected = "Fast"),
                                    radioButtons("Denoising", label = "Denoise image ?",inline = T,
                                                 choices = c("YES","NO"),selected = "NO"),
                                    numericInput("NucWindow",label='Window size for Nucleus',value = 50, step=5,min=1),
                                    numericInput("w1OFF",label='Offset for Nucleus',value = 0.1, step=0.01)),
                                  conditionalPanel(
                                    condition = "input.Denoising == 'YES'",
                                    numericInput('RO.size',label='Denoising filter size', value = 9, step = 2, min = 3)
                                  ),
                                  numericInput('adj.step1',label='Extrema smoothing',value = 2, step=1,min=1,max=8),
                                  radioButtons("AutoAd1", label = "Auto-adjust",inline = T,
                                               choices = c("YES","NO"),selected = "YES"),
                                  uiOutput("hRm1"),
                                  (actionButton('Test1',label='TEST SETTINGS'))
                                )
                                
                         ),
                         
                         column(6, align='center',
                                br(),
                                EBImage::displayOutput("LookUp1",width='600px',height='600px'),
                                br(),br(),
                                htmlOutput("NucStatus", inline=F)
                         ))
                       
              ),
              
              tabPanel("Cytoplasm",
                       fluidRow(
                         column(3,align='center',
                                br(),
                                wellPanel(
                                  radioButtons("Cyto", label = "Cytoplasm segmentation on :",inline = T,
                                               choices = c('Compt 1','Compt 2','Both'),selected = 'Compt 1'),
                                  radioButtons("SegCytoMet", label = "Cytoplasm segmentation method :",inline = T,
                                               choices = c('Automated','Adaptive','Both'),selected = 'Automated'),
                                  conditionalPanel(
                                    condition = "input.SegCytoMet == 'Automated' | input.SegCytoMet == 'Both'" ,
                                    numericInput('adj',label='Adjustment for cytoplasm segmentation',value = 1, step = 0.01, min = -2, max = 2)
                                  ),
                                  conditionalPanel(
                                    condition = "input.SegCytoMet == 'Adaptive' | input.SegCytoMet == 'Both'" ,
                                    numericInput("CytoWindow","Window size for Cytoplasm:", value = 50, step = 5, min = 0),
                                    numericInput('CytoOFF',label='Offset for Cytoplasm:',value = 0.1, step=0.01)
                                  ),
                                  (actionButton('Test2',label='TEST SETTINGS'))
                                )
                         ),
                         
                         column(6, align='center',
                                br(),
                                EBImage::displayOutput("LookUp2",width='600px',height='600px'),
                                br(),br(),
                                htmlOutput("CytoStatus", inline=F)
                         ))
                       
              ),
              
              tabPanel(title="Compartment 1",
                       value="Compartment_1",
                       fluidRow(
                         column(3,align='center',
                                br(),
                                wellPanel(
                                  radioButtons("auto2", label = "Automated segmentation :",inline = T,
                                               choices = c("YES","NO"),selected = "YES"),
                                  conditionalPanel(
                                    condition = "input.auto2 == 'NO'",
                                    numericInput('w2OFF',label='Offset for Comparment1',value = 0.1, step=0.01)),
                                  numericInput('TopSize2',label='Top Hat filter size',value = 35, step=2),
                                  numericInput('adj.step2',label='Extrema smoothing',value = 6, step=1,min=1,max=8),
                                  radioButtons("AutoAd2", label = "Auto-adjust",inline = T,
                                               choices = c("YES","NO"),selected = "YES"),
                                  uiOutput("hRm2"),
                                  (actionButton('Test3',label='TEST SETTINGS'))
                                )),
                         column(6, align='center',
                                br(),
                                EBImage::displayOutput("LookUp3",width='600px',height='600px'),
                                br(),br(),
                                htmlOutput("Cpt1Status", inline=F)
                         ))
              ),
              
              tabPanel(title="Compartment 2",
                       value="Compartment_2",
                       fluidRow(
                         column(3,align='center',
                                br(),
                                wellPanel(
                                  radioButtons("auto3", label = "Automated segmentation :",inline = T,
                                               choices = c("YES","NO"),selected = "YES"),
                                  conditionalPanel(
                                    condition = "input.auto3 == 'NO'",
                                    numericInput('w3OFF',label='Offset for Comparment2',value = 0.1, step=0.01)),
                                  numericInput('TopSize3',label='Top Hat filter size',value = 35, step=2),
                                  numericInput('adj.step3',label='Extrema smoothing',value = 6, step=1,min=1,max=8),
                                  radioButtons("AutoAd3", label = "Auto-adjust",inline = T,
                                               choices = c("YES","NO"),selected = "YES"),
                                  uiOutput("hRm3"),
                                  (actionButton('Test4',label='TEST SETTINGS'))
                                )),
                         column(6, align='center',
                                br(),
                                EBImage::displayOutput("LookUp4",width='600px',height='600px'),
                                br(),br(),
                                htmlOutput("Cpt2Status", inline=F)
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
                                EBImage::displayOutput("LookUp5",width='600px',height='600px')
                         )
                       )
              )
  )
)

