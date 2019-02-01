library(shiny)
library(shinyjs)
library(MetaxpR)

library(gtools)
library(rChoiceDialogs)
library(RODBC)
library(EBImage)
library(MorphR)

library(doParallel)
library(doSNOW)
library(foreach)

server = function(input, output, session) {
  
  session$onSessionEnded(function(){
    stopApp()
  })
  
  OS = Sys.info()[['sysname']]
  
  #===================================================================================================================================================
  ## Buttons activation logistics
  
  observe({
    if (is.null(ConInf$Welldat)){
      shinyjs::disable(id='ColSet')
    }else{
      shinyjs::enable(id='ColSet')
    }
  })
  
  #===================================================================================================================================================
  ## PlateMap information
  ConInf = reactiveValues(DrugID = data.frame(PlateID = rep('PlateID',15),WellID=paste0('Well',1:15), Drug = paste('Drug',c(1:15)),
                                              Conc = c(sample.int(100,size=15)),Unit = rep('uM',15)),
                          Welldat = NULL, DB=NULL)
  Plates <- reactiveValues()
  TestImage = reactiveValues(Im = NULL, size = NULL, status = 1)
  
  output$PlateMap = renderDataTable(ConInf$DrugID, options = list(pageLength = 15,lengthMenu = c(15,20,25)))
  AdjIm = reactiveValues(Auto1 = c(0,1), Auto2 = c(0,1), Auto3 = c(0,1))
  Settings.status = reactiveValues(text = "No error detected", color = '#76EE00')
  #
  
  ## Plate selection
  observe({
    if(input$UseSQL == 'YES'){
      cnx = names(odbcDataSources(type='system'))
      if(length(cnx!=0)){
        output$ODBC = renderUI({
          selectizeInput("SERVER","Database :", choices = lapply(names(odbcDataSources(type='system')),function(x)x), 
                         multiple = F, selected=tail(names(odbcDataSources(type='system')),n=1))
        })
        observeEvent(input$SERVER,{
          ConInf$DB = GetMDCInfo(input$SERVER,Unix.diff = c('//H','/media/H')) #Unix.diff can be replaced with your own config
          output$PlateIDs = renderUI({
            Plates = unique(as.numeric((ConInf$DB)$PlateID))
            selectizeInput("SPlates","Plate selection :", choices = Plates[order(Plates, decreasing = T)], multiple = input$MulPlates=='YES', selected=Plates[length(Plates[!is.na(Plates)])])
          })
        })
      }else{
        output$ODBC = renderUI({
          textOutput("SERVER")
        })
        output$SERVER = renderText("No database detected")
      }
    }})
  
  observeEvent(input$folder,{
    if(OS == 'Darwin'){
      Plates$location = normalizePath(path.expand(as.character(system("osascript -e 'set thisPOSIXPath to (the POSIX path of (choose folder with prompt \"Select image location\"))'", intern = T))))
    }else{
      Plates$location = normalizePath(path.expand(rchoose.dir(caption = "Select image location")), winslash = '/')
    }
    #
    if(length(Plates$location) == 0){ Plates$location = normalizePath(path.expand(getwd()), winslash = '/') }
    paths = list.dirs(Plates$location, recursive = T, full.names = T)
    Plates$folders = na.omit(basename(paths))
    Plates$paths = paths[!is.na(basename(paths))]
    output$PlateIDs = renderUI({
      selectizeInput("SPlates","Plate selection :", choices = Plates$folders, multiple = input$MulPlates == 'YES', selected=Plates$folders[1])
    })
  })
  
  
  observeEvent(input$savefolder,{
    if(OS == 'Darwin'){
      Plates$resultsloc = paste(normalizePath(path.expand(as.character(system("osascript -e 'set thisPOSIXPath to (the POSIX path of (choose folder with prompt \"Select location for saving your results\"))'", intern = T)))), 'Results', sep = '/')
    }else{
      Plates$resultsloc = paste(normalizePath(path.expand(rchoose.dir(caption = "Select location for saving your results")), winslash = '/'), 'Results', sep = '/')
    }
    
    if(length(Plates$resultsloc) == 0){
      Plates$resultsloc = paste(normalizePath(path.expand(getwd()), winslash = '/'), 'Results', sep = '/')
    }
    updateTextInput(session, inputId = 'savefolder.str', value = Plates$resultsloc)
  })
  
  ## Import images------------------------------------------------
  
  observeEvent(input$ImpImg,{
    Plates$PlateIDs = input$SPlates
    withProgress({
      
      setProgress(message = "Creating clusters...")
      
      try({
        ncores = detectCores()-1
        if(ncores>length(input$SPlates)){
          ncores=length(input$SPlates)
        }
        cl = parallel::makeCluster(ncores)
        invisible(clusterEvalQ(cl, c(library(MetaxpR),library(gtools),library(reshape2),library(ColocalizR),library(doParallel))))
        TimeCourse = (input$TimeCourse == 'YES')
        PlateIDs = input$SPlates
        
        setProgress(message = "Loading images info...")
        if(input$UseSQL=='YES'){
          SERVER = input$SERVER
          DB = ConInf$DB
          ConInf$Welldat = do.call('smartbind',parLapply(cl=cl, PlateIDs, function(x) getImInfo(as.numeric(x), SQL.use = T, SERVER = SERVER, DB = DB, TimeCourse = TimeCourse)))
        }else{
          folders = Plates$folders
          paths = Plates$paths
          ConInf$Welldat = do.call('smartbind',parLapply(cl=cl,PlateIDs, function(x) getImInfo(x, SQL.use = F, PlateLoc = paths[grep(x, folders, fixed = T)], TimeCourse = TimeCourse)))
        }
        stopCluster(cl);rm(cl)
      })
      setProgress(message = "Images loaded!")
    })
  })
  
  ## Import plateMaps--------------------------------------------
  
  observe({
    PM = data.frame()
    if(input$PLATEMAP == 'YES'){
      platemaps <- reactive({
        platemaps <- input$platemaps
        platemaps$datapath <- gsub("\\\\", "/", platemaps$datapath)
        platemaps
      })
      
      if(is.null(input$platemaps)) return(NULL)
      PM = do.call('smartbind', lapply(platemaps()$datapath, function(x) read.csv(x, header = T, sep=',', stringsAsFactors = T)))
      
    }
    if(nrow(PM)!=0){
      ConInf$DrugID <- PM
    }
  })
  
  #Initialize test conditions--------------------------------
  observe({
    output$SampPlate = renderUI({
      Plate1 = (ConInf$Welldat)$PlateID[1]
      textInput("SampPlate1", "Plate :", Plate1)
    })
    output$SampTime = renderUI({
      Time1 = 1
      textInput("SampTime1", "Time point :", Time1)
    })
    output$SampWell = renderUI({
      Well1 = (ConInf$Welldat)$Well[1]
      textInput("SampWell1", "Well :", Well1)
    })
    output$SampSite = renderUI({
      Site1 = 1
      numericInput("SampSite1", "Site :", Site1)
    })
  })
  observe({  
    if(!is.null(ConInf$Welldat)){
      TestImage$Im = as.character(ConInf$Welldat$MyIm[which(ConInf$Welldat$PlateID == input$SampPlate1 & ConInf$Welldat$TimePoint ==  input$SampTime1 &
                                                              ConInf$Welldat$Well == input$SampWell1 & ConInf$Welldat$Site == input$SampSite1 & ConInf$Welldat$Channel == 'w1')])
      
      if(length(TestImage$Im)==0){
        TestImage$status = 0
      }else{
        TestImage$status = 1
      }
    }
  })
  
  #----------------------------------------------------------------
  ##CPU threads optimization
  
  observeEvent(input$CPU.OP,{
    withProgress({
      updateNumericInput(session, inputId = 'UsedCores', value = ResConfig(input$CellIm=='YES'))
      setProgress(message='Optimized')
    })
  })
  
  ##Check if settings are valid : Basic error handling
  observe({
    if(is.null(ConInf$Welldat)){
      Settings.status$text="No images loaded"
      Settings.status$color="#FF0000"
    }else{
      if(!is.null(ConInf$Welldat)&length(unique((ConInf$Welldat)$Channel))<3 & (input$CellIm=='YES')){
        Settings.status$text="Need nucleus for cell-by-cell analysis"
        Settings.status$color="#FF0000"
      }else{
        if(!is.null(ConInf$Welldat)&length(unique((ConInf$Welldat)$Channel))<3 & (input$BlueChannel!=3)){
          Settings.status$text='Nucleus must be set to "3" when only two channels are present'
          Settings.status$color="#FF0000"  
        }else{
          if(input$MulPlates=='YES' & input$TimeCourse=='YES'){
            Settings.status$text='Cannot run a timecourse analysis with multiple plates'
            Settings.status$color="#FF0000"   
          }else{
            if(length(unique(c(input$BlueChannel,input$GreenChannel,input$RedChannel)))!=3){
              Settings.status$text='The 3 channels cannot have the same number'
              Settings.status$color="#FF0000" 
            }else{
              if(!is.null(ConInf$Welldat) & TestImage$status==0){
                Settings.status$text = "Images not found";Settings.status$color = "#FF0000"
              }else{
                if(input$savefolder==0){
                  Settings.status$text="Warning : No output folder selected"
                  Settings.status$color="FF9900"
                }else{
                  Settings.status$text = "No obvious error detected";Settings.status$color = '#76EE00'
                }
              }
            }
          }
        }
      }
    }
  })
  output$SetStatus = renderText({paste(paste0('<font color=\"',Settings.status$color,'\"><b>'), Settings.status$text, "</b></font>") })
  #
  
  ##Initialize image signal ranges
  sapply(1:3,function(x)eval(parse(text=sprintf("output$hRm%d = renderUI({
    sliderInput('Rm%d', 'Adjust image', min=0, max=1, step=0.01, value=AdjIm$Auto%d)})",x,x,x))))
  
  #==========================================================================================================================================================
  ## Image display for testing channel parameters
  
  Thumb = readImage('Thumb.jpg')
  ThumbIm = reactiveValues(I = c(lapply(1:5,function(x)Thumb),c(0,0)))
  #
  
  observeEvent(input$Test1 | input$Test2 | input$Test3 | input$Test4, {
    if(input$Test1==0 && input$Test2==0 && input$Test3==0 && input$Test4==0){
      return()
    }
    if(is.null(TestImage$Im) | length(TestImage$Im)==0){
      return()
    }else{
      invisible(sapply(c("button","input"), function(x) shinyjs::disable(selector = x)))
      withProgress({
        setProgress(message='Segmenting images...')
        isolate({
          TempI = coloc.Sgl(MyImCl = ConInf$Welldat, Plate = input$SampPlate1, Time = input$SampTime1, Well= input$SampWell1, Site = input$SampSite1, Blue = as.numeric(input$BlueChannel), Green = as.numeric(input$GreenChannel), Red = as.numeric(input$RedChannel), auto2 = (input$auto2 == 'YES'), 
                            auto3 = (input$auto3 == 'YES'), Cyto = input$Cyto, Nuc.rm = (input$Nucrm == 'YES'), TopSize2 = input$TopSize2, TopSize3 = input$TopSize3,  w1OFF = input$w1OFF,w2OFF = input$w2OFF,w3OFF = input$w3OFF, Nuc.denoising = (input$Denoising=='YES'), RO.size = input$RO.size,  
                            NucWindow = as.numeric(input$NucWindow), SegCyto.met = input$SegCytoMet, CytoOFF = input$CytoOFF, CytoWindow = as.numeric(input$CytoWindow), FullIm = T, TEST=T, getCell = (input$CellIm == 'YES'), adj.step1 = input$adj.step1, adj.step2 = input$adj.step2, adj.step3 = input$adj.step3, 
                            adj = as.numeric(input$adj), getRange = c((input$AutoAd1=='YES'),(input$AutoAd2=='YES'),(input$AutoAd3=='YES')),Rm1 = input$Rm1, Rm2 = input$Rm2, Rm3= input$Rm3)
        })
        ThumbIm$I = TempI[c(1:5,9:10)]
        AdjIm$Auto1 = TempI[[6]];AdjIm$Auto2 = TempI[[7]];AdjIm$Auto3 = TempI[[8]]
        #
        setProgress(message='Done !')
      })
    }
    invisible(sapply(c("button","input"), function(x) shinyjs::enable(selector = x)))
  })
  #
  observe({
    sapply(1:5,function(x)eval(parse(text=sprintf("output$LookUp%d <-renderDisplay({
          display(ThumbIm$I[[%d]])})",x,x))))
  })
  
  ## Give an idea of PCC/SOC
  output$SampPCC = renderTable({
    data.frame(PCC = round(ThumbIm$I[[6]],2), SOC = round(ThumbIm$I[[7]],2))
  })
  
  #===============================================================================================================================================================
  ## Launch Image Analysis
  observeEvent(input$launcher,{
    invisible(sapply(c("button","input"), function(x) shinyjs::disable(selector = x)))
    #Parameters to be exported in foreach 
    isolate({
      MyImCl.FOR = ConInf$Welldat
      UniID = unique(MyImCl.FOR$GlobalID)
      #--
      args.coloc = list(MyImCl = MyImCl.FOR, Blue = as.numeric(input$BlueChannel), Green = as.numeric(input$GreenChannel),Red = as.numeric(input$RedChannel), auto2 = (input$auto2 == 'YES'), auto3 = (input$auto3 == 'YES'),
                        Cyto = input$Cyto, Nuc.rm = (input$Nucrm == 'YES'), TopSize2 = input$TopSize2, TopSize3 = input$TopSize3, Nuc.denoising = (input$Denoising=='YES'), RO.size = input$RO.size,  
                        NucWindow = as.numeric(input$NucWindow), SegCyto.met = input$SegCytoMet, CytoOFF = input$CytoOFF, CytoWindow = as.numeric(input$CytoWindow), TEST=F,getCell = (input$CellIm == 'YES'),
                        w1OFF = input$w1OFF,w2OFF = input$w2OFF,w3OFF = input$w3OFF, adj = as.numeric(input$adj), adj.step1 = input$adj.step1, adj.step2 = input$adj.step2, adj.step3 = input$adj.step3, add.features = (input$ExpFea == 'YES'), 
                        writeSeg = (input$ExpSeg == 'YES'), writePDF = (input$ExpPDF == 'YES'), path = as.character(input$savefolder.str), getRange = rep(TRUE,3))
      #--
      args.miscs = list(ExportResults=input$ExportResults,as.FCS=(input$ExportFCS == 'YES'))
      
    })
    #Log files
    log.sum = cbind.data.frame(NbPlates = length(unique(MyImCl.FOR$PlateID)), NucOffset = input$w1OFF, RmNucFromMask =(input$Nucrm == 'YES'), CytoSeg =  input$Cyto, AdjCyto = input$adj, AutoSeg.Cyto1 = (input$auto2 == 'YES'),
                               Cyto1Offset = ifelse((input$auto2 == 'YES'), NA, input$w2OFF), TopHatCyto1 = input$TopSize2, AutoSeg.Cyto2 = (input$auto2 == 'YES'), TopHatCyto2 = input$TopSize3,
                               Cyto2Offset = ifelse((input$auto3 == 'YES'), NA, input$w3OFF), HistoAdjust.Nuc = input$adj.step1, HistoAdjust.Cyto1 = input$adj.step2, HistoAdjust.Cyto2 = input$adj.step3)
    
    withProgress({
      
      # Creating folder architecture
      setProgress(message = "Creating folders...")
      MyImCl.FOR$paths = paste(Plates$resultsloc, MyImCl.FOR$PlateID, MyImCl.FOR$TimePoint, MyImCl.FOR$Well, sep = '/')
      invisible(sapply(unique(MyImCl.FOR$paths), function(x) dir.create(x, recursive = T)))
      
      #------------------------------------------------------------
      setProgress(message = "Creating clusters...")
      isolate({cl = makeSOCKcluster(input$UsedCores)})
      registerDoSNOW(cl)
      opts = list(progress=function(n)(setProgress(value=n)))
      #------------------------------------------------------------
      setProgress(value=1, message = "Treating images")
      t1=Sys.time()
      
      Summary = foreach(ID = UniID,j=icount(), .packages = c("EBImage","tiff","gtools","ColocalizR","shiny","MorphR","MetaxpR"),.inorder=FALSE,
                        .combine = 'smartbind',.options.snow=opts) %dopar% {
                          IDs = unlist(strsplit(ID,split='_'));names(IDs) = c('P', 'TI', 'W', 'S')
                          try({do.call(coloc.Sgl,c(args.coloc,list(Plate = IDs['P'], Time = IDs['TI'], Well = IDs['W'], Site = IDs['S'])))})
                        }
      Summary$ObjNum = as.numeric(Summary$ObjNum)
      Summary=Summary[order(Summary$PlateID, Summary$Time,Summary$SiteID,Summary$WellID),,drop=F]
      Summary$GlobalID = paste(Summary$PlateID, Summary$Time, Summary$Well, sep = '_')
      
      setProgress(message = "Exporting (can take a while)")
      
      if(args.miscs$ExportResults == 'CSV'){
        parLapply(cl = cl, unique(Summary$PlateID), function(x){
          PSummary = Summary[which(Summary$PlateID == x),]
          write.table(PSummary[!is.na(PSummary$WellID),], paste(args.coloc$path, paste0(x,'/Results.csv'), sep = '/'), sep=',',row.names=F)
        })
      }else{
        save(Summary, file = paste0(args.coloc$path,'/Results.RData'))
      }
      
      if(args.coloc$getCell & args.miscs$as.FCS){
        clusterEvalQ(cl, library("flowCore"))
        if(!input$TimeCourse == 'YES'){
          parLapply(cl = cl, unique(Summary$PlateID), function(x){ #Several plates, one timepoint
            PSummary = Summary[which(Summary$PlateID==x),]
            PSummary$PCC=(PSummary$PCC+1)*50; PSummary$ICQ=(PSummary$ICQ+0.5)*100;
            PSummary$MOC=PSummary$MOC*100 ;PSummary$SOC=PSummary$SOC*100;
            PSummary$SOCR=PSummary$SOCR*100;PSummary$SOCG=PSummary$SOCG*100;
            #
            MX2FCS(dat = PSummary[-which(PSummary$ObjNum==0),],
                   coi = setdiff(colnames(PSummary),c("ObjNum","PlateID","Time","WellID","SiteID","GlobalID")),
                   export = paste(args.coloc$path,x,Summary$Time[1],sep = '/'),wlab='WellID')
          })
        }else{
          parLapply(cl = cl,unique(Summary$Time),function(x){ #One plate, several timepoints
            TSummary = Summary[which(Summary$Time==x),]
            TSummary$PCC=(TSummary$PCC+1)*50; TSummary$ICQ=(TSummary$ICQ+0.5)*100;
            TSummary$MOC=TSummary$MOC*100 ;TSummary$SOC=TSummary$SOC*100;
            TSummary$SOCR=TSummary$SOCR*100;TSummary$SOCG=TSummary$SOCG*100;
            #
            MX2FCS(dat = TSummary[-which(TSummary$ObjNum==0),],
                   coi = setdiff(colnames(TSummary),c("ObjNum","PlateID","Time","WellID","SiteID","GlobalID")),
                   export = paste(args.coloc$path,Summary$PlateID[1],x,sep = '/'),wlab='WellID')
          })
        }
      }
      
      stopCluster(cl);rm(cl)
      setProgress(message = "End of treatments!")
      #---------------------------------------------------------
      
      t2=Sys.time()
      write.csv(cbind.data.frame(log.sum, Time = difftime(t2,t1, units = 'mins')),paste(args.coloc$path,'Parameters.csv', sep = '/'),row.names=F)
      rm(Summary)
      gc()
    },min=1,max=length(UniID))
    invisible(sapply(c("button","input"), function(x) shinyjs::enable(selector = x)))
  })
  
  #==================================================================================================================================================================
  ## Exit
  
  observeEvent(input$stop,{
    stopApp()
  })
}
