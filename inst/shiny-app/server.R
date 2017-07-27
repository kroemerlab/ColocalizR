library(shiny)
library(MiXR)

library(gtools)
library(rChoiceDialogs)
library(tiff)
library(EBImage)
library(RODBC)

library(doParallel)
library(doSNOW)
library(foreach)

server = function(input, output, session) {
  
  session$onSessionEnded(function(){
    stopApp()
  })
  
  OS = Sys.info()[['sysname']]
  #===================================================================================================================================================
  ## PlateMap information
  ConInf = reactiveValues(DrugID = data.frame(PlateID = rep('PlateID',15),WellID=paste0('Well',1:15), Drug = paste('Drug',c(1:15)),
                                              Conc = c(sample.int(100,size=15)),Unit = rep('uM',15)),
                          Welldat = NULL, DB=NULL)
  Plates <- reactiveValues()
  #
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
          ConInf$DB = GetMDCInfo(input$SERVER,Unix.diff = c('//H','/media/h')) #Unix.diff can be replaced with your own config
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
    Plates$folders = as.numeric(basename(paths))[!is.na(as.numeric(basename(paths)))]
    Plates$paths = paths[!is.na(as.numeric(basename(paths)))]
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
    Plates$PlateIDs = as.numeric(input$SPlates) 
    withProgress({
      
      setProgress(message = "Creating clusters...")
      ncores = detectCores()-1
      if(ncores>length(as.numeric(input$SPlates))){
        ncores=length(as.numeric(input$SPlates))
      }
      cl = parallel::makeCluster(ncores)
      invisible(clusterEvalQ(cl, c(library(MiXR),library(gtools),library(reshape2),library(ColocalizR),library(doParallel))))
      TimeCourse = (input$TimeCourse == 'YES')
      PlateIDs = as.numeric(input$SPlates)
      
      setProgress(message = "Loading images info...")
      if(input$UseSQL=='YES'){
        SERVER = input$SERVER
        DB = ConInf$DB
        ConInf$Welldat = do.call('smartbind',parLapply(cl=cl, PlateIDs, function(x) getImInfo(x, SQL.use = T, SERVER = SERVER, DB = DB, TimeCourse = TimeCourse)))
      }else{
        folders = Plates$folders
        paths = Plates$paths
        ConInf$Welldat = do.call('smartbind',parLapply(cl=cl,PlateIDs, function(x) getImInfo(x, SQL.use = F, PlateLoc = paths[grep(x, folders, fixed = T)], TimeCourse = TimeCourse)))
      }
      stopCluster(cl);rm(cl)
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
    if(input$savefolder==0){
      Settings.status$text="Warning : No output folder selected"
      Settings.status$color="FF9900"
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
              Settings.status$text = "No obvious error detected";Settings.status$color = '#76EE00'
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
  
  ##
  Adj = reactiveValues(adj_value = 1)
  observe({
    if(input$Cyto == 'Both'){ 
      observeEvent((input$adj1 | input$adj2),{
        inputs = c(input$adj1,input$adj2)
        if(all(!is.null(inputs) & !is.na(inputs))){
          if(any(inputs != Adj$adj_value)){
            Adj$adj_value = inputs[which(inputs != Adj$adj_value)]
            updateNumericInput(session, inputId = 'adj1', value = Adj$adj_value)
            updateNumericInput(session, inputId = 'adj2', value = Adj$adj_value)
          }
        }
      })
    }else if(input$Cyto == 'Compt 1'){ 
      Adj$adj_value = input$adj1
    }else if(input$Cyto == 'Compt 2'){ 
      Adj$adj_value = input$adj2
    }
  })
  
  Zoom = reactiveValues(zoom_value = 1)
  observeEvent((input$zoom1 | input$zoom2 | input$zoom3 | input$zoom4),{
    Inputs = c(input$zoom1,input$zoom2,input$zoom3,input$zoom4)
    if(any(Inputs != Zoom$zoom_value)){
      Zoom$zoom_value = unique(Inputs[which(Inputs != Zoom$zoom_value)])
      invisible(sapply(1:4,function(x)eval(parse(text=sprintf("updateSliderInput(session, 'zoom%d', value = Zoom$zoom_value,
                        min=1, max=5, step=0.5)",x)))))
    }
  })
  
  #==========================================================================================================================================================
  ## Image display for testing channel parameters
  
  Thumb = readImage('Thumb.jpg')
  ThumbIm = reactiveValues(I = c(lapply(1:4,function(x)Thumb),c(0,0)))
  #
  observeEvent({
    input$Test3
    input$Test2
    input$Test1},{
      withProgress({
        setProgress(message='Segmenting images...')
        TempI = coloc.Sgl(MyImCl = ConInf$Welldat, Plate = input$SampPlate1, Time = input$SampTime1, Well= input$SampWell1, Site = input$SampSite1, Blue = as.numeric(input$BlueChannel), Green = as.numeric(input$GreenChannel), Red = as.numeric(input$RedChannel), auto2 = (input$auto2 == 'YES'), 
                          auto3 = (input$auto3 == 'YES'), Cyto = input$Cyto, Nuc.rm = (input$Nucrm == 'YES'), TopSize2 = input$TopSize2, TopSize3 = input$TopSize3,  w1OFF = input$w1OFF,w2OFF = input$w2OFF,w3OFF = input$w3OFF, Nuc.denoising = (input$Denoising=='YES'), RO.size = input$RO.size,  
                          FullIm = T, TEST=T, getCell = (input$CellIm == 'YES'), adj.step1 = input$adj.step1, adj.step2 = input$adj.step2, adj.step3 = input$adj.step3, adj = Adj$adj_value, getRange = c((input$AutoAd1=='YES'),(input$AutoAd2=='YES'),(input$AutoAd3=='YES')), 
                          Rm1 = input$Rm1, Rm2 = input$Rm2, Rm3= input$Rm3)
        ThumbIm$I = TempI[c(1:4,8:9)]
        AdjIm$Auto1 = TempI[[5]];AdjIm$Auto2 = TempI[[6]];AdjIm$Auto3 = TempI[[7]]
        #
        setProgress(message='Done !')
      })
    })
  #
  observe({
    Center = dim(ThumbIm$I[[1]])/2
    xEx1 = ((Center[1]*(1-(1/Zoom$zoom_value))):(Center[1]*(1+(1/Zoom$zoom_value))))
    yEx1 = ((Center[2]*(1-(1/Zoom$zoom_value))):(Center[2]*(1+(1/Zoom$zoom_value))))
    
    isolate({
      sapply(1:4,function(x)eval(parse(text=sprintf("output$LookUp%d <-renderPlot({
          display((ThumbIm$I[[%d]])[xEx1,yEx1,], method='raster')
        })",x,x))))
    })
  })
  
  ## Give an idea of PCC/SOC
  output$SampPCC = renderTable({
    data.frame(PCC = round(ThumbIm$I[[5]],2), SOC = round(ThumbIm$I[[6]],2))
  })
  
  #===============================================================================================================================================================
  ## Launch Image Analysis
  observeEvent(input$launcher,{
    
    #Parameters to be exported in foreach 
    isolate({
      MyImCl.FOR = ConInf$Welldat
      #--
      UniID = unique(MyImCl.FOR$GlobalID)
      #--
      Blue.FOR = as.numeric(input$BlueChannel); Green.FOR = as.numeric(input$GreenChannel);Red.FOR = as.numeric(input$RedChannel); auto2.FOR = (input$auto2 == 'YES');auto3.FOR = (input$auto3 == 'YES');
      Cyto.FOR = input$Cyto;Nuc.rm.FOR = (input$Nucrm == 'YES'); Nucdenoising.FOR = (input$Denoising == 'YES'); ROsize.FOR = input$RO.size ; TopSize2.FOR = input$TopSize2; TopSize3.FOR = input$TopSize3;
      w1OFF.FOR = input$w1OFF;w2OFF.FOR = input$w2OFF;w3OFF.FOR = input$w3OFF; getCell.FOR = (input$CellIm == 'YES'); as.FCS.FOR = (input$ExportFCS == 'YES'); adj.FOR = Adj$adj_value; 
      adj.step1.FOR = input$adj.step1; adj.step2.FOR = input$adj.step2; adj.step3.FOR = input$adj.step3; Rm1.FOR = input$Rm1; Rm2.FOR = input$Rm2; Rm3.FOR= input$Rm3; 
      ExportResults.FOR = input$ExportResults ; ExpSeg.FOR = (input$ExpSeg == 'YES') ; ExpPDF.FOR = (input$ExpPDF == 'YES'); ExpFea.FOR = (input$ExpFea == 'YES') ; path.FOR = as.character(input$savefolder.str)
      #--
    })
    #Log files
    log.sum = cbind.data.frame(NbPlates = length(unique(MyImCl.FOR$PlateID)), NucOffset = w1OFF.FOR, NucInMask = Nuc.rm.FOR, CytoSeg = Cyto.FOR, AdjCyto = adj.FOR, AutoSeg.Cyto1 = auto2.FOR, Cyto1Offset = ifelse(auto2.FOR, NA, w2OFF.FOR), TopHatCyto1 = TopSize2.FOR,
                               AutoSeg.Cyto2 = auto2.FOR, TopHatCyto2 = TopSize3.FOR,  Cyto2Offset = ifelse(auto3.FOR, NA, w3OFF.FOR),  HistoAdjust.Nuc = adj.step1.FOR, HistoAdjust.Cyto1 = adj.step2.FOR, HistoAdjust.Cyto2 = adj.step3.FOR)
    
    withProgress({
      
      # Creating folder architecture
      setProgress(message = "Creating folders...")
      MyImCl.FOR$paths = paste(Plates$resultsloc, MyImCl.FOR$PlateID, MyImCl.FOR$TimePoint, MyImCl.FOR$Well, sep = '/')
      invisible(sapply(unique(MyImCl.FOR$paths), function(x) dir.create(x, recursive = T)))
      
      #------------------------------------------------------------
      setProgress(message = "Creating clusters...")
      cl <- makeSOCKcluster(input$UsedCores)
      registerDoSNOW(cl)
      opts <- list(progress=function(n)(setProgress(value=n)))
      #------------------------------------------------------------
      setProgress(value=1, message = "Treating images")
      t1=Sys.time()
      
      Summary = foreach(ID = UniID,j=icount(), .packages = c("EBImage","tiff","gtools","ColocalizR","shiny"),.inorder=FALSE,
                        .combine = 'smartbind',.options.snow=opts) %dopar% {
                          IDs = unlist(strsplit(ID,split='_'));names(IDs) = c('P', 'TI', 'W', 'S')
                          try({
                            coloc.Sgl(MyImCl = MyImCl.FOR, Plate = IDs['P'], Time = IDs['TI'], Well = IDs['W'], Site = IDs['S'], Blue = Blue.FOR, Green = Green.FOR,Red = Red.FOR, auto2 = auto2.FOR, auto3 = auto3.FOR,
                                      Cyto = Cyto.FOR,Nuc.rm = Nuc.rm.FOR, TopSize2 = TopSize2.FOR, TopSize3 = TopSize3.FOR, Nuc.denoising = Nucdenoising.FOR, RO.size =  ROsize.FOR, TEST=F,getCell = getCell.FOR,
                                      w1OFF = w1OFF.FOR,w2OFF = w2OFF.FOR,w3OFF = w3OFF.FOR, adj = adj.FOR, adj.step1 = adj.step1.FOR, adj.step2 = adj.step2.FOR, adj.step3 = adj.step3.FOR, add.features = ExpFea.FOR, 
                                      writeSeg = ExpSeg.FOR, writePDF = ExpPDF.FOR, path = path.FOR, getRange = rep(TRUE,3))
                          })
                        }
      colnames(Summary)[c(2:5)]=c('PlateID', 'Time','WellID','SiteID');Summary$ObjNum = as.numeric(Summary$ObjNum)
      Summary=Summary[order(Summary$PlateID, Summary$Time,Summary$SiteID,Summary$WellID),,drop=F]
      Summary$GlobalID = paste(Summary$PlateID, Summary$Time, Summary$Well, sep = '_')
      
      setProgress(message = "Exporting (can take a while)")
      
      if(ExportResults.FOR == 'CSV'){
        parLapply(cl = cl, unique(Summary$PlateID), function(x){
          PSummary = Summary[which(Summary$PlateID == x),]
          write.table(PSummary[!is.na(PSummary$WellID),], paste(path.FOR, paste0(x,'/Results.csv'), sep = '/'), sep=',',row.names=F)
        })
      }else{
        save(Summary, file = paste0(path.FOR,'/Results.RData'))
      }
      
      if(getCell.FOR==T & as.FCS.FOR==T){
        clusterEvalQ(cl, library("flowCore"))
        parLapply(cl = cl, unique(Summary$GlobalID), function(x){
          WSummary = Summary[which(Summary$GlobalID == x),]
          FlowFrame = new("flowFrame",exprs=as.matrix(WSummary[-which(WSummary$ObjNum==0),c(12:ncol(WSummary))]))
          write.FCS(FlowFrame, paste(path.FOR, paste(gsub('_','/', x), 'fcs', sep='.'), sep = '/'))
        })
        parLapply(cl = cl, unique(Summary$PlateID), function(x){
          PSummary = Summary[which(Summary$PlateID == x),]
          FlowFrame = new("flowFrame",exprs=as.matrix(PSummary[-which(PSummary$ObjNum==0),c(12:ncol(PSummary))]))
          write.FCS(FlowFrame, paste(path.FOR, paste0(x, '/GlobalFCS.fcs'), sep = '/'))
        })
      }
      
      stopCluster(cl);rm(cl)
      setProgress(message = "End of treatments!")
      
      #---------------------------------------------------------
      
      t2=Sys.time()
      write.csv(cbind.data.frame(log.sum, Time = difftime(t2,t1, units = 'mins')),paste(path.FOR,'Parameters.csv', sep = '/'),row.names=F)
      rm(Summary)
      gc()
    },min=1,max=length(UniID))
  })
  
  #==================================================================================================================================================================
  ## Emergency exit
  
  observeEvent(input$stop,{
    
    stopApp()
    if(length(grep('mytempdir',ls()))!=0){
      unlink(mytempdir, recursive =T, force=T)
    }
  })
}
