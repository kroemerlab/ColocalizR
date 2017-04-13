
server = function(input, output, session) {
  
  session$onSessionEnded(function(){
    stopApp()
  })
  
  ## Load functions
  MyFuncs = list.files('Funcs',full.names = T)
  for(i in 1:length(MyFuncs)){
    source(MyFuncs[i])
  }
  
  OS = Sys.info()[['sysname']]
  #===================================================================================================================================================
  
  ## PlateMap information
  
  ConInf <- reactiveValues(DrugID = data.frame(PlateID = rep('PlateID',15),WellID=paste0('Well',1:15), Drug = paste('Drug',c(1:15)),
                                               Conc = c(sample.int(100,size=15)),Unit = rep('uM',15)),Welldat = NULL)
  #
  output$PlateMap = renderDataTable(ConInf$DrugID, options = list(pageLength = 15,lengthMenu = c(15,20,25)))
  #
  AdjIm = reactiveValues(Auto1 = c(0,1), Auto2 = c(0,1), Auto3 = c(0,1))
  #
  Settings.status = reactiveValues(text = "No error detected", color = '#76EE00')
  
  ## Plate selection
  Plates <- reactiveValues()
  
  observe({
    if(input$MulPlates == 'YES'){
      multi = T
    }else{
      multi = F
    }
    
    if(input$UseSQL == 'YES'){
    output$Server = renderUI({
      selectizeInput("SERVER","Server to use :", choices = list("YODA-SERVER" = "YODA-SERVER","VIDEO-SERVER" = "VIDEO-SERVER","HTS-SERVER" = "HTS-SERVER","WORKING-SERVER" = "WORKING-SERVER"), 
                     multiple = F, selected="YODA-SERVER")
    })
    output$PlateIDs = renderUI({
      Plates = unique(as.numeric(GetMDCInfo(input$SERVER)$PlateID))
      selectizeInput("SPlates","Plate selection :", choices = Plates[order(Plates, decreasing = T)], multiple = multi, selected=Plates[length(Plates[!is.na(Plates)])])
    })
    }else{
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
          selectizeInput("SPlates","Plate selection :", choices = Plates$folders, multiple = multi, selected=Plates$folders[1])
        })
      })
    }
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
  })
  
  ## Import images
  
  observeEvent(input$ImpImg,{
    withProgress({
      setProgress(message = "Loading images info...")

      cl = makeCluster(detectCores())
      clusterExport(cl, c('getImInfo'))
      clusterEvalQ(cl, c(library(MiXR),library(pbapply),library(gtools),library(reshape2),library(naturalsort),library(doParallel)))
      Plates$PlateIDs = as.numeric(input$SPlates)
      
      PlateIDs = as.numeric(input$SPlates)
      TimeCourse = (input$TimeCourse == 'YES')
        if(input$UseSQL == 'YES'){
          SERVER = input$SERVER
          imInfo = do.call('smartbind',parLapply(cl = cl, PlateIDs, function(x) getImInfo(x, SQL.use = T, SERVER = SERVER, TimeCourse = TimeCourse)))
        }else if(input$UseSQL == 'NO'){
          folders = Plates$folders
          paths = Plates$paths
          imInfo = do.call('smartbind',parLapply(cl = cl, PlateIDs, function(x) getImInfo(x, SQL.use = F, PlateLoc = paths[grep(x, folders, fixed = T)], TimeCourse = TimeCourse)))
        }
      stopCluster(cl)

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
      ConInf$Welldat <- imInfo
      ConInf$DrugID <- PM
    })

    #Initialize test conditions
    output$SampPlate = renderUI({
      Plate1 = Plates$PlateIDs[1]
      textInput("SampPlate1", "Plate :", Plate1)
    })
    
    output$SampTime = renderUI({
      Time1 = 1
      textInput("SampTime1", "Time point :", Time1)
    })
    
    output$SampWell = renderUI({
      Well1 = imInfo$Well[1]
      textInput("SampWell1", "Well :", Well1)
    })
    
    output$SampSite = renderUI({
      Site1 = 1
      numericInput("SampSite1", "Site :", Site1)
    }) 
    
  })

  
  ##Check if settings are valid : Basic error handling
  observe({
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
            Settings.status$text = "No error detected";Settings.status$color = '#76EE00'
          }
        }
      }
    }
  })
  output$SetStatus = renderText({paste(paste0('<font color=\"',Settings.status$color,'\"><b>'), Settings.status$text, "</b></font>") })
  #
  
  ##Initialize image signal ranges
  output$hRm1 = renderUI({
    sliderInput("Rm1", "Adjust image", 
                min=0, max=1, step=0.01, value=AdjIm$Auto1)
  }) 
  output$hRm2 = renderUI({
    sliderInput("Rm2", "Adjust image", 
                min=0, max=1, step=0.01, value=AdjIm$Auto2)
  }) 
  output$hRm3 = renderUI({
    sliderInput("Rm3", "Adjust image", 
                min=0, max=1, step=0.01, value=AdjIm$Auto3)
  }) 
  
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
      updateSliderInput(session, "zoom1", value = Zoom$zoom_value,
                        min=1, max=5, step=0.5)
      updateSliderInput(session, "zoom2", value = Zoom$zoom_value,
                        min=1, max=5, step=0.5)
      updateSliderInput(session, "zoom3", value = Zoom$zoom_value,
                        min=1, max=5, step=0.5)
      updateSliderInput(session, "zoom4", value = Zoom$zoom_value,
                        min=1, max=5, step=0.5)
    }
  })

  #===================================================================================================================================================
  
  ## Image display for testing channel parameters
  
  ThumbIm = reactiveValues(I = list(readImage('Thumb.jpg'),
                                    readImage('Thumb.jpg'),
                                    readImage('Thumb.jpg'),
                                    readImage('Thumb.jpg'),
                                    0,0))
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
        #
        AdjIm$Auto1 = TempI[[5]]
        AdjIm$Auto2 = TempI[[6]]
        AdjIm$Auto3 = TempI[[7]]
        #
        setProgress(message='Done !')
      })
    })
  #
  observe({
    Center = dim(ThumbIm$I[[1]])/2
    xEx1 = ((Center[1]*(1-(1/Zoom$zoom_value))):(Center[1]*(1+(1/Zoom$zoom_value))))
    yEx1 = ((Center[2]*(1-(1/Zoom$zoom_value))):(Center[2]*(1+(1/Zoom$zoom_value))))

    output$LookUp1 <-renderPlot({
      display((ThumbIm$I[[1]])[xEx1,yEx1,], method='raster')
    })
    output$LookUp2 <-renderPlot({
      display((ThumbIm$I[[2]])[xEx1,yEx1,], method='raster')
    })
    output$LookUp3 <-renderPlot({
      display((ThumbIm$I[[3]])[xEx1,yEx1,], method='raster')
    })
    output$LookUp4 <-renderPlot({
      display((ThumbIm$I[[4]])[xEx1,yEx1,], method='raster')
    })

  })

  
  ## Give an idea of PCC
  
  output$SampPCC = renderTable({
    data.frame(PCC = round(ThumbIm$I[[5]],2), SOC =round(ThumbIm$I[[6]],2))
  })
  
  #=====================================================================================================================================================================
  ## Launch Image Analysis
  
  observeEvent(input$launcher,{
    
    #Parameters to be exported in foreach  
    MyImCl.FOR = ConInf$Welldat
    #--
    Plate = unique(MyImCl.FOR$PlateID)
    Time = unique(MyImCl.FOR$TimePoint)
    UniWell = unique(MyImCl.FOR$Well)
    UniSite = unique(MyImCl.FOR$Site)
    UniID = unique(MyImCl.FOR$GlobalID)
    #--
    Blue.FOR = as.numeric(input$BlueChannel); Green.FOR = as.numeric(input$GreenChannel);Red.FOR = as.numeric(input$RedChannel); auto2.FOR = (input$auto2 == 'YES');auto3.FOR = (input$auto3 == 'YES');
    Cyto.FOR = input$Cyto;Nuc.rm.FOR = (input$Nucrm == 'YES'); Nucdenoising.FOR = (input$Denoising == 'YES'); ROsize.FOR = input$RO.size ; TopSize2.FOR = input$TopSize2; TopSize3.FOR = input$TopSize3;
    w1OFF.FOR = input$w1OFF;w2OFF.FOR = input$w2OFF;w3OFF.FOR = input$w3OFF; getCell.FOR = (input$CellIm == 'YES'); as.FCS.FOR = (input$ExportFCS == 'YES'); adj.FOR = Adj$adj_value; 
    adj.step1.FOR = input$adj.step1; adj.step2.FOR = input$adj.step2; adj.step3.FOR = input$adj.step3; Rm1.FOR = input$Rm1; Rm2.FOR = input$Rm2; Rm3.FOR= input$Rm3; 
    ExportResults.FOR = input$ExportResults ; ExpSeg.FOR = (input$ExpSeg == 'YES') ; ExpFea.FOR = (input$ExpFea == 'YES') ; path.FOR = as.character(Plates$resultsloc)
    #--
    
    #Log file
    log.file = cbind.data.frame(NbPlates = length(Plate), NucOffset = w1OFF.FOR, NucInMask = Nuc.rm.FOR, CytoSeg = Cyto.FOR, AdjCyto = adj.FOR, AutoSeg.Cyto1 = auto2.FOR, Cyto1Offset = ifelse(auto2.FOR, NA, w2OFF.FOR), TopHatCyto1 = TopSize2.FOR,
                                AutoSeg.Cyto2 = auto2.FOR, TopHatCyto2 = TopSize3.FOR,  Cyto2Offset = ifelse(auto3.FOR, NA, w3OFF.FOR),  HistoAdjust.Nuc = adj.step1.FOR, HistoAdjust.Cyto1 = adj.step2.FOR, HistoAdjust.Cyto2 = adj.step3.FOR)
    
    withProgress({
      
      setProgress(message = "Optimizing Performance...")#########
      OS = Sys.info()[['sysname']]
      if(OS == 'Windows'){
        RAM = shell('wmic OS get FreePhysicalMemory /Value',intern=T)
        RAM = RAM[grep('FreePhysicalMemory', RAM)]
        RAM = as.numeric(gsub('FreePhysicalMemory=','',RAM))
      }else if(OS == 'Linux'){
        RAM = as.numeric(system(" awk '/MemFree/ {print $2}' /proc/meminfo", intern=T))
      }else if(OS == 'Darwin'){
        RAM = as.numeric(unlist(strsplit(system("sysctl hw.memsize", intern = T), split = ' '))[2])/1e3
      }
      Cores = detectCores()
      if(getCell.FOR==T){
        Core2RAM = 1.5e06 # Assuming one core uses 1.5GB of RAM
      }else{
        Core2RAM = 8e05  # Assuming one core uses 0.8GB of RAM
      }
      MaxCores = floor(RAM/Core2RAM) 
      
      if(MaxCores>=Cores){
        UsedCores = Cores
      }else{
        UsedCores = MaxCores
      }
      
      # Creating folder architecture
      MyImCl.FOR$paths = paste(Plates$resultsloc, MyImCl.FOR$PlateID, MyImCl.FOR$TimePoint, MyImCl.FOR$Well, sep = '/')
      folders = lapply(unique(MyImCl.FOR$paths), function(x) dir.create(x, recursive = T))
      rm(folders)
      #------------------------------------------------------------
      setProgress(message = "Creating clusters...")
      
      if(OS == 'Windows'){
        cl <- makeCluster(UsedCores)
      }else{
        cl <- makeCluster(UsedCores, type='FORK')
      }
      registerDoParallel(cl, cores = UsedCores)
      opts =list(chunkSize=1)
      #------------------------------------------------------------
      setProgress(message = "Treating images")
      t1=Sys.time()

      Summary <- foreach(ID = UniID,j=icount(), .packages = c("EBImage","tiff","reshape2","gtools"),.inorder=FALSE,
                           .combine = 'smartbind',.export=c('coloc.Sgl', 'ReconsOpening', 'geodilate', 'AutoAdjust')) %dopar% {
                             IDs = reshape2::colsplit(ID, pattern = '_', names = c('P', 'TI', 'W', 'S'))
                             try({
                               coloc.Sgl(MyImCl = MyImCl.FOR, Plate = IDs$P, Time = IDs$TI, Well = IDs$W, Site = IDs$S, Blue = Blue.FOR, Green = Green.FOR,Red = Red.FOR, auto2 = auto2.FOR, auto3 = auto3.FOR,
                                         Cyto = Cyto.FOR,Nuc.rm = Nuc.rm.FOR, TopSize2 = TopSize2.FOR, TopSize3 = TopSize3.FOR, Nuc.denoising = Nucdenoising.FOR, RO.size =  ROsize.FOR, TEST=F,getCell = getCell.FOR,
                                         w1OFF = w1OFF.FOR,w2OFF = w2OFF.FOR,w3OFF = w3OFF.FOR, adj = adj.FOR, adj.step1 = adj.step1.FOR, adj.step2 = adj.step2.FOR, adj.step3 = adj.step3.FOR, add.features = ExpFea.FOR, 
                                         writeSeg = ExpSeg.FOR, path = path.FOR, getRange = c(T,T,T))
                             })
                           }
      colnames(Summary)[c(2:5)]=c('PlateID', 'Time','WellID','SiteID')
      Summary$ObjNum = as.numeric(Summary$ObjNum)
      Summary = Summary[order(Summary$PlateID, Summary$Time,Summary$SiteID,Summary$WellID),,drop=F]
      Summary$GlobalID = paste(Summary$PlateID, Summary$Time, Summary$Well, sep = '_')
      
      ####
      parLapply(cl = cl, unique(Summary$PlateID), function(x){
        PSummary = Summary[which(Summary$PlateID == x),]
        if(any(!is.na(Summary$PCC))){
          pdf(paste(path.FOR, x,'BoxPlot.pdf', sep = '/'),w=12,h=7)
          bp = boxplot(as.numeric(PCC) ~ GlobalID, data= PSummary[!is.na(PSummary$WellID),], ylim=c(-1,1), xaxt='n', outline=F, col='green',xpd=T,xaxs='i',yaxs='i',bty='n')
          stripchart(as.numeric(PCC) ~ GlobalID, data= PSummary,vertical=T,pch=1,add=T,cex=0.3)
          text(c(1:length(bp$n)), rep(-1,length(bp$n)), gsub('_|NA',' ',bp$names),cex=0.7, srt=45,adj=c(1.1,1.1),xpd=T)
          dev.off()
        }

        if(ExportResults.FOR == 'CSV'){
          write.table(PSummary[!is.na(PSummary$WellID),], paste(path.FOR, paste0(x,'/Results.csv'), sep = '/'), sep=',',row.names=F)
        }else{
          ExpSummary = PSummary[!is.na(PSummary$WellID),]
          save(ExpSummary, file = paste(path.FOR,paste0(x,'/Results.RData'), sep = '/'))
        }
      })
      
       if(getCell.FOR==T & as.FCS.FOR==T){
         clusterEvalQ(cl, library("flowCore"))
         parLapply(cl = cl, unique(Summary$GlobalID), function(x){
                                                          WSummary = Summary[which(Summary$GlobalID == x),]
                                                          FlowFrame = new("flowFrame",exprs=as.matrix(WSummary[-which(WSummary$ObjNum==0),c(12:17)]))
                                                          write.FCS(FlowFrame, paste(path.FOR, paste(gsub('_','/', x), 'fcs', sep='.'), sep = '/'))
                                                      })
         parLapply(cl = cl, unique(Summary$PlateID), function(x){
                                                          PSummary = Summary[which(Summary$PlateID == x),]
                                                          FlowFrame = new("flowFrame",exprs=as.matrix(PSummary[-which(PSummary$ObjNum==0),c(12:17)]))
                                                          write.FCS(FlowFrame, paste(path.FOR, paste0(x, '/GlobalFCS.fcs'), sep = '/'))
                                                      })
       }
         

  
      stopCluster(cl)
      setProgress(message = "End of treatments!")

      #---------------------------------------------------------

      setProgress(message = "Exporting data...")
      t2=Sys.time()
      write.csv(cbind.data.frame(log.file, Time = difftime(t2,t1, units = 'mins')),paste(path.FOR,'log.csv', sep = '/'),sep=',',row.names=F)
      rm(Summary)
      gc()
    })
  })
  
  #==============================================================================================================================================================
  ## Emergency exit
  
  observeEvent(input$stop,{
    
    stopApp()
    if(length(grep('mytempdir',ls()))!=0){
      unlink(mytempdir, recursive =T, force=T)
    }
  })
}
