
server = function(input, output, session) {
  
  session$onSessionEnded(function(){
    stopApp()
  })
  
  ## Load functions
  MyFuncs = list.files('Funcs',full.names = T)
  for(i in 1:length(MyFuncs)){
    source(MyFuncs[i])
  }

  #===================================================================================================================================================
  
  ## PlateMap information
  
  ConInf <- reactiveValues(DrugID = data.frame(PlateID = rep('MyPlateID',15),WellID=rep('MyWellID',15), Drug = paste('Drug',c(1:15)), Conc = c(rep(10,15)),Unit = rep('uM',15)),
                           Welldat = NULL)
  
  output$PlateMap = renderDataTable(ConInf$DrugID, options = list(pageLength = 15,lengthMenu = c(15,20,25)))
  
  ## Plate selection
  Plates <- reactiveValues()
  
  observeEvent(input$folder,{
    
    Plates$location = normalizePath(path.expand(choose.dir(caption = "Select image location")), winslash = '/')
 
    #
    if(input$MulPlates == 'YES'){
      multi = T
      Plates$folders = list.dirs(Plates$location, recursive = F, full.names = F)
      Plates$paths = list.dirs(Plates$location, recursive = F, full.names = T)
    }else{
      multi = F
      Plates$paths = Plates$location
      Plates$folders = tail(unlist(strsplit(Plates$location,'/')),n=1)
      if(Plates$folders == 'NA'){
        Plates$folders = ''
      }
    }
    output$PlateIDs = renderUI({
      selectizeInput("SPlates","Plate selection :", choices = Plates$folders, multiple = multi, selected=Plates$folders[1])
    })
  })
  
  observeEvent(input$savefolder,{
    Plates$resultsloc = normalizePath(path.expand(choose.dir(caption = "Select location for saving")), winslash = '/')
  })
  
  ## Import images
  observeEvent(input$ImpImg,{
    withProgress({
      setProgress(message = "Loading images info...")
      
      Plates$PlateIDs = as.numeric(input$SPlates)
      
      imInfo = data.frame()
      PM = data.frame()
      
      for(i in 1:length(Plates$PlateIDs)){
        PlateLocation = Plates$paths[grep(Plates$PlateIDs[i], Plates$folders, fixed = T)]
        imInfo = rbind(imInfo, getImInfo(PlateID = Plates$PlateIDs[i], PlateLoc = PlateLocation, TimeCourse=(input$TimeCourse=='YES')))
      }
      
      if(input$PLATEMAP == 'YES'){
        platemaps <- reactive({
          platemaps <- input$platemaps
          platemaps$datapath <- gsub("\\\\", "/", platemaps$datapath)
          platemaps
        })
        
        if(is.null(input$platemaps)) return(NULL)
        
        for (i in 1:nrow(platemaps())){
          PM = rbind(PM, read.csv(platemaps()$datapath[i], header = T, sep=',', stringsAsFactors = T))
        }
      }
      ConInf$Welldat <- imInfo
      ConInf$DrugID <- PM
    })
    
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
    
    print(imInfo)
    
  })
  
  
  #===================================================================================================================================================
  
  ## Image display for testing channel parameters
  
  ThumbIm = reactiveValues(I = list(readImage('Thumb.jpg'),readImage('Thumb.jpg'),readImage('Thumb.jpg'),readImage('Thumb.jpg'), 0))
  
  #
  observeEvent(input$Test1,{
    
    withProgress({
      
      setProgress(message='Segmenting images...')
      
      ThumbIm$I = coloc.Sgl(MyImCl = ConInf$Welldat, Plate = input$SampPlate1, Time = input$SampTime1, Well= input$SampWell1, Site = input$SampSite1, Blue = as.numeric(input$BlueChannel), Green = as.numeric(input$GreenChannel), Red = as.numeric(input$RedChannel), auto2 = (input$auto2 == 'YES'), auto3 = (input$auto3 == 'YES'),
                            Cyto = input$Cyto, Nuc.rm = (input$Nucrm == 'YES'), TopSize2 = input$TopSize2, TopSize3 = input$TopSize3,  w1OFF = input$w1OFF,w2OFF = input$w2OFF,w3OFF = input$w3OFF, TEST=T, Rm1 = input$Rm1, Rm2 = input$Rm2, Rm3= input$Rm3,
                            W_Name = W_Name, getCell = (input$CellIm == 'YES'), adj.step1 = input$adj.step1, adj.step2 = input$adj.step2, adj.step3 = input$adj.step3)
      
      setProgress(message='Done !')
      
    })
    
  })
  #
  observeEvent(input$Test2,{
    
    withProgress({
      
      setProgress(message='Segmenting images...')
      
      ThumbIm$I = coloc.Sgl(MyImCl = ConInf$Welldat, Plate = input$SampPlate1, Time = input$SampTime1, Well= input$SampWell1, Site = input$SampSite1, Blue = as.numeric(input$BlueChannel), Green = as.numeric(input$GreenChannel), Red = as.numeric(input$RedChannel), auto2 = (input$auto2 == 'YES'), auto3 = (input$auto3 == 'YES'),
                            Cyto = input$Cyto, Nuc.rm = (input$Nucrm == 'YES'), TopSize2 = input$TopSize2, TopSize3 = input$TopSize3,  w1OFF = input$w1OFF,w2OFF = input$w2OFF,w3OFF = input$w3OFF, TEST=T, Rm1 = input$Rm1, Rm2 = input$Rm2, Rm3= input$Rm3,
                            W_Name = W_Name, getCell = (input$CellIm == 'YES'), adj.step1 = input$adj.step1, adj.step2 = input$adj.step2, adj.step3 = input$adj.step3)
      
      setProgress(message='Done !')
      
    })
    
  })
  #
  observeEvent(input$Test3,{
    
    withProgress({
      
      setProgress(message='Segmenting images...')
      
      ThumbIm$I = coloc.Sgl(MyImCl = ConInf$Welldat, Plate = input$SampPlate1, Time = input$SampTime1, Well= input$SampWell1, Site = input$SampSite1, Blue = as.numeric(input$BlueChannel), Green = as.numeric(input$GreenChannel), Red = as.numeric(input$RedChannel), auto2 = (input$auto2 == 'YES'), auto3 = (input$auto3 == 'YES'),
                            Cyto = input$Cyto, Nuc.rm = (input$Nucrm == 'YES'), TopSize2 = input$TopSize2, TopSize3 = input$TopSize3,  w1OFF = input$w1OFF,w2OFF = input$w2OFF,w3OFF = input$w3OFF, TEST=T, Rm1 = input$Rm1, Rm2 = input$Rm2, Rm3= input$Rm3,
                            W_Name = W_Name, getCell = (input$CellIm == 'YES'), adj.step1 = input$adj.step1, adj.step2 = input$adj.step2, adj.step3 = input$adj.step3)
      
      setProgress(message='Done !')
      
    })
    
  })
  
  #
  output$LookUp1 <-renderPlot({
    
    display(ThumbIm$I[[1]], method='raster')
    
  })
  output$LookUp2 <-renderPlot({
    
    display(ThumbIm$I[[2]], method='raster')
    
  })
  output$LookUp3 <-renderPlot({
    
    display(ThumbIm$I[[3]], method='raster')
    
  })
  output$LookUp4 <-renderPlot({
    
    display(ThumbIm$I[[4]], method='raster')
    
  })
  
  ## Give an idea of PCC
  
  output$SampPCC = renderText({
    round(ThumbIm$I[[5]],2)
  })
  
  #=====================================================================================================================================================================
  ## Launch Image Analysis
  
  observeEvent(input$launcher,{
    
    mytempdir = paste0(Plates$resultsloc,'/ID',paste(Plates$PlateIDs, collapse = '_'),'_',gsub(':','-',(gsub(' CEST','',Sys.time()))),'/')
    dir.create(mytempdir)
    
    #Parameters to be exported in foreach  
    MyImCl.FOR = ConInf$Welldat
    #--
    Plate = unique(MyImCl.FOR$PlateID)
    Time = unique(MyImCl.FOR$TimePoint)
    UniWell = unique(MyImCl.FOR$Well)
    UniSite = unique(MyImCl.FOR$Site)
    #--
    Blue.FOR = as.numeric(input$BlueChannel); Green.FOR = as.numeric(input$GreenChannel);Red.FOR = as.numeric(input$RedChannel); auto2.FOR = (input$auto2 == 'YES');auto3.FOR = (input$auto3 == 'YES');
    Cyto.FOR = input$Cyto;Nuc.rm.FOR = (input$Nucrm == 'YES'); TopSize2.FOR = input$TopSize2; TopSize3.FOR = input$TopSize3;w1OFF.FOR = input$w1OFF;w2OFF.FOR = input$w2OFF;w3OFF.FOR = input$w3OFF;
    getCell.FOR = (input$CellIm == 'YES'); as.FCS.FOR = (input$ExportFCS == 'YES');adj.step1.FOR = input$adj.step1; adj.step2.FOR = input$adj.step2; adj.step3.FOR = input$adj.step3;
    Rm1.FOR = input$Rm1; Rm2.FOR = input$Rm2; Rm3.FOR= input$Rm3; ExportResults.FOR = input$ExportResults ; ExpSeg.FOR = (input$ExpSeg == 'YES')
    #--
    
    #Log file
    log.file = cbind.data.frame(NbPlates = length(Plate), TopHatCyto1 = TopSize2.FOR,TopHatCyto2 = TopSize3.FOR, NucOffset = w1OFF.FOR, Cyto1Offset = w2OFF.FOR, Cyto2Offset = w3OFF.FOR, HistoAdjust.Nuc = adj.step1.FOR, 
                                HistoAdjust.Cyto1 = adj.step2.FOR, HistoAdjust.Cyto2 = adj.step3.FOR)
    
    withProgress({
      
      setProgress(message = "Optimizing Performance...")#########
      OS = .Platform$OS.type
      if(OS == 'windows'){
        RAM = shell('wmic OS get FreePhysicalMemory /Value',intern=T)
        RAM = RAM[grep('FreePhysicalMemory', RAM)]
        RAM = as.numeric(gsub('FreePhysicalMemory=','',RAM))
      }else{
        RAM = as.numeric(system(" awk '/MemFree/ {print $2}' /proc/meminfo", intern=T))
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
      #------------------------------------------------------------
      setProgress(message = "Creating clusters...")
      
      cl <- makeCluster(UsedCores)
      registerDoParallel(cl, cores = UsedCores)
      opts =list(chunkSize=1)
      #------------------------------------------------------------
      setProgress(message = "Treating images")
      t1=Sys.time()
      foreach(P = Plate, p=icount() ,.packages = c("EBImage","tiff","reshape","parallel","foreach","doParallel","iterators"),
              .combine = 'rbind',.options.nws=opts,.export='coloc.Sgl') %do% {
                
                P_Name = paste0(mytempdir,P)
                dir.create(P_Name)
                
                Summary<-foreach(TI = Time,t=icount() ,.packages = c("EBImage","tiff","reshape","parallel","foreach","doParallel","iterators"),
                                 .combine = 'rbind',.options.nws=opts,.export='coloc.Sgl') %do% {
                                   
                                   T_Name = paste0(P_Name,'/TimePoint_',TI)                    
                                   dir.create(T_Name)
                                   
                                   foreach(W = UniWell,i=icount() ,.packages = c("EBImage","tiff","reshape","parallel","foreach","doParallel","iterators","flowCore"),
                                           .combine = 'rbind',.options.nws=opts,.export='coloc.Sgl') %dopar% {
                                             
                                             W_Name = paste0(T_Name,'/',W)
                                             dir.create(W_Name)
                                             
                                             WSummary <- foreach(S = UniSite,j=icount(), .packages = c("EBImage","tiff","reshape"),.inorder=FALSE,
                                                                 .combine = 'rbind',.export='coloc.Sgl') %do% { 
                                                                   
                                                                   try({
                                                                     coloc.Sgl(MyImCl = MyImCl.FOR, Plate = P, Time = TI, Well= W, Site = S, Blue = Blue.FOR, Green = Green.FOR,Red = Red.FOR, auto2 = auto2.FOR, auto3 = auto3.FOR,
                                                                               Cyto = Cyto.FOR,Nuc.rm = Nuc.rm.FOR, TopSize2 = TopSize2.FOR, TopSize3 = TopSize3.FOR, w1OFF = w1OFF.FOR,w2OFF = w2OFF.FOR,w3OFF = w3OFF.FOR,
                                                                               adj.step1 = adj.step1.FOR, adj.step2 = adj.step2.FOR, adj.step3 = adj.step3.FOR, TEST=F,W_Name = W_Name,getCell = getCell.FOR, writeSeg = ExpSeg.FOR)
                                                                   })
                                                                   
                                                                 }
                                             
                                             if(getCell.FOR==T & as.FCS.FOR==T){
                                               FlowFrame = new("flowFrame",exprs=as.matrix(WSummary[-which(WSummary$ObjNum==0),c(12:17)]))
                                               write.FCS(FlowFrame, paste0(T_Name,'/T',TI,'_',W,'.fcs'))
                                             }
                                             print(WSummary)
                                             
                                           }
                                 }
                
                if(getCell.FOR==T & as.FCS.FOR==T){
                  FlowFrame = new("flowFrame",exprs=as.matrix(Summary[-which(Summary$ObjNum==0),c(12:17)]))
                  write.FCS(FlowFrame, paste0(P_Name,'/GlobalFCS.fcs'))
                }
                
                colnames(Summary)[c(1:4)]=c('PlateID', 'Time','WellID','SiteID')
                Summary = Summary[order(Summary$PlateID, Summary$Time,Summary$SiteID,Summary$WellID),,drop=F]
                Summary$GlobalID = paste(Summary$PlateID, Summary$Time, Summary$WellID, sep='_')
                
                ##
                pdf(paste0(P_Name,'/BoxPlot.pdf'),w=12,h=7)
                bp = boxplot(as.numeric(PCC) ~ GlobalID, data= Summary[!is.na(Summary$WellID),], ylim=c(-1,1), xaxt='n', outline=F, col='green',xpd=T,xaxs='i',yaxs='i',bty='n')
                stripchart(as.numeric(PCC) ~ GlobalID, data= Summary,vertical=T,pch=1,add=T,cex=0.3)
                text(c(1:length(bp$n)), rep(-1,length(bp$n)), gsub('_|NA',' ',bp$names),cex=0.7, srt=45,adj=c(1.1,1.1),xpd=T)
                dev.off()
                ##
                if(ExportResults.FOR == 'CSV'){
                  write.table(Summary[!is.na(Summary$WellID),c(1:11)],paste0(P_Name,'/Results.csv'),sep=',',row.names=F)
                }else{
                  ExpSummary = Summary[!is.na(Summary$WellID),c(1:11)]
                  save(ExpSummary, file = paste0(P_Name,'/Results.RData'))
                }
              }
      
      stopCluster(cl)
      setProgress(message = "End of treatments!")
      
      #---------------------------------------------------------
      
      setProgress(message = "Exporting data...")
      t2=Sys.time()
      write.csv(cbind.data.frame(log.file, Time = difftime(t2,t1, units = 'mins')),paste0(mytempdir,'log.csv'),sep=',',row.names=F)
      
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
