MX.getAllSQLInfo = function(SERVER='YODA-SERVER'){
 
  options(warn=-1)
  library(RODBC)
  library(reshape)
  library(naturalsort)
  library(stringi)

  #Query DB ###########################################################################################################
  
  ch <- odbcConnect(SERVER, uid = "sa", pwd = "moldev") #Open DataBase
  
  #Data 
  Cor1 <- sqlFetch(ch,'ASSAY_PLATES')  ##Correspondence between MIC and Plate ID
  Cor2 <- sqlFetch(ch,'PLATES') ##All plates Information
  Cor3 <- sqlQuery(ch, "SELECT ASSAY_ID,ASSAY_NAME, SETTINGS_NAME, TABLE_ID FROM ASSAYS")#Custom Modules information
  Cor4 <- sqlFetch(ch,'PLATE_PROPERTY')##Get experiment set
  
  #Images location
  FileLoc <- sqlFetch(ch, 'FILE_LOCATION')
  
  odbcClose(ch)
  
  #Add plate ID to image location table
  FileLoc = FileLoc[!is.na(FileLoc$DIRECTORY),]
  FileLoc$ImPath = gsub('[\\]','/',FileLoc$DIRECTORY)
  FileLoc$ImPath = gsub('/TimePoint_.*','',FileLoc$ImPath)
    
  Plates = FileLoc$ImPath
  Plates = gsub('/TimePoint_.*','',Plates)
  
  #
  PlateIDs <- as.numeric(sapply(strsplit(Plates,"/"), tail, 1))
  #

  FileLoc$Plate.ID = PlateIDs
  FileLoc = FileLoc[!duplicated(FileLoc$Plate.ID),]
  #
  FileLoc$ImFullPath = paste0(gsub('[\\]','/',FileLoc$SERVER_NAME),FileLoc$ImPath)
  
  #####################################################################################################################
  #Create table !
  
  mixCor = data.frame(Cor1,Cor3) #Now every "Cor" has a plate ID reference
  Cor2 = Cor2[order(Cor2$PLATE_ID),,drop=F]
  Cor4 = Cor4[order(Cor4$PLATE_ID),,drop=F]
  FileLoc = FileLoc[order(FileLoc$Plate.ID),,drop=F]
  #
  mixCor = mixCor[order(mixCor$PLATE_ID),,drop=F]
  
  #Plates with measurements (easy !)
  Info = merge(mixCor, Cor2, by.x='PLATE_ID',by.y='PLATE_ID', all=T)
  Info = merge(Info, Cor4, by.x = 'PLATE_ID', by.y = 'PLATE_ID',all=T)
  Info = merge(Info, FileLoc, by.x = 'PLATE_ID', by.y = 'Plate.ID',all=T)
  
  ##
  Info = Info[,colnames(Info)[-c(2,5,8,9,10,11,12,18,20,21,22,23,24,26,27,28,29,30,31,32)]]
  colnames(Info) = c('PlateID','MeasurementID','BeingDeleted','AssayName','SettingsName','BarCode','Identifier',
                     'PlateName','PlateDesc','TimeCreated','TimeModified','ExpSet','PlateLoc')

  Info = Info[order(Info$PlateID, Info$MeasurementID),,drop=F]
  
  ###########################################################################################################################
  return(Info)
  
}