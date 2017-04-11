library(MiXR)
library(reshape2)
library(pbapply)
library(gtools)
library(naturalsort)

getImInfo = function(PlateID, SQL.use = F, SERVER = 'YODA-SERVER', PlateLoc='...', TimeCourse=F){
  
  if(SQL.use){
    Plate = GetMDCInfo(SERVER)[which(GetMDCInfo(SERVER)$PlateID == PlateID),]
    PlateLoc = unique(Plate$PlateLoc)
    
    OS = .Platform$OS.type
    if(OS == 'unix'){
      PlateLoc = paste0('/media/', gsub('//', '', PlateLoc))
      PlateLoc = gsub('Hcs-screen10', 'hcs-screen10', PlateLoc)
    }
  }
  #Get Images
  ###################################################################################################################################################
  if(TimeCourse==F){
    MyIm = list.files(PlateLoc,recursive=F)
  }else{
    MyIm = list.files(PlateLoc,recursive=T)
  }
  if(length(grep('thumb|Thumbs|HTD|Thumb',MyIm))!=0){
    MyIm = MyIm[-grep('thumb|Thumbs|HTD|Thumb',MyIm)]
  }
  
  has.site = grepl("_s._w1",MyIm[1])
  ########################################
  cl = makeCluster(detectCores())
  clusterEvalQ(cl, library("reshape2"))
  
  MyImCl = do.call('smartbind', parLapply(cl = cl, MyIm, function(x) {if(has.site){
    if(TimeCourse){
      data.frame(TimePoint = dirname(x), colsplit(basename(x), pattern = '_', names = c('Name', 'Well', 'Site', 'Channel')), stringsAsFactors = F)
    }else{
      data.frame(TimePoint = 1, colsplit(basename(x), pattern = '_', names = c('Name', 'Well', 'Site', 'Channel')), stringsAsFactors = F)
    }
  }else{
    if(TimeCourse){
      data.frame(TimePoint = dirname(x), colsplit(basename(x), pattern = '_', names = c('Name', 'Well', 'Channel')), Site = 1, stringsAsFactors = F)
    }else{
      data.frame(TimePoint = 1, colsplit(x, pattern='_',names=c('Name','Well','Channel')), Site = 1, stringsAsFactors = F)
    }
  }}))
  
  stopCluster(cl)
  
  #######################################
  
  MyImCl$TimePoint = as.numeric(gsub('TimePoint_','',MyImCl$TimePoint))
  MyImCl$Site = as.numeric(gsub('s','',MyImCl$Site))
  MyImCl$Channel = gsub('(w[1-3]).*', '\\1', MyImCl$Channel)
  MyImCl = cbind(MyIm = paste(PlateLoc, MyIm, sep='/'),MyImCl)
  MyImCl$PlateID = PlateID
  MyImCl = MyImCl[order(MyImCl$PlateID, MyImCl$TimePoint,MyImCl$Well, MyImCl$Site),,drop=F]
  
  MyImCl$GlobalID = paste(MyImCl$PlateID, MyImCl$TimePoint, MyImCl$Well, MyImCl$Site, sep = '_')
  
  ########################################
  
  return(MyImCl)
  
  ###################################################################################################################################################
  
}
