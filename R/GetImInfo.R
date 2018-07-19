#' Get Information about images
#'
#' @keywords internal
#' @import MetaxpR
#' @import reshape2
#' @import gtools
#' @export

getImInfo = function(PlateID, SQL.use = F, SERVER = 'MDCStore',DB, PlateLoc='...', TimeCourse=F){
  
  if(SQL.use){
    if(missing(DB)){DB = GetMDCInfo(SERVER,Unix.diff=c('//H','/media/h'))}
    Plate = DB[which(DB$PlateID == PlateID),]
    PlateLoc = unique(Plate$PlateLoc)
  }
  
  #---------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(!TimeCourse){
    MyIm = list.files(PlateLoc,recursive=F)
  }else{
    MyIm = list.files(PlateLoc,recursive=T)
  }
  
  MyIm = MyIm[!grepl('thumb|Thumbs|HTD|Thumb',MyIm)]
  has.site = grepl("_s._w1",MyIm[1])
  #---------------------------------------------------------------------------------------------------------------------------------------------------
  
  MyImCl = do.call('smartbind', lapply(MyIm, function(x) {if(has.site){
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

  #---------------------------------------------------------------------------------------------------------------------------------------------------
  
  MyImCl$TimePoint = as.numeric(gsub('TimePoint_','',MyImCl$TimePoint))
  MyImCl$Site = as.numeric(gsub('s','',MyImCl$Site))
  MyImCl$Channel = gsub('(w[1-3]).*', '\\1', MyImCl$Channel)
  MyImCl = cbind(MyIm = paste(PlateLoc, MyIm, sep='/'),MyImCl)
  MyImCl$PlateID = PlateID
  MyImCl = MyImCl[order(MyImCl$PlateID, MyImCl$TimePoint,MyImCl$Well, MyImCl$Site),,drop=F]
  MyImCl$GlobalID = paste(MyImCl$PlateID, MyImCl$TimePoint, MyImCl$Well, MyImCl$Site, sep = '_')
  
  #---------------------------------------------------------------------------------------------------------------------------------------------------
  
  return(MyImCl)

}
