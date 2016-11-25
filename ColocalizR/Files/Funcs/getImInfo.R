library(EBImage)
library(tiff)
library(naturalsort)
library(reshape)
library(rChoiceDialogs)

library("parallel")
library("foreach")
library("doParallel")
library("iterators")

getImInfo = function(PlateID = NA, PlateLoc='...', TimeCourse=F){
  
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
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl, cores = detectCores()-1)
  
  opts =list(chunkSize=2)
  MyImCl<-foreach(IM = MyIm,i=icount() ,.packages = c("EBImage","tiff","reshape","parallel","foreach","doParallel","iterators"),
                  .combine = 'rbind',.options.nws=opts) %dopar% {
                    if(has.site){                 
                      if(TimeCourse ==T){
                        cbind(TimePoint=as.numeric(gsub('TimePoint_','',unlist(strsplit(IM,'/'))[1])),colsplit(unlist(strsplit(IM,'/'))[2],split='_',names=c('Name','Well','Site','Channel')))
                      }else{
                        cbind(TimePoint=1,colsplit(IM,split='_',names=c('Name','Well','Site','Channel')))
                      }
                    }else{
                      if(TimeCourse ==T){
                        cbind(TimePoint=as.numeric(gsub('TimePoint_','',unlist(strsplit(IM,'/'))[1])),colsplit(unlist(strsplit(IM,'/'))[2],split='_',names=c('Name','Well','Channel')),Site=1)
                      }else{
                        cbind(TimePoint=1,colsplit(IM,split='_',names=c('Name','Well','Channel')),Site=1)
                      }  
                    }
                  }
  stopCluster(cl)
  #######################################
  
  MyImCl$Site = as.numeric(gsub('s','',MyImCl$Site))
  MyImCl$Channel = gsub('.TIF|.tif|.tiff','',MyImCl$Channel)
  MyImCl$Channel = gsub('w1.*', 'w1', MyImCl$Channel)
  MyImCl$Channel = gsub('w2.*', 'w2', MyImCl$Channel)
  MyImCl$Channel = gsub('w3.*', 'w3', MyImCl$Channel)
  MyImCl = cbind(MyIm = paste(PlateLoc, MyIm, sep='/'),MyImCl)
  MyImCl$PlateID = PlateID
  MyImCl = MyImCl[order(MyImCl$PlateID, MyImCl$TimePoint,MyImCl$Well, MyImCl$Site),,drop=F]
  
  ########################################
  
  MyImCl$Well = as.character(MyImCl$Well)
  
  return(MyImCl)
  
  ###################################################################################################################################################
  
}
