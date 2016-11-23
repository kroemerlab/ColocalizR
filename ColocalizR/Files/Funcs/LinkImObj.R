####################################################################################################################################################
library(EBImage)
library(tiff)

library("parallel")
library("foreach")
library("doParallel")
library("iterators")
#####################################################################################################################################################

LinkImObj = function(Obj, Lab, Feat){
  
  
  Link = tapply(as.numeric(Obj[which(Lab!=0)]), as.numeric(Lab[which(Lab!=0)]), function(x) unique(x))
  FeatList = list()
  
  for(i in 1:length(Link)){
    
    FeatList = c(FeatList, list(Feat[(1:(dim(Feat)[1])) %in% Link[[i]],]))
    
  }
  
  names(FeatList) = sort(unique(as.numeric(Lab[which(Lab!=0)])))
  return(FeatList)

  
}