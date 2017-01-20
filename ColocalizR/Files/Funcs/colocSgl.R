library(flowCore)
library(EBImage)
library(tiff)

options(warn=-1)

coloc.Sgl = function(MyImCl,W_Name, Plate,Time,Well,Site ,Blue=1,Green=2, Red=3, auto1=T,auto2=T, auto3=T, Cyto = 'Compt 1',Nuc.rm = T, TopSize1 = 29, TopSize2 = 29, TopSize3 = 29,  w1OFF = 0.1, w2OFF = 0.15, w3OFF = 0.1, adj = 1,
                     Rm1=c(0,0.3), Rm2=c(0,0.1),Rm3=c(0,0.1), adj.step1 = 2, adj.step2 = 2, adj.step3 = 2, getCell = F,writeSeg = T, TEST = F, FullIm = F,Extract = 1000, getRange = c(F,F,F)){  
  
  #==================
  
  if(missing(MyImCl)){
    return('Need images information to be imported!')
  }
  
  if(missing(Well)){
    return('Need WellID to keep going!')
  }
  
  MyCols = colorRampPalette(c('white','black','darkblue','blue','green','red','orange','gold','yellow'))
  
  #==================
  
  #################################################################################################################  
  #IMAGE IMPORT
  Im = MyImCl[which(MyImCl$PlateID==Plate & MyImCl$TimePoint==Time & MyImCl$Well==Well & MyImCl$Site==Site),]
  #
  nch = dim(Im)[1]
  if(nch<2 | nch>3){
    return('Need at least two and at most three channels !')
  }
  if(nch<3 & getCell==T){
    return('Cannot do cell-by-cell analysis without nuclear stain !')
  }
  
  if(TEST==T && FullIm==F){
    Center = dim(readTIFF(as.character(Im$MyIm[1]), info=F))/2
    xEx = (Center[1]-Extract/2):(Center[1]+Extract/2)
    yEx = (Center[2]-Extract/2):(Center[2]+Extract/2)
  }
  
  IMList = list()  #Contains original images
  RAWList = list() #Contains adjusted images
  PList = list() #Contains measured images
  MList = list() #Mask list
  #
  wOFF = c(w1OFF, w2OFF, w3OFF)
  wOFF = wOFF[order(c(Blue, Green,Red))]
  #
  TopSize = c(TopSize1, TopSize2, TopSize3)
  TopSize = TopSize[order(c(Blue, Green,Red))]
  #
  autoSeg = c(auto1, auto2, auto3)
  autoSeg = autoSeg[order(c(Blue, Green,Red))]
  #
  adj.step = c(adj.step1, adj.step2, adj.step3)
  adj.step = adj.step[order(c(Blue, Green,Red))]
  
  if(unique(getRange)!=F){
    AutoRange = list()
  }
  ################################################################################################################
  
  for(i in 1:nch){
    IP = as.character(Im$MyIm[which(Im$Channel==paste0('w',i))])
    if(length(IP) == 1){
      RAW = readTIFF(IP, info=F)
      if(TEST==T & FullIm == F){
        RAW = RAW[xEx,yEx]
      }
      rm(IP)
      
      #############################################################################################################
      
      if(unique(getRange)!=F){
        AutoRange = c(AutoRange, list(AutoAdjust(RAW,step=0.001)))
      }
      #############################################################################################################
      #Gets Masks 
      
      LOG = RAW
      #---
      q = quantile(LOG,probs=seq(0,1,10^(-(adj.step[i]))))
      R = c(head(q,n=2)[2],tail(q,n=2)[1])
      #---
      LOG[which(LOG==0)] = min(LOG[which(LOG!=0)])
      LOG = EBImage::normalize(log10(LOG))
      
      #-------------------------------------------------
      IMList = c(IMList,list(RAW))
      RAW = EBImage::normalize(RAW,inputRange=R)
      #
      P = RAW-gblur(RAW,50)
      P[which(P<0)] = 0
      P = EBImage::normalize(P)
      
      PList = c(PList,list(P))
      
      ##-------------------------------------------------
      #Here we get LOGs and cytosol Mask
      
      if(Cyto =='Compt 1' & i == Green){
        CytoIm = (medianFilter(LOG,10))
      }else if(Cyto=='Compt 2' & i == Red){
        CytoIm = (medianFilter(LOG,10))
      }else if(Cyto == 'Both'){
        if(i == Green){
          LOGG = LOG
        }
        if(i == Red){
          LOGR = LOG
        }
        if(i == max(nch)){
          CytoIm = (medianFilter(EBImage::normalize((LOGG+LOGR)/2),10))
          rm(list = c('LOGR','LOGG'))
        }
      }
      if(length(grep('CytoIm',objects()))!=0){
        CMask = CytoIm>(otsu(CytoIm))
      }
      
      if((Nuc.rm==T|getCell==T) & i==Blue){
        T1 = opening(closing(fillHull(thresh(RAW,w=50,h=50, offset = w1OFF)),makeBrush(5,'disc')),makeBrush(5,'disc'))
        
        if(getCell==T){
          seed = erode(T1, makeBrush(19,'box'))
          NMask = propagate(T1, bwlabel(seed), mask = T1)
        }
      }
      
      ##-------------------------------------------------
      #Let's get uncleaned masks
      
      TOP = whiteTopHat(LOG,makeBrush(TopSize[i],'disc'))
      TOP[which(TOP<0)]=0
      TOP = EBImage::normalize(TOP)
      
      if(autoSeg[i]==T){
        C = TOP>(otsu(TOP)*adj)
      }else{
        C = thresh(TOP,TopSize[i],TopSize[i],wOFF[i])
      }
      
      MList = c(MList,list(C))
      RAWList = c(RAWList,list(RAW))   
      
      rm(list=c('LOG','TOP','C','RAW'))
      gc()
    }
  }
  
  #############################################################################################################
  if(nch>2){
    
    if(Nuc.rm==T){
      CMask = ceiling((CMask+T1)/2 - T1)
    }
    
    if(getCell==T){
      OMask = propagate(CMask,NMask,mask = ceiling((CMask+NMask)/2))
      rm(T1)
    }
  }
  #
  for(i in 1:nch){
    MList[[i]] = floor((MList[[i]]+CMask)/2)
    MList[[i]] = round(gblur(MList[[i]],0.6,3)) # Remove lonely pixels
  }
  
  #COLOC
  COLOC = floor((MList[[Green]]+MList[[Red]])/2)
  
  #UNION
  UNION = ceiling((MList[[Green]]+MList[[Red]])/2) ##OR operation
  
  if(getCell==T){
    COLOC = COLOC*OMask
    UNION = UNION*OMask 
  }
  
  ################################################################################################################
  #Converting images to arrays & calculate coeff
  
  pixG = as.numeric(PList[[Green]])
  pixR = as.numeric(PList[[Red]])
  
  ## Mask#
  pixM = as.numeric(UNION)
  
  #Grayscale In Mask
  pixGM = pixG[which(pixM!=0)]
  pixRM = pixR[which(pixM!=0)]
  
  ################################################################################################################
  
  lim=c(0,1)
  ####
  
  ICQ.calc = function(MAT){
    
    if(dim(MAT)[2]!=2){
      return('Must have two columns!')
    }
    Av1 = mean(MAT[,1],na.rm=T)
    Av2 = mean(MAT[,2], na.rm=T)
    
    return(((length(which(lapply(MAT[,1],function(x) (Av1-x))>=0 & (lapply(MAT[,2],function(x) (Av2-x)))>=0)) + length(which(lapply(MAT[,1],function(x) (Av1-x))<0 & (lapply(MAT[,2],function(x) (Av2-x)))<0))) / length(MAT[,1]))-0.5)
    
  }
  ###
  MOC.calc = function(MAT){
    
    if(dim(MAT)[2]!=2){
      return('Must have two columns!')
    }
    
    return((sum(MAT[,1]*MAT[,2]))/sqrt((sum(MAT[,1]^2))*(sum(MAT[,2]^2))))
    
  }
  ###
  
  ##############################################################################################################################
  #Correlation Values calculation on whole image
  
  PCC = cor(pixGM, pixRM)
  SOC = length(which(COLOC!=0))/(length(which(UNION!=0)))
  
  if(TEST==T){
    
    if(getRange[Blue] & nch>2){
      Rm1 = AutoRange[[Blue]]
    }
    if(getRange[Green]){
      Rm2 = AutoRange[[Green]]
    }
    if(getRange[Red]){
      Rm3 = AutoRange[[Red]]
    }
    
    if(getCell==T){
      
      if(nch>2){
        testB = paintObjects(NMask,rgbImage(blue=(EBImage::normalize(IMList[[Blue]],inputRange=Rm1))), col = 'yellow')
      }
      testG = paintObjects(OMask,paintObjects(MList[[Green]],rgbImage(green=(EBImage::normalize(IMList[[Green]],inputRange=Rm2))), col = 'red'),col = 'yellow')
      testR = paintObjects(OMask,paintObjects(MList[[Red]],rgbImage(red=(EBImage::normalize(IMList[[Red]],inputRange=Rm3))),col = 'green'),col='yellow')
    }else{
      if(nch>2 & Nuc.rm==T){
        testB = paintObjects(T1,rgbImage(blue=(EBImage::normalize(IMList[[Blue]],inputRange=Rm1))), col = 'yellow')
      }
      if(nch>2 & Nuc.rm==F){
        testB = rgbImage(blue=(EBImage::normalize(IMList[[Blue]],inputRange=Rm1)))
      }
      testG = paintObjects(CMask,paintObjects(MList[[Green]],rgbImage(green=(EBImage::normalize(IMList[[Green]],inputRange=Rm2))), col = 'red'),col = 'yellow')
      testR = paintObjects(CMask,paintObjects(MList[[Red]],rgbImage(red=(EBImage::normalize(IMList[[Red]],inputRange=Rm3))),col = 'green'),col='yellow')
    }
    
    if(nch>2){
      testRGB = paintObjects(MList[[Red]],paintObjects(MList[[Green]], rgbImage(blue = EBImage::normalize(IMList[[Blue]],inputRange=Rm1), green=EBImage::normalize(IMList[[Green]], inputRange=Rm2),red =EBImage::normalize(IMList[[Red]], inputRange=Rm3)),col=c('green',NA)),col=c('red',NA))  
    }else{
      testB = EBImage::as.Image(matrix(0,nrow=dim(MList[[Green]])[1],ncol=dim(MList[[Green]])[2]))
      testRGB = paintObjects(MList[[Red]],paintObjects(MList[[Green]], rgbImage(green=EBImage::normalize(IMList[[Green]], inputRange=Rm2),red =EBImage::normalize(IMList[[Red]], inputRange=Rm3)),col=c('green',NA)),col=c('red',NA))
    }
    
    return(list(testB, testG, testR, testRGB, Rm1, Rm2, Rm3, PCC, SOC))
    
  }else{
    
    ICQ = ICQ.calc(cbind(pixGM,pixRM))
    MOC = MOC.calc(cbind(pixGM,pixRM)) #Manders Overlap coefficient
    
    SOCR = length(which(COLOC!=0))/length(which(MList[[Red]]==1))  #Manders Overlap coeff for red component
    SOCG = length(which(COLOC!=0))/length(which(MList[[Green]]==1))  #Manders Overlap coeff for green component
    
    if(getCell==T){
      #Correlation calculation cell by cell
      
      PixTable = data.frame(Object = pixM[which(pixM!=0)], C2=pixGM, C3 = pixRM, COLOC = as.numeric(COLOC)[which(pixM!=0)])
      
      ObjNum = as.numeric(tapply(PixTable$Object, PixTable$Object, function(x)unique(x)))
      
      PCCi = as.numeric(unlist(by(PixTable[,c(2,3)],PixTable$Object,FUN=function(x)cor(x[,1],x[,2]))))
      ICQi = as.numeric(unlist(by(PixTable[,c(2,3)],PixTable$Object,FUN=function(x) ICQ.calc(x))))
      MOCi = as.numeric(unlist(by(PixTable[,c(2,3)],PixTable$Object,FUN=function(x) MOC.calc(x))))
      
      SOCi = as.numeric(tapply(PixTable$COLOC, PixTable$Object, function(x)length(x[which(x!=0)])/length(x)))
      SOCRi = as.numeric(unlist(by(PixTable[,c(3,4)],PixTable$Object,FUN=function(x)length(which(x[,2]!=0))/length(which(x[,1]!=0)))))
      SOCGi = as.numeric(unlist(by(PixTable[,c(2,4)],PixTable$Object,FUN=function(x)length(which(x[,2]!=0))/length(which(x[,1]!=0)))))
      
      CellTable = data.frame(Plate, Time, Well, Site, ObjNum, PCC = PCCi, ICQ = ICQi,SOC= SOCi, SOCR = SOCRi,SOCG = SOCGi,MOC = MOCi, PCC_FCS = 50*(PCCi + 1), ICQ_FCS = 100*(ICQi + 0.5), SOC_FCS = 100*SOCi, SOCR_FCS = 100*SOCRi, SOCG_FCS = 100*SOCGi, MOC_FCS = 100*MOCi)
      
    }
    
    
    ################################################################################################################################
    
    pdf(paste0(W_Name,'/',Well,'_s',Site,'_PixelProfiling.pdf'),w=10,h=10)
    smoothScatter(pixGM, pixRM, nrpoints = 0, colramp=MyCols, main=paste(Well,'-',Site, sep=' '),nbin=512,
                  bandwidth=0.00005,xaxs='i',yaxs='i',xlab='Channel 2',ylab='Channel 3',useRaster=T,xlim=lim,ylim=lim)
    legend('topright',bty='n',cex=0.6,legend=c(paste('PCC=',round(PCC,2)),paste('ICQ=',round(ICQ,2)),paste('MOC=',round(MOC,2)),paste('SOC=',round(SOC,2))),text.col='red')
    dev.off()
    
    if(Site==1 & writeSeg==T){
      
      if(nch==3){ 
        writeImage(paintObjects(MList[[Red]],paintObjects(MList[[Green]], rgbImage(blue = RAWList[[Blue]]*CMask,green=RAWList[[Green]]*CMask,red =RAWList[[Red]]*CMask),col=c('green',NA)),col=c('red',NA)),
                   paste0(W_Name,'/',Well,'_s',Site,'_ImSeg1.tif'))
      }else{
        writeImage(paintObjects(MList[[Red]],paintObjects(MList[[Green]], rgbImage(green=RAWList[[Green]]*CMask,red =RAWList[[Red]]*CMask),col=c('green',NA)),col=c('red',NA)),
                   paste0(W_Name,'/',Well,'_s',Site,'_ImSeg1.tif'))
      }
    }
    
    #Write results in arrays
    if(getCell==F){
      return(c(Plate,Time,Well,Site,PCC,ICQ,SOC,SOCR,SOCG,MOC, PCC_FCS = 50*(PCC + 1), ICQ_FCS = 100*(ICQ + 0.5), SOC_FCS = 100*SOC, SOCR_FCS = 100*SOCR, SOCG_FCS = 100*SOCG, MOC_FCS = 100*MOC))
    }else{
      return(rbind(CellTable,data.frame(Plate,Time,Well,Site,ObjNum=0,PCC,ICQ,SOC,SOCR,SOCG,MOC, PCC_FCS = 50*(PCC + 1), ICQ_FCS = 100*(ICQ + 0.5), SOC_FCS = 100*SOC, SOCR_FCS = 100*SOCR, SOCG_FCS = 100*SOCG, MOC_FCS = 100*MOC)))  
    }
  }
}
