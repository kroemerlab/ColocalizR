#' Main function calculating colocalization parameters
#'
#' @keywords internal
#' @import flowCore
#' @import EBImage
#' @import tiff
#' @import gtools
#' @import MiXR
#' @import doParallel

#' @export

coloc.Sgl = function(MyImCl, Plate,Time,Well,Site ,Blue=1,Green=2, Red=3, auto1=T,auto2=T, auto3=T, Cyto = 'Compt 2',Nuc.rm = T, Seg.method = 'Fast', TopSize1 = 29, TopSize2 = 29, TopSize3 = 29,  w1OFF = 0.1, w2OFF = 0.15, w3OFF = 0.1, Nuc.denoising = T,
                     RO.size = 25, Rm1=c(0,0.3), Rm2=c(0,0.1),Rm3=c(0,0.1), adj = 1, adj.step1 = 2, adj.step2 = 2, adj.step3 = 2, getCell = F, add.features = F, writeSeg = T, writePDF = F, TEST = F, FullIm = F, getRange = c(F,F,F), path = '...'){  
  
  #================================================================================================================
  
  if(missing(MyImCl)){
    return('Need images information to be imported!')
  }
  
  if(missing(Well)){
    return('Need WellID to keep going!')
  }
  
  MyCols = colorRampPalette(c('white','black','darkblue','blue','green','red','orange','gold','yellow'))
  
  #================================================================================================================
  #IMAGE IMPORT
  Im = MyImCl[which(MyImCl$PlateID==Plate & MyImCl$TimePoint==Time & MyImCl$Well==Well & MyImCl$Site==Site),]
  #
  nch = nrow(Im)
  if(nch<2 | nch>3){
    return('Need at least two and at most three channels !')
  }
  if(nch<3 & getCell){
    return('Cannot do cell-by-cell analysis without nuclear stain !')
  }
  if(getCell){
    noNucMask = F
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
  
  if(any(getRange)){
    AutoRange = list()
  }
  #================================================================================================================
  #IMAGE TREATMENT
  
  for(i in 1:nch){
    IP = as.character(Im$MyIm[which(Im$Channel==paste0('w',i))])
    if(length(IP) == 1){
      RAW = suppressWarnings({readImage(IP)});rm(IP)
      #--------------------------------------------------------------------------------------------------------------
      
      if(any(getRange)){
        if(getRange[i]){
          AutoRange = c(AutoRange, list(AutoAdjust(RAW,step=0.001)))
        }else{
          AutoRange = c(AutoRange, NA)
        }
      }
      #---------------------------------------------------
      #Gets Masks 
      
      #---
      q = quantile(RAW,probs=seq(0,1,10^(-(adj.step[i]))))
      R = c(head(q,n=2)[2],tail(q,n=2)[1]);rm(q);gc()
      #---
      LOG = NormalizeIm(log10(RAW+1))

      
      #-------------------------------------------------
      IMList = c(IMList,list(RAW))
      RAW = NormalizeIm(RAW,inputRange=R)
      #
      P = RAW-gblur(RAW,50)
      P[which(P<0)] = 0
      P = NormalizeIm(P)
      
      PList = c(PList,list(P));rm(P);gc()
      
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
        if(i == 3){
          CytoIm = (medianFilter(NormalizeIm((LOGG+LOGR)/2),10))
          rm(list = c('LOGR','LOGG'))
        }
      }
      if(length(grep('CytoIm',objects()))!=0){
        CMask = CytoIm>(adj*otsu(CytoIm))
      }
      
      if((Nuc.rm|getCell) & i==Blue){
        if(getCell){
          Nuc = RAW
          if(Nuc.denoising){
            Nuc = ReconsOpening(Nuc, makeBrush(RO.size,'disc'))
          }
          T1 = opening(closing(fillHull(thresh(Nuc,w=50,h=50, offset = w1OFF)),makeBrush(5,'disc')),makeBrush(5,'disc'))
          
          if(Seg.method == 'Fast'){
            seed = erode(T1, makeBrush(19,'box'))
            NMask = propagate(T1, bwlabel(seed), mask = T1)
          }else if(Seg.method == 'Robust'){
            NMask = watershed(distmap(opening(thresh(Nuc,w=50,h=50, offset = w1OFF),makeBrush(15,'disc'))))
          }
          rm(list=c('Nuc','seed'));gc()
        }
      }
      
      ##-------------------------------------------------
      #Let's get uncleaned masks
      
      TOP = whiteTopHat(LOG,makeBrush(TopSize[i],'disc'))
      TOP[which(TOP<0)]=0
      TOP = NormalizeIm(TOP)


      if(autoSeg[i]){
        C = TOP>(otsu(TOP))
      }else{
        C = thresh(TOP,TopSize[i],TopSize[i],wOFF[i])
      }
      
      MList = c(MList,list(C))
      RAWList = c(RAWList,list(RAW))   
      
      rm(list=c('LOG','TOP','C','RAW'))
      gc()
    }
  }
  
  #-------------------------------------------------------
  
  if(nch>2){
    if(Nuc.rm){
      CMask = ceiling((CMask+T1)/2 - T1)
    }
    
    if(getCell){
      OMask = propagate(CMask,NMask,mask = ceiling((CMask+NMask)/2))
      rm(T1)
      
      if(all(OMask == 0)){
        getCell = F ; Nuc.rm = F ; noNucMask = T
      }
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
  
  if(getCell){
    COLOC = COLOC*OMask
    UNION = UNION*OMask 
  }
  
  #================================================================================================================
  #IMAGE ANALYSIS
  
  #Converting images to arrays & calculate coeff
  pixG = as.numeric(PList[[Green]])
  pixR = as.numeric(PList[[Red]])
  
  ## Mask#
  pixM = as.numeric(UNION)
  
  #Grayscale In Mask
  pixGM = pixG[which(pixM!=0)]
  pixRM = pixR[which(pixM!=0)]
  
  #----------------------------------------------------------------------------------------------------------------
  #Correlation Values, calculation on whole image
  lim=c(0,1)
  PCC = cor(pixGM, pixRM)
  SOC = length(which(COLOC!=0))/(length(which(UNION!=0)))
  
  if(getRange[Blue] & nch>2){
    Rm1 = AutoRange[[Blue]]
  }
  if(getRange[Green]){
    Rm2 = AutoRange[[Green]]
  }
  if(getRange[Red]){
    Rm3 = AutoRange[[Red]]
  }
  
  if(TEST){
    
    if(getCell){
      
      testR = paintObjects(OMask,paintObjects(MList[[Red]],rgbImage(red=(NormalizeIm(IMList[[Red]],inputRange=Rm3))),col = 'green'),col='yellow')
      testG = paintObjects(OMask,paintObjects(MList[[Green]],rgbImage(green=(NormalizeIm(IMList[[Green]],inputRange=Rm2))), col = 'red'),col = 'yellow')
      
      
      if(Nuc.rm){
        testB = paintObjects(NMask,rgbImage(blue=(NormalizeIm(IMList[[Blue]],inputRange=Rm1))), col = 'yellow')
        testRGB = paintObjects(NMask,paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(OMask, rgbImage(blue = NormalizeIm(IMList[[Blue]],inputRange=Rm1), green=NormalizeIm(IMList[[Green]], inputRange=Rm2),red =NormalizeIm(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                                               col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), col = 'blue', thick = T)
      }else{
        testB = rgbImage(blue=(NormalizeIm(IMList[[Blue]],inputRange=Rm1)))
        testRGB = paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(OMask, rgbImage(green=NormalizeIm(IMList[[Green]], inputRange=Rm2),red =NormalizeIm(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                            col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F)
      }
    }else{
      
      testR = paintObjects(CMask,paintObjects(MList[[Red]],rgbImage(red=(NormalizeIm(IMList[[Red]],inputRange=Rm3))),col = 'green'),col='yellow')
      testG = paintObjects(CMask,paintObjects(MList[[Green]],rgbImage(green=(NormalizeIm(IMList[[Green]],inputRange=Rm2))), col = 'red'),col = 'yellow')
      
      if(nch>2){
        if(Nuc.rm){
          testB = paintObjects(T1,rgbImage(blue=(NormalizeIm(IMList[[Blue]],inputRange=Rm1))), col = 'yellow')
          testRGB = paintObjects(T1,paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(CMask, rgbImage(blue = NormalizeIm(IMList[[Blue]],inputRange=Rm1), green=NormalizeIm(IMList[[Green]], inputRange=Rm2),red =NormalizeIm(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                                              col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), col = 'blue', thick = T)
        }else{
          testB = rgbImage(blue=(NormalizeIm(IMList[[Blue]],inputRange=Rm1)))
          testRGB = paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(CMask, rgbImage(green=NormalizeIm(IMList[[Green]], inputRange=Rm2),red =NormalizeIm(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                              col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F)
        }
      }else{
        testB = rgbImage(blue = matrix(0,nrow=dim(MList[[Green]])[1],ncol=dim(MList[[Green]])[2]))
        testRGB = paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(CMask, rgbImage(green=NormalizeIm(IMList[[Green]], inputRange=Rm2),red =NormalizeIm(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                            col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F)
      }
    }
    return(list(testB, testG, testR, testRGB, Rm1, Rm2, Rm3, PCC, SOC))
    
  }else{
    
    ICQ = ICQ.calc(cbind(pixGM,pixRM))
    MOC = MOC.calc(cbind(pixGM,pixRM)) #Manders Overlap coefficient
    SOCR = length(which(COLOC!=0))/length(which(MList[[Red]]==1))  #Surface Overlap coeff for red component
    SOCG = length(which(COLOC!=0))/length(which(MList[[Green]]==1))  #Surface Overlap coeff for green component
    
    #----------------------------------------------------------------------------------------------------------------
    #Correlation values calculation on a cell-by-cell basis
    
    if(getCell){
      
      PixTable = data.frame(Object = pixM[which(pixM!=0)], C2=pixGM, C3 = pixRM, COLOC = as.numeric(COLOC)[which(pixM!=0)])
      ObjNum = as.numeric(tapply(PixTable$Object, PixTable$Object, function(x)unique(x)))
      
      if(length(ObjNum) == 0){
        ObjNum = NA ; PCCi = NA ; ICQi = NA ; SOCi = NA ; SOCRi = NA ; SOCGi = NA ; MOCi = NA
      }else{
        PCCi = as.numeric(unlist(by(PixTable[,c(2,3)],PixTable$Object,FUN=function(x)cor(x[,1],x[,2]))))
        ICQi = as.numeric(unlist(by(PixTable[,c(2,3)],PixTable$Object,FUN=function(x) ICQ.calc(x))))
        MOCi = as.numeric(unlist(by(PixTable[,c(2,3)],PixTable$Object,FUN=function(x) MOC.calc(x))))
        
        SOCi = as.numeric(tapply(PixTable$COLOC, PixTable$Object, function(x)length(x[which(x!=0)])/length(x)))
        SOCRi = as.numeric(unlist(by(PixTable[,c(3,4)],PixTable$Object,FUN=function(x)length(which(x[,2]!=0))/length(which(x[,1]!=0)))))
        SOCGi = as.numeric(unlist(by(PixTable[,c(2,4)],PixTable$Object,FUN=function(x)length(which(x[,2]!=0))/length(which(x[,1]!=0)))))
      }
      CellTable = data.frame(ObjNum, PlateID = Plate,Time = Time,WellID = Well,SiteID = Site,PCC = PCCi, ICQ = ICQi,SOC= SOCi, SOCR = SOCRi,SOCG = SOCGi,MOC = MOCi,
                             PCC_FCS = 50*(PCCi + 1), ICQ_FCS = 100*(ICQi + 0.5), SOC_FCS = 100*SOCi, SOCR_FCS = 100*SOCRi, SOCG_FCS = 100*SOCGi, MOC_FCS = 100*MOCi)
      
      rm(list=c('PCCi','ICQi','MOCi','SOCi','SOCRi','SOCGi','PixTable','ObjNum','COLOC','UNION','pixG','pixM','pixR','pixRM','pixGM'))
      gc()
      
      ## Features ##
      if(add.features){
        # Nucleus
        NucIDs = sort(unique(as.numeric(NMask[which(NMask!=0)])))
        Nuc_Data = data.frame(CellID = NucIDs,data.frame(computeFeatures(x = NMask,ref = IMList[[Blue]], methods.ref=,'computeFeatures.basic',methods.noref=c('computeFeatures.moment','computeFeatures.shape'),xname='Nucleus',refnames='Nuc')))
        
        # Cytoplasm
        CytoIDs = sort(unique(as.numeric(OMask[which(OMask!=0)])))
        Cyto_Data = data.frame(CellID = CytoIDs,data.frame(computeFeatures(x = OMask,ref= IMList[[Green]],methods.ref=,'computeFeatures.basic',methods.noref=c('computeFeatures.moment','computeFeatures.shape'), xname='Cyto',refnames='Cpt1')),
                               data.frame(computeFeatures(x = OMask,ref= IMList[[Red]],methods.ref=,'computeFeatures.basic',methods.noref=c('computeFeatures.moment','computeFeatures.shape'), xname='Cyto',refnames='Cpt2')))

        
        # Dots
        GDots_Mask = OMask*MList[[Green]]
        GDotsIDs = sort(unique(as.numeric(GDots_Mask[which(GDots_Mask!=0)])))
        GDots_Data = tapply(as.numeric(GDots_Mask[-which(GDots_Mask==0)]),as.numeric(GDots_Mask[-which(GDots_Mask==0)]),function(x)length(x))
        GDots_Data = data.frame(CellID = as.numeric(names(GDots_Data)), GDots.Surface = as.numeric(GDots_Data))
        
        RDots_Mask = OMask*MList[[Red]]
        RDotsIDs = sort(unique(as.numeric(RDots_Mask[which(RDots_Mask!=0)])))
        RDots_Data = tapply(as.numeric(RDots_Mask[-which(RDots_Mask==0)]),as.numeric(RDots_Mask[-which(RDots_Mask==0)]),function(x)length(x))
        RDots_Data = data.frame(CellID = as.numeric(names(RDots_Data)), RDots.Surface = as.numeric(RDots_Data))
        
        CellFeatures = Reduce(function(x, y) merge(x, y,by = 'CellID', all=TRUE), list(Nuc_Data, Cyto_Data, GDots_Data, RDots_Data))
        CellTable = merge(CellTable, CellFeatures, by.x ='ObjNum', by.y = 'CellID', all.x = T, all.y = F)
        rm(list=grep(paste0(c('IDs','Data','Dots_Mask','CellFeatures'), collapse = '|'),ls()), value = T)
        gc()
      }
    }
    
    
    #================================================================================================================
    #DATA EXPORT
    if(writePDF){
      if(getCell){
        if(!is.na(CellTable$ObjNum)){
          try({
            pdf(paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_PixelProfiling.pdf'),w=10,h=10)
            smoothScatter(pixGM, pixRM, nrpoints = 0, colramp=MyCols, main=paste(Well,'-',Site, sep=' '),nbin=512,
                          bandwidth=0.00005,xaxs='i',yaxs='i',xlab='Channel 2',ylab='Channel 3',useRaster=T,xlim=lim,ylim=lim)
            legend('topright',bty='n',cex=0.6,legend=c(paste('PCC=',round(PCC,2)),paste('ICQ=',round(ICQ,2)),paste('MOC=',round(MOC,2)),paste('SOC=',round(SOC,2))),text.col='red')
            dev.off()
          })
        }
      }else{
        if(!is.na(PCC)){
          try({
            pdf(paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_PixelProfiling.pdf'),w=10,h=10)
            smoothScatter(pixGM, pixRM, nrpoints = 0, colramp=MyCols, main=paste(Well,'-',Site, sep=' '),nbin=512,
                          bandwidth=0.00005,xaxs='i',yaxs='i',xlab='Channel 2',ylab='Channel 3',useRaster=T,xlim=lim,ylim=lim)
            legend('topright',bty='n',cex=0.6,legend=c(paste('PCC=',round(PCC,2)),paste('ICQ=',round(ICQ,2)),paste('MOC=',round(MOC,2)),paste('SOC=',round(SOC,2))),text.col='red')
            dev.off()
          })
        }
      }
    }
    
    
    if(Site==1 & writeSeg){
      if(getCell){
        if(nch>2){ 
          writeImage(paintObjects(NMask,paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(OMask, rgbImage(blue = NormalizeIm(IMList[[Blue]],inputRange=Rm1), green=NormalizeIm(IMList[[Green]], inputRange=Rm2),red =NormalizeIm(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                                                  col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), col = 'blue', thick = T), paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_ImSeg1.tif'))
        }else{
          writeImage(paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(OMask, rgbImage(green=NormalizeIm(IMList[[Green]], inputRange=Rm2),red =NormalizeIm(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                               col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_ImSeg1.tif'))
        }
      }else{
        if(nch>2  & !noNucMask){ 
          writeImage(paintObjects(T1,paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(CMask, rgbImage(blue = NormalizeIm(IMList[[Blue]],inputRange=Rm1), green=NormalizeIm(IMList[[Green]], inputRange=Rm2),red =NormalizeIm(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                                               col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), col = 'blue', thick = T), paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_ImSeg1.tif'))
        }else{
          writeImage(paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(CMask, rgbImage(green=NormalizeIm(IMList[[Green]], inputRange=Rm2),red =NormalizeIm(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                               col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_ImSeg1.tif'))
        }
      }
    }
    
    rm(list=c(grep('Mask|List', ls(),value=T),'CytoIm'))
    gc()
    
    #Write results in arrays
    if(!getCell){
      return(data.frame(ObjNum=0,PlateID = Plate,Time = Time,WellID = Well,SiteID = Site,PCC = PCC,ICQ = ICQ,SOC = SOC,SOCR = SOCR,SOCG = SOCG,MOC = MOC,
               PCC_FCS = 50*(PCC + 1), ICQ_FCS = 100*(ICQ + 0.5), SOC_FCS = 100*SOC, SOCR_FCS = 100*SOCR, SOCG_FCS = 100*SOCG, MOC_FCS = 100*MOC))
    }else{
      return(smartbind(CellTable,data.frame(ObjNum=0,PlateID = Plate,Time = Time,WellID = Well,SiteID = Site,PCC = PCC,ICQ = ICQ,SOC = SOC,SOCR = SOCR,SOCG = SOCG,MOC = MOC,
                                            PCC_FCS = 50*(PCC + 1), ICQ_FCS = 100*(ICQ + 0.5), SOC_FCS = 100*SOC, SOCR_FCS = 100*SOCR, SOCG_FCS = 100*SOCG, MOC_FCS = 100*MOC)))  
    }
  }
}
