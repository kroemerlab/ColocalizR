#' Main function calculating colocalization parameters
#'
#' @keywords internal
#' @import flowCore
#' @import EBImage
#' @import tiff
#' @import gtools
#' @import MetaxpR
#' @import doParallel
#' @import MorphR

#' @export

coloc.Sgl = function(MyImCl, Plate,Time,Well,Site ,Blue=1,Green=2, Red=3, auto1=T,auto2=T,auto3=T, Cyto = 'Compt 2',Nuc.rm = T, Seg.method = 'Fast', TopSize1 = 29, TopSize2 = 29, TopSize3 = 29,  w1OFF = 0.1, w2OFF = 0.15, w3OFF = 0.1, Nuc.denoising = T, NucWindow = 50, RO.size = 25,
                     SegCyto.met = 'Automated', CytoWindow = 50, CytoOFF = 0.1, Rm1=c(0,0.3), Rm2=c(0,0.1),Rm3=c(0,0.1), adj = 1, adj.step1 = 2, adj.step2 = 2, adj.step3 = 2, getCell = F, add.features = F, writeSeg = T, writePDF = F, TEST = F, FullIm = F, getRange = c(F,F,F), path = '...'){  
  
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
      #--------------------------------------------------
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
      q = quantile(RAW,probs=seq(0,1,10**(-(adj.step[i]))))
      R = c(head(q,n=2)[2],tail(q,n=2)[1]);rm(q)
      #---
      LOG = MetaxpR::Normalize(log10(RAW+1))
      
      #-------------------------------------------------
      IMList = c(IMList,list(RAW))
      RAW = MetaxpR::Normalize(RAW,inputRange=R)
      #
      P = RAW-gblur(RAW,50)
      P[which(P<0)] = 0
      P = MetaxpR::Normalize(P)
      
      PList = c(PList,list(P));rm(P)
      
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
          CytoIm = (medianFilter(MetaxpR::Normalize((LOGG+LOGR)/2),10))
          rm(list = c('LOGR','LOGG'))
        }
      }
      if(length(grep('CytoIm',objects()))!=0){
        if(SegCyto.met == 'Automated'){
          CMask = CytoIm>(adj*otsu(CytoIm))
        }else if(SegCyto.met == 'Adaptive'){
          CMask = Reduce(function(x,y) {x|y}, lapply(seq(0.5,1.5,0.25)*as.numeric(CytoWindow), function(x) thresh(CytoIm, w = x, h = x, offset = CytoOFF)))
        }else if(SegCyto.met == 'Both'){
          CMask = Reduce(function(x,y) {x|y}, c(list(CytoIm>(adj*otsu(CytoIm))),lapply(seq(0.5,1.5,0.25)*as.numeric(CytoWindow), function(x) thresh(CytoIm, w = x, h = x, offset = CytoOFF))))
        }
      }
      
      if((Nuc.rm|getCell) & i==Blue){
        TH = thresh(RAW,w=NucWindow,h=NucWindow, offset = w1OFF)
        
        if(Nuc.denoising){
          TH = MorphR::MorphRecons(f = geoerode,reference = MorphR::MorphRecons(f = geodilate,reference = TH,element = makeBrush(RO.size,'disc'),oc='open'),
                                   element = makeBrush(RO.size,'disc'),oc='close')
        }
        T1 = opening(closing(fillHull(TH),makeBrush(5,'disc')),makeBrush(5,'disc'))

        if(Seg.method == 'Fast'){
          seed = erode(T1, makeBrush(15,'box'))
          NMask = propagate(T1, bwlabel(seed), mask = T1)
          rm(seed)
        }else if(Seg.method == 'Robust'){
          NMask = watershed(distmap(opening(TH,makeBrush(15,'disc'))))
        }
        rm(TH)
      }
      
      ##------------------------------------------------- 
      #Let's get uncleaned masks
      
      TOP = whiteTopHat(MorphR::LowPass(LOG),makeBrush(TopSize[i],'disc'))
      TOP[which(TOP<0)]=0
      TOP = MetaxpR::Normalize(TOP)
      
      if(autoSeg[i]){
        C = TOP>(otsu(TOP))
      }else{
        C = thresh(TOP,TopSize[i],TopSize[i],wOFF[i])
      }
      
      MList = c(MList,list(C))
      RAWList = c(RAWList,list(RAW))   
      
      rm(list=c('LOG','TOP','C','RAW'))
    }
  }
  gc()
  #-------------------------------------------------------
  
  if(nch>2){
    if(Nuc.rm){
      CMask = CMask & !T1
    }
    if(getCell){
      OMask = propagate(CMask,NMask,mask = CMask|NMask)
      rm(T1)
      
      if(all(OMask == 0)){
        getCell = F;Nuc.rm = F
      }
    }
  }
  #
  for(i in 1:nch){
    MList[[i]] = MList[[i]] & CMask
  }
  
  #COLOC
  COLOC =  MList[[Green]] & MList[[Red]]
  
  #UNION
  UNION =  MList[[Green]] | MList[[Red]]
  
  if(getCell){
    COLOC = COLOC*OMask
    UNION = UNION*OMask 
  }
  
  #================================================================================================================
  #IMAGE ANALYSIS
  
  #Converting images to arrays & calculate coeff
  pixG = c(PList[[Green]])
  pixR = c(PList[[Red]])
  
  ## Mask#
  pixM = c(UNION)
  
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
      
      testR = paintObjects(OMask,paintObjects(MList[[Red]],rgbImage(red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'green'),col='yellow')
      testG = paintObjects(OMask,paintObjects(MList[[Green]],rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2))), col = 'red'),col = 'yellow')
      
      if(Nuc.rm){
        testB = paintObjects(NMask,rgbImage(blue=(MetaxpR::Normalize(IMList[[Blue]],inputRange=Rm1))), col = 'yellow')
        
        if(Cyto == 'Compt 1'){
          testC = paintObjects(NMask,paintObjects(OMask,rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2))),col = 'yellow'),col='yellow')
        }else if(Cyto == 'Compt 2'){
          testC = paintObjects(NMask,paintObjects(OMask,rgbImage(red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'yellow'),col='yellow')
        }else if(Cyto == 'Both'){
          testC = paintObjects(NMask,paintObjects(OMask,rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2)), red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'yellow'),col='yellow')
        }
        testRGB = paintObjects(NMask,paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(OMask, rgbImage(blue = MetaxpR::Normalize(IMList[[Blue]],inputRange=Rm1), green=MetaxpR::Normalize(IMList[[Green]], inputRange=Rm2),red =MetaxpR::Normalize(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                                             col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), col = 'blue', thick = T)
      }else{
        testB = paintObjects(NMask,rgbImage(blue=(MetaxpR::Normalize(IMList[[Blue]],inputRange=Rm1))), col = 'yellow')
        
        if(Cyto == 'Compt 1'){
          testC = paintObjects(OMask,rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2))),col = 'yellow')
        }else if(Cyto == 'Compt 2'){
          testC = paintObjects(OMask,rgbImage(red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'yellow')
        }else if(Cyto == 'Both'){
          testC = paintObjects(OMask,rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2)), red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'yellow')
      }
    
        testRGB = paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(OMask, rgbImage(green=MetaxpR::Normalize(IMList[[Green]], inputRange=Rm2),red =MetaxpR::Normalize(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                            col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F)
      }
    }else{
      
      testR = paintObjects(CMask,paintObjects(MList[[Red]],rgbImage(red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'green'),col='yellow')
      testG = paintObjects(CMask,paintObjects(MList[[Green]],rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2))), col = 'red'),col = 'yellow')
      
      if(nch>2){
        if(Nuc.rm){
          testB = paintObjects(T1,rgbImage(blue=(MetaxpR::Normalize(IMList[[Blue]],inputRange=Rm1))), col = 'yellow')
          
          if(Cyto == 'Compt 1'){
            testC = paintObjects(T1,paintObjects(CMask,rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2))),col = 'yellow'),col='yellow')
          }else if(Cyto == 'Compt 2'){
            testC = paintObjects(T1,paintObjects(CMask,rgbImage(red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'yellow'),col='yellow')
          }else if(Cyto == 'Both'){
            testC = paintObjects(T1,paintObjects(CMask,rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2)), red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'yellow'),col='yellow')
          }
          testRGB = paintObjects(T1,paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(CMask, rgbImage(blue = MetaxpR::Normalize(IMList[[Blue]],inputRange=Rm1), green=MetaxpR::Normalize(IMList[[Green]], inputRange=Rm2),red =MetaxpR::Normalize(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                                              col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), col = 'blue', thick = T)
        }else{
          testB = rgbImage(blue=(MetaxpR::Normalize(IMList[[Blue]],inputRange=Rm1)))
          
          if(Cyto == 'Compt 1'){
            testC = paintObjects(CMask,rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2))),col = 'yellow')
          }else if(Cyto == 'Compt 2'){
            testC = paintObjects(CMask,rgbImage(red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'yellow')
          }else if(Cyto == 'Both'){
            testC = paintObjects(CMask,rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2)), red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'yellow')
          }
          testRGB = paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(CMask, rgbImage(green=MetaxpR::Normalize(IMList[[Green]], inputRange=Rm2),red =MetaxpR::Normalize(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                              col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F)
        }
      }else{
        testB = rgbImage(blue = matrix(0,nrow=dim(MList[[Green]])[1],ncol=dim(MList[[Green]])[2]))
        
        if(Cyto == 'Compt 1'){
          testC = paintObjects(CMask,rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2))),col = 'yellow')
        }else if(Cyto == 'Compt 2'){
          testC = paintObjects(CMask,rgbImage(red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'yellow')
        }else if(Cyto == 'Both'){
          testC = paintObjects(CMask,rgbImage(green=(MetaxpR::Normalize(IMList[[Green]],inputRange=Rm2)), red=(MetaxpR::Normalize(IMList[[Red]],inputRange=Rm3))),col = 'yellow')
        }
        testRGB = paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(CMask, rgbImage(green=MetaxpR::Normalize(IMList[[Green]], inputRange=Rm2),red =MetaxpR::Normalize(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                            col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F)
      }
    }
    return(list(CB=testB, CC=testC, CG=testG, CR=testR, CRGB=testRGB, RB=Rm1,RG=Rm2,RR=Rm3,PCC=PCC,SOC=SOC))
    
  }else{
    
    ICQ = ICQ.calc(cbind(pixGM,pixRM))
    MOC = MOC.calc(cbind(pixGM,pixRM)) #Manders Overlap coefficient
    SOCR = length(which(COLOC!=0))/length(which(MList[[Red]]==1))  #Surface Overlap coeff for red component
    SOCG = length(which(COLOC!=0))/length(which(MList[[Green]]==1))  #Surface Overlap coeff for green component
    
    #----------------------------------------------------------------------------------------------------------------
    #Correlation values calculation on a cell-by-cell basis
    
    if(getCell){
      
      PixTable = data.frame(Object = pixM[which(pixM!=0)], C2=pixGM, C3 = pixRM, COLOC = c(COLOC)[which(pixM!=0)])
      ObjNum = sort(unique(PixTable$Object))
      
      if(length(ObjNum) == 0){
        CellTable = data.frame(ObjNum = NA,PCC = NA, ICQ = NA, SOC = NA, SOCR = NA, SOCG = NA, MOC = NA)
      }else{
        CellTable = data.frame(ObjNum,do.call('rbind',by(PixTable[,c('C2','C3','COLOC')],PixTable$Object,function(x)data.frame(PCC=cor(x[,1],x[,2]),
                                                                                                                               ICQ = ICQ.calc(x[1:2]),
                                                                                                                               MOC = MOC.calc(x[1:2]),
                                                                                                                               SOC = length(which(x[,3]!=0))/length(x[,3]),
                                                                                                                               SOCR = length(which(x[,3]!=0))/length(which(x[,2]!=0)),
                                                                                                                               SOCG = length(which(x[,3]!=0))/length(which(x[,1]!=0))
        ))))
        
      }
      CellTable = data.frame(PlateID = Plate,Time = Time,WellID = Well,SiteID = Site, CellTable)
      rm(list=c('PixTable','ObjNum','COLOC','UNION','pixG','pixM','pixR'))
      gc()
      
      ## Features ##
      if(add.features){
        # Nucleus
        NucIDs = sort(unique(c(NMask[which(NMask!=0)])))
        Nuc_Data = data.frame(CellID = NucIDs,data.frame(computeFeatures(x = NMask,ref = IMList[[Blue]], methods.ref=,'computeFeatures.basic',methods.noref=c('computeFeatures.moment','computeFeatures.shape'),xname='Nucleus',refnames='Nuc')))
        
        # Cytoplasm
        CytoIDs = sort(unique(c(OMask[which(OMask!=0)])))
        Cyto_Data = data.frame(CellID = CytoIDs,data.frame(computeFeatures(x = OMask,ref= IMList[[Green]],methods.ref=,'computeFeatures.basic',methods.noref=c('computeFeatures.moment','computeFeatures.shape'), xname='Cyto',refnames='Cpt1')),
                               data.frame(computeFeatures(x = OMask,ref= IMList[[Red]],methods.ref=,'computeFeatures.basic',methods.noref=c('computeFeatures.moment','computeFeatures.shape'), xname='Cyto',refnames='Cpt2')))
        
        
        # Dots
        GDots_Mask = OMask*MList[[Green]]
        GDotsIDs = sort(unique(c(GDots_Mask[which(GDots_Mask!=0)])))
        GDots_Data = tapply(c(GDots_Mask[-which(GDots_Mask==0)]),c(GDots_Mask[-which(GDots_Mask==0)]),function(x)length(x))
        GDots_Data = data.frame(CellID = c(names(GDots_Data)), GDots.Surface = c(GDots_Data))
        
        RDots_Mask = OMask*MList[[Red]]
        RDotsIDs = sort(unique(c(RDots_Mask[which(RDots_Mask!=0)])))
        RDots_Data = tapply(c(RDots_Mask[-which(RDots_Mask==0)]),c(RDots_Mask[-which(RDots_Mask==0)]),function(x)length(x))
        RDots_Data = data.frame(CellID = c(names(RDots_Data)), RDots.Surface = c(RDots_Data))
        
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
            pdf(paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_PixelProfiling.pdf'),w=7,h=7)
            smoothScatter(pixGM, pixRM, nrpoints = 0, colramp=MyCols, main=paste(Well,'-',Site, sep=' '),nbin=512,
                          bandwidth=0.00005,xaxs='i',yaxs='i',xlab='Channel 2',ylab='Channel 3',useRaster=T,xlim=lim,ylim=lim)
            legend('topright',bty='n',cex=0.6,legend=c(paste('PCC=',round(PCC,2)),paste('ICQ=',round(ICQ,2)),paste('MOC=',round(MOC,2)),paste('SOC=',round(SOC,2))),text.col='red')
            dev.off()
          })
        }
      }else{
        if(!is.na(PCC)){
          try({
            pdf(paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_PixelProfiling.pdf'),w=7,h=7)
            smoothScatter(pixGM, pixRM, nrpoints = 0, colramp=MyCols, main=paste(Well,'-',Site, sep=' '),nbin=512,
                          bandwidth=0.00005,xaxs='i',yaxs='i',xlab='Channel 2',ylab='Channel 3',useRaster=T,xlim=lim,ylim=lim)
            legend('topright',bty='n',cex=0.6,legend=c(paste('PCC=',round(PCC,2)),paste('ICQ=',round(ICQ,2)),paste('MOC=',round(MOC,2)),paste('SOC=',round(SOC,2))),text.col='red')
            dev.off()
          })
        }
      }
    }
    rm(list=c('pixRM','pixGM'))
    ##
    if(Site==1 & writeSeg){
      if(getCell){
        if(nch>2){ 
          writeImage(paintObjects(NMask,paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(OMask, rgbImage(blue = MetaxpR::Normalize(IMList[[Blue]],inputRange=Rm1), green=MetaxpR::Normalize(IMList[[Green]], inputRange=Rm2),red =MetaxpR::Normalize(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                                                  col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), col = 'blue', thick = T), paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_ImSeg1.tif'),compression = 'LZW')
        }else{
          writeImage(paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(OMask, rgbImage(green=MetaxpR::Normalize(IMList[[Green]], inputRange=Rm2),red =MetaxpR::Normalize(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                               col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_ImSeg1.tif'),compression = 'LZW')
        }
      }else{
        if(nch>2 & Nuc.rm){
          writeImage(paintObjects(T1,paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(CMask, rgbImage(blue = MetaxpR::Normalize(IMList[[Blue]],inputRange=Rm1), green=MetaxpR::Normalize(IMList[[Green]], inputRange=Rm2),red =MetaxpR::Normalize(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                                               col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), col = 'blue', thick = T), paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_ImSeg1.tif'),compression = 'LZW')
        }else{
          writeImage(paintObjects(MList[[Green]], paintObjects(MList[[Red]],paintObjects(CMask, rgbImage(green=MetaxpR::Normalize(IMList[[Green]], inputRange=Rm2),red =MetaxpR::Normalize(IMList[[Red]], inputRange=Rm3)), col = 'yellow', thick = T),                                                                        
                                                               col = c('red', 'dark red'), opac = c(1,0.3), thick = F), col = c('green', 'darkgreen'), opac = c(1,0.3), thick = F), paste0(path, '/',Plate,'/',Time,'/',Well,'/',Well,'_s',Site,'_ImSeg1.tif'),compression = 'LZW')
        }
      }
    }
    rm(list=c(grep('Mask|List', ls(),value=T),'CytoIm'))
    gc()
    #Write results in arrays
    if(!getCell){
      return(data.frame(ObjNum=0,PlateID = Plate,Time = Time,WellID = Well,SiteID = Site,PCC = PCC,ICQ = ICQ,SOC = SOC,SOCR = SOCR,SOCG = SOCG,MOC = MOC))
    }else{
      return(smartbind(CellTable,data.frame(ObjNum=0,PlateID = Plate,Time = Time,WellID = Well,SiteID = Site,PCC = PCC,ICQ = ICQ,SOC = SOC,SOCR = SOCR,SOCG = SOCG,MOC = MOC)))  
    }
  }
}
