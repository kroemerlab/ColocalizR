#' ICQ calculation
#'
#' @keywords internal
#' @export
ICQ.calc = function(MAT){
  if(ncol(MAT)!=2){
    return('Must have two columns!')
  }
  Av = apply(MAT,2,function(x)mean(x,na.rm=T))
  NM1 = sum(sapply(MAT[,1],function(x) (Av[1]-x))>=0 & (sapply(MAT[,2],function(x) (Av[2]-x)))>=0)
  NM2 = sum(sapply(MAT[,1],function(x) (Av[1]-x))<0 & (sapply(MAT[,2],function(x) (Av[2]-x)))<0)

  return(((NM1+NM2)/nrow(MAT))-0.5)
}


#' MOC calculation
#'
#' @keywords internal
#' @export
MOC.calc = function(MAT){
  
  if(ncol(MAT)!=2){
    return('Must have two columns!')
  }
  return((sum(MAT[,1]*MAT[,2]))/sqrt((sum(MAT[,1]**2))*(sum(MAT[,2]**2))))
}
