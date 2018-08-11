#' ICQ calculation
#'
#' @keywords internal
#' @export
ICQ.calc = function(MAT){
  
  if(dim(MAT)[2]!=2){
    return('Must have two columns!')
  }
  Av1 = mean(MAT[,1],na.rm=T)
  Av2 = mean(MAT[,2], na.rm=T)
  
  return(((length(which(lapply(MAT[,1],function(x) (Av1-x))>=0 & (lapply(MAT[,2],function(x) (Av2-x)))>=0)) + length(which(lapply(MAT[,1],function(x) (Av1-x))<0 & (lapply(MAT[,2],function(x) (Av2-x)))<0))) / length(MAT[,1]))-0.5)
  
}

#' MOC calculation
#'
#' @keywords internal
#' @export
MOC.calc = function(MAT){
  
  if(dim(MAT)[2]!=2){
    return('Must have two columns!')
  }
  
  return((sum(MAT[,1]*MAT[,2]))/sqrt((sum(MAT[,1]^2))*(sum(MAT[,2]^2))))
  
}
