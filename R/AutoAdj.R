#' Returns automated inputRange
#'
#' @keywords internal
#' @import EBImage

#' @export

AutoAdjust = function(I = choose.files(), step=0.01){
  
  if(!is.character(I)){
    MAT = I
  }else{
    MAT = suppressWarnings(readImage(I))
  }
  
  qt = quantile(MAT, probs = seq(0,1,step))
  return(as.numeric(c(qt[2], qt[length(qt)-1])))
}