#' Geomorphometric operations
#'
#' @keywords internal
#' @import EBImage
#' @export

#-------------------------------------------------
geodilate = function(marker, reference, element){
  return(pmin(dilate(marker, element), reference))
}
#-------------------------------------------------
geoerode = function(marker, reference, element){
  return(pmax(erode(marker, element), reference))
}
#-------------------------------------------------
ReconsDilation = function(marker, reference, element){    
  recons = marker
  result = reference
  
  while(any(recons != result)){
    result = recons
    recons = geodilate(result, reference, element)
  }
  return(result)
}
#-------------------------------------------------
ReconsErosion = function(marker, reference, element){    
  recons = marker
  result = reference
  
  while(any(recons != result)){
    result = recons
    recons = geoerode(result, reference, element)
  }
  return(result)
}
#-------------------------------------------------
ReconsOpening = function(image, element){
  
  marker = erode(image, element)
  
  recons = marker
  result = image
  
  while(any(recons != result)){
    result = recons
    recons = geodilate(result, image, element)
  }
  return(result)
}
#-------------------------------------------------
ReconsClosing = function(image, element){
  
  marker = dilate(image, element)
  
  recons = marker
  result = image
  
  while(any(recons != result)){
    result = recons
    recons = geoerode(result, image, element)
  }
  return(result)
}













