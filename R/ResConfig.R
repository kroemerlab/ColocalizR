#' Calculate ressources allocated to app
#'
#' @keywords internal
#' @import doParallel
#' @export

ResConfig = function(CellSeg = T){
  
  OS = Sys.info()[['sysname']]
  if(OS == 'Windows'){
    RAM = shell('wmic OS get FreePhysicalMemory /Value',intern=T)
    RAM = RAM[grep('FreePhysicalMemory', RAM)]
    RAM = as.numeric(gsub('FreePhysicalMemory=','',RAM))
  }else if(OS == 'Darwin'){
    RAM = as.numeric(unlist(strsplit(system("sysctl hw.memsize", intern = T), split = ' '))[2])/1e3
  }else{
    RAM = as.numeric(system(" awk '/MemFree/ {print $2}' /proc/meminfo", intern=T))
  }
  Cores = detectCores()
  if(CellSeg){
    Core2RAM = 2e06 # Assuming one core uses 2GB of RAM
  }else{
    Core2RAM = 1e06  # Assuming one core uses 1GB of RAM
  }
  MaxCores = floor(RAM/Core2RAM) 
  
  if(MaxCores>=Cores){
    UsedCores = Cores
  }else{
    UsedCores = MaxCores
  }
  
  return(UsedCores)
  
}