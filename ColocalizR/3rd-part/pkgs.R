#Add libPath,not requiring superuser privileges#
Rdir = paste(shell("echo %USERPROFILE%", intern=T),'\\Documents\\R',sep='')
if(!dir.exists(Rdir)){dir.create(Rdir)}

libdir = paste(Rdir,'\\win-library',sep='')
if(!dir.exists(libdir)){dir.create(libdir)}
libdir = paste(libdir,'\\3.3',sep='')
if(!dir.exists(libdir)){dir.create(libdir)}
#
.libPaths(new = libdir)
shell(paste('setx R_LIBS_USER ',libdir,sep=''))
#
source("http://bioconductor.org/biocLite.R")
biocLite(pkgs=c('EBImage','flowCore'), ask=F)
install.packages(pkgs = c('shiny','tiff','reshape','RODBC','foreach','doParallel','stringi','naturalsort','rChoiceDialogs'),
                 lib = libdir, repos = "http://cloud.r-project.org")
