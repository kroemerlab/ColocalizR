# ColocalizR
ColocalizR is an R-based image-analysis application developed for the quantification of co-localization in high-throughput, on a cell-by-cell basis.

## Requirements
ColocalizR can be used on Windows, Linux and MacOS X, with at least 8GB of RAM and a quad-core (or higher) last-generation CPU. Internet connection is required for installation.

## Installation
To use this app, you need to have [R](https://cran.r-project.org/), [RStudio](https://www.rstudio.com/products/rstudio/download/) and [Java](https://www.java.com/fr/) installed.

### Windows
Click on the previous links to download and install all the required program files. 

### MacOS X
Install the program files in the same way than in Windows configuration. On Mac, Java is not required but you have to install Xcode either by downloading it via the AppStore or by installing it directly via the terminal.
The easiest way on on MacOS 10.9 or later is to type in a terminal:
```
xcode-select --install
```
It will be asked you to install Command Line Tools, you just have to click on *Install* in the pop-up window.

For previous versions of MacOS, you can refer to [these instructions](https://www.moncefbelyamani.com/how-to-install-xcode-homebrew-git-rvm-ruby-on-mac/) to install basic programmation tools on Mac (including Xcode).

### Linux
On Linux, you can either install the program files in the same way than in the two previous configurations or use command lines to do it. If you choose to do it in command lines, type in a terminal the commands below to correctly install and configure R, RStudio and Java. 

#### Install R
##### Ubuntu
```sh
grep -q -F "deb http://cran.rstudio.com/bin/linux/ubuntu *UbuntuVersion*/" /etc/apt/sources.list || sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu *UbuntuVersion*/" >> /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
```
##### Debian
```sh
grep -q -F "deb http://cran.rstudio.com/bin/linux/debian *DebianVersion*-*RVersion*/" /etc/apt/sources.list || sudo echo "deb http://cran.rstudio.com/bin/linux/debian *DebianVersion*-*RVersion*/" >> /etc/apt/sources.list


```sh
sudo apt update
sudo apt install r-base
```
#### Install RStudio
```sh
wget https://download1.rstudio.org/rstudio-*RStudioRelease*-amd64.deb
sudo dpkg -i rstudio-*RStudioRelease*-amd64.deb
sudo apt-get install -f  # install missing dependencies
```
#### Install Java
```sh
sudo apt-get install -y default-jre
sudo apt-get install -y default-jdk
sudo R CMD javareconf
```

If OpenSSL and libcurl libraries are not already installed, run these lines to do it.
```sh
sudo apt-get install libcurl4-openssl-dev
sudo apt-get install libssl-dev
```

## R Packages
You have then to install some packages before running the application. To do so, open RStudio and copy/paste the lines below in the console :
```R
install.packages(pkgs=c('rJava','devtools','R6','shiny','tiff','reshape','reshape2','RODBC','foreach',
'doParallel','stringi','rChoiceDialogs','gtools'), repos = "http://cloud.r-project.org")

source("http://bioconductor.org/biocLite.R")
biocLite(pkgs=c('EBImage','flowCore'), ask=F)

devtools::install_github(c("kroemerlab/MetaxpR","kroemerlab/ColocalizR"))
```
## Application
Once all is installed, you can run this line in the console to download and launch the app. 
```library(ColocalizR);Launcher()``` 

For more informations, please refer to the [User Manual](https://github.com/kroemerlab/ColocalizR/blob/master/ColocalizR%20-%20User%20Manual.pdf).
