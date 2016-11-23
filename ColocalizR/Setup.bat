@echo off
echo Installing R 3.3.2...
"3rd part\R-3.3.2-win.exe"
echo Installing Java...
if exist "C:\Program Files\R\R-3.3.2\bin\x64\R.exe" ("3rd part\jre-8u111-windows-x64.exe") ELSE ("3rd part\jre-8u111-windows-i586.exe")
echo Downloading and installing R packages...
if exist "C:\Program Files\R\R-3.3.2\bin\x64\R.exe" ("C:\Program Files\R\R-3.3.2\bin\x64\R.exe" CMD BATCH "3rd part\pkgs.R") ELSE ("C:\Program Files\R\R-3.3.2\bin\i386\R.exe" CMD BATCH "3rd part\pkgs.R")
echo Installation complete!
echo You can now use Launcher to start the app
pause
