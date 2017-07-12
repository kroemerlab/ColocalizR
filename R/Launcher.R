#' Shiny App launcher
#' @export

Launcher = function(x, ...){

  shiny::runApp(appDir = system.file("shiny-app",package='ColocalizR'),...)
  
}