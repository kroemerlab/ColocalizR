#' Set shiny progress bar in foreach loop
#'
#' @keywords internal
#' @import shiny
#' @import gtools
#' @export

pb.Combine <- function(){
  count <- 0
  function(...) {
    count <<- count + length(list(...)) - 1
    setProgress(value=count)
    Sys.sleep(0.01)
    flush.console()
    smartbind(...)
  }
}
