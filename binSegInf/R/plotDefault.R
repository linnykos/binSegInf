defaultColors <- function(vec = NA, transparency = 0){
  stopifnot(is.numeric(vec))
  stopifnot(is.numeric(transparency), length(transparency) == 1)
  
  if(all(is.na(vec))) vec <- .defaultColorRange()
  trans <- round(transparency / 255)
  
  color.numeric <- .defaultColorRange()
  color.rgb <- c(
    grDevices::rgb(0, 0, 0, trans, max = 255), #black
    grDevices::rgb(255, 0, 0, trans, max = 255), #red
    grDevices::rgb(0, 205, 0, trans, max = 255), #green3
    grDevices::rgb(0, 0, 255, trans, max = 255)) #blue
  
  plyr::mapvalues(vec, from = color.numeric, to = color.rgb, warn_missing = F)
}

.defaultColorRange <- function(){
  1:4
}

.demoColors <- function(){
  col.vec <- defaultColors()
  defaultPlotDefaults()
  graphics::plot(x = 1:length(col.vec), y = rep(1, length(col.vec)), 
    col = col.vec, cex = 4)
}


#' Default Font
#'
#' @return a character of the font
#' @export
defaultFont <- function(){
  "sans"
}

#' Default Graphical Parameter (Cex) 
#'
#' @return a numeric of the cex
#' @export
defaultCexDefault <- function(){1.5}

#' Default Graphical Parameter (Pch)
#'
#' @return a numeric of the pch
#' @export
defaultPchDefault <- function(){16}


defaultPlotDefaults <- function(){
  res <- names(grDevices::dev.cur())
  if(res != "null device" & res != "pdf") grDevices::graphics.off()
  graphics::par(pch = defaultPchDefault(), cex = defaultCexDefault(), 
    family = defaultFont(), mar = c(4,4,1,1))
  
  invisible()
}
