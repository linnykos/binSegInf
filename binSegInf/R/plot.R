plot.CpVector <- function(x, col.dots = grDevices::rgb(0.5,0.5,0.5), 
  col.line = 1, ...){
  
  graphics::plot(x$data, col = col.dots, ...)
  
  lis <- .splitChangepoints(length(x$data), x$jump.height, x$jump.idx)
  
  for(i in 1:length(lis)){
    graphics::lines(lis[[i]]$x, lis[[i]]$y, col = col.line, ...)
  }
  
  invisible()
}

.splitChangepoints <- function(n, jump.height, jump.idx){
  jump.idx <- c(0, jump.idx, n)
  
  lis <- vector("list", length(jump.idx) - 2)
  for(i in 2:(length(jump.idx))){
    x <- (jump.idx[i-1]+1):jump.idx[i]
    y <- rep(jump.height[i-1], length(x))
    lis[[i-1]] <- list(x = x, y = y)
  }
  
  lis
}