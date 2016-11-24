##' makes QQ line on existing plot.
make.qq.line = function(p, pcol='red', pch=16){
  ## set.seed(0)
  ## unif.p = runif(sum(!is.na(p)),0,1)
  unif.p = seq(from=0,to=1,length=sum(!is.na(p)))
  a = qqplot(x=unif.p, y=p, plot.it=FALSE)
  points(x=a$y, y=a$x, col = pcol, pch=pch)
}


##' makes QQ plot background without any lines.
make.qqplot.background = function(plot.title=""){
  lcol.diag = "lightgrey"
  plot(NA,ylim=c(0,1),xlim=c(0,1), axes=F,  xlab="", ylab="")
  axis(1);axis(2)
  
  ## Embellishments
  mtext("Observed",2,padj=-4)
  mtext("Expected",1,padj=4)
  abline(0,1,col=lcol.diag)
  title(plot.title)
}

##' makes QQ plot legend on existing plot.
make.qqplot.legend = function(deltas,pch=16,pcols=rep("red",length(deltas)), where.legend="bottomright"){
  legend(x=where.legend, col=pcols, pch=rep(pch,2),
         legend = sapply(c(bquote(delta == .(deltas[1])), 
                           bquote(delta == .(deltas[2])),
                           bquote(delta == .(deltas[3])),
                           bquote(delta == .(deltas[4]))), as.expression))
  title(main=bquote(atop(Segment~test~p-values)))
  
}

