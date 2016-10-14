## Define some helper functions
do.one.sim = function(delta,nsim,sigma,numsteps,n=12){
  p = rep(NA,nsim)
  theta = c(rep(0,n/2), rep(delta,n/2))
  for(isim in 1:nsim){
    y = theta + rnorm(n,0,sigma)
    a = binSegInf:::binseg.by.size(y, numsteps,verbose=FALSE)
    if(a$B[1]==n/2){
      v = binSegInf:::make.v(a$B[1],a$B,a$Z,n)
      p[isim]=binSegInf:::poly.pval(y,a$G,a$u,v,sigma)$pv
    }
  }
  return(p)
}

make.qq.line = function(p, pch='red', pcol=16){
  lcol.diag = "lightgrey"
  unif.p = runif(sum(!is.na(p)),0,1)
  a = qqplot(x=unif.p, y=p, plot.it=FALSE)
  points(x=a$y, y=a$x, col = pcol, pch=pch)
}

make.qqplot.background = function(){
  lcol.diag = "lightgrey"
  plot(NA,ylim=c(0,1),xlim=c(0,1), axes=F,  xlab="", ylab="")
  axis(1);axis(2)
  
  mtext("Observed",2,padj=-4)
  mtext("Expected",1,padj=4)
  abline(0,1,col=lcol.diag)
}

make.qqplot.legend = function(deltas,pch=16,pcols=rep("red",length(deltas))){
  legend("bottomright", col=pcols, pch=rep(pch,2),
         legend = sapply(c(bquote(delta == .(deltas[1])), 
                           bquote(delta == .(deltas[2])),
                           bquote(delta == .(deltas[3])),
                           bquote(delta == .(deltas[4]))), as.expression))
  title(main=bquote(atop(Segment~test~p-values)))
  
}