## Read in data
a = read.table("../data/PAGED_p_TCGA_182_OvGBM_SNP_N_GenomeWideSNP_6_D04_896216.grch38.seg.txt",header=TRUE)
b = read.table("../data/PAGED_p_TCGA_182_OvGBM_SNP_N_GenomeWideSNP_6_D04_896216.nocnv_grch38.seg.txt", header=TRUE)

par(mfrow=c(2,3))
for(chromnum in 1:6){

    adat = a[a[,"Chromosome"]%in%chromnum,]
    adat[,"Chromosome"]
    plot(adat[,"Segment_Mean"])
    ## abline(v=cumsum(sapply(1:22, function(ii) sum(a[,"Chromosome"]==ii))))
    
    
    ## bdat
    bdat = b[b[,"Chromosome"]%in%chromnum,]
    bdat[,"Chromosome"]
    
    plot(adat[,"Segment_Mean"]~adat[,"Start"])
    points(bdat[,"Segment_Mean"]~bdat[,"Start"],col='red')

    ## Plot segments
    plot(NA, xlim = range(bdat[,"Start"]), ylim = c(-4,4))
    Map(cnvpoint, adat[,"Segment_Mean"], adat[,"Start"], adat[,"End"], rep("black",nrow(adat)))
    Map(cnvline, adat[,"Segment_Mean"], adat[,"Start"], adat[,"End"], rep("black",nrow(adat)))
    Map(cnvpoint, bdat[,"Segment_Mean"], bdat[,"Start"], bdat[,"End"], rep("red",nrow(bdat)))
    Map(cnvline, bdat[,"Segment_Mean"], bdat[,"Start"], bdat[,"End"], rep("red",nrow(bdat)))
    
}

plot(a[,"Start"])

cnvline <- function(y,start,end,col){
    segments(x0=start,y0=y,x1=end,y1=y,col=col)
}

cnvpoint <- function(y,start,end,col){
    points(x = mean(start:end),y=y,col=col)
}
