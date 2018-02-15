## Synopsis: make the plot

outputdir = "../output"
filename = "uniform.pdf"
pdf(file=file.path(outputdir, filename), width=5, height=5)
mar = c(4,4,2,2)
nsim=1000
qqunif(runif(nsim,0,1), pch=16)
legend("topleft", pch=16, legend="null p-values")
graphics.off()
