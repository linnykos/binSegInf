## Setup
outputdir = "output"
library(RColorBrewer)
library(igraph)
library(stringdist)
source('funs.R')

## Setting 1
n = 4
set.seed(0)
y = c(rep(0,n/2),rep(5,n/2)) + rnorm(n,0,1)
thresh = .1


## Setting 2
n = 60
sd = .5
mn = c(rep(-3, n/4), rep(2, n/4), rep(-1, n/4), rep(1, n/4))
set.seed(1)
y = mn + rnorm(n,0,sd)

slist = elist = blist = Blist = zlist = Zlist = matrix(NA, nrow = n, ncol = 2^8)
thresh = .5


## Do binary segmentation
binseg(s = 1, e = n, j = 0, k = 1, thresh = thresh, y = y, n = n)

#' Function to trim a matrix from the right and bottom, ridding of all-NA rows/columns.
trim = function(mat, type = c("rowcol","row")){
    type = match.arg(type)
    if(is.null(dim(mat))){ mat = mat[1:max(which(!is.na(mat)))]; return(mat)}
    last.j = max(which(!(apply(mat,1,function(myrow) return(all(is.na(myrow)))))))
    mat = mat[1:last.j,,drop=F]
    if(type=="rowcol"){
        last.j = max(which(!(apply(mat,2,function(mycol) return(all(is.na(mycol)))))))
        mat = mat[,1:last.j,drop=F]
    }
    return(mat)
}

## blist = trim(blist)
## zlist = trim(zlist)
## Blist = trim(Blist)
## Zlist = trim(Zlist)
## slist = trim(slist)
## elist = trim(elist)
   
bzlist = blist * zlist
basislist = list()
ii = 1
nrowGammat=0
for(my.j in 1:nrow(blist)){
    bsublist = blist[my.j,]
    bsublist = bsublist[!is.na(bsublist)]
    bsublist.allprev = collapse.prev(blist[1:my.j,])
    for(b in bsublist){
        bsublist.allprev = bsublist.allprev[bsublist.allprev!=b]
        all.other.b = unique(c(0,bsublist.allprev,n))
        my.se = get.closest(b, all.other.b)

        ## Calculate basis
        basislist[[ii]] = haarbasis(s = my.se[1], e = my.se[2], b = b, n = n, y = y, type="basis")
        
        ii=ii+1
    }
}


## Check orthogonality
check.orth.basis(basislist=basislist)


## Show results
basislist
blist = blist[1:last.j, 1:(2^last.j)]
zlist = zlist[1:last.j, 1:(2^last.j)]
bzlist = bzlist[1:last.j, 1:(2^last.j)]


## Prediction
X = do.call(cbind, basislist)
g = lm(y ~ X)
yhat = predict(g)


#####################
## Make data plot ###
#####################
tempcols = brewer.pal(3,"Set1")
pcol.dat = "grey50"
pch.dat = 16
lcol.mn = tempcols[1]
lty.mn = 1
lwd.mn = 2
lcol.pred = tempcols[2]
lty.pred = 1
lwd.pred = 2
lcol.bkpt = "grey80"
lty.bkpt = 2
lwd.bkpt = 1

ylim = range(y)*1.5
xlim = c(0,n)
xlab = "index"
ylab = ""
mar = c(4.5,4.5,1.5,1.5)
w = 7; h=5

pdf(file = file.path(outputdir, "example-data.pdf"), width = w, height = h)
par(mar=mar)
plot(NA, ylim = ylim, xlim = xlim, axes = F, xlab = xlab, ylab = ylab)
axis(1); axis(2);
points(y, col = pcol.dat, pch = pch.dat)
lines(mn, col = lcol.mn, lty = lty.mn, lwd = lwd.mn)
lines(yhat, col = lcol.pred, lty = lty.pred, lwd = lwd.pred)

## Make empty data frames for collecting data.
edgedat = data.frame(matrix(NA, nrow = 1000, ncol = 2))
nodedat = data.frame(matrix(NA, nrow = 1000, ncol = 2))

names(edgedat) = c("from", "to")
names(nodedat) = c("name", "value")

empty.nodedat = nodedat[0,]
empty.edgedat = edgedat[0,]


node.count = 0
set.seed(0)
ru = seq(from=-1.5,to=1.5,by=.3)
for(rownum in 1:nrow(blist)){
    my.pts = bzlist[rownum,]
    for(my.pt in my.pts[!is.na(my.pts)]){
        
        ii = which(my.pts == my.pt)
        node.count = node.count+1

        ## Add to plot
        abline(v = abs(my.pt), col = lcol.bkpt, lty = lty.bkpt, lwd = lwd.bkpt)
        bpower = paste0(rownum-1, ", ", ii)
        text(x = abs(my.pt), y = max(y)*1.5 + ru[node.count], label = bquote(b^{.(bpower)}))

        ## Add to edge data frame
        edgedat[node.count,"from"] = paste0(max(0,rownum-2),"_", floor((ii+1)/2))
        edgedat[node.count,"to"]   = paste0(rownum-1,"_", ii)
        ## edgedat[node.count,"text"]   = paste0("b[",(bpower),"]=", abs(my.pt), "\n s[",(bpower),"]=", (if(my.pt>0)"+" else "-" ))

        ## Add to node data frome
        nodedat[node.count,"name"] = paste0(rownum-1 ,"_", ii)
        nodedat[node.count,"value"] = paste0("b0[",(bpower),"]=", abs(my.pt), "\n z0[",(bpower),"]=", (if(my.pt>0)"+" else "-" ))
    }
}

edgedat = edgedat[1:node.count,]
nodedat = nodedat[1:node.count,]

legend("bottomright", lty = c(lty.mn, lty.pred), lwd = c(lwd.mn, lwd.pred),
       col = c(lcol.mn, lcol.pred), legend = c("Mean", "Standard Binseg"))
graphics.off()


## Helper function to extract j and k from node name (string)
extract.jk = function(nodename, what=c("j","k")){
    what = match.arg(what)
    where.underscore = as.numeric(gregexpr(pattern="_", nodename))
    if(what=="j"){
        return(as.numeric(substr(nodename, 1 , (where.underscore-1))))
    } else {
        return(as.numeric(substr(nodename, (where.underscore+1), nchar(nodename))))
    }
    ## Test cases
    ## extract.jk("11_556","j")
    ## extract.jk("11_556", "k")
}

node.count = 0

## Also add the children nodes directly after termination 
for(ii in 1:nrow(nodedat)){

    mynode = nodedat[ii,"name"]
    
    my.j = extract.jk(mynode, "j")
    my.k = extract.jk(mynode, "k")

    my.left.j = my.j+1
    my.left.k = 2*my.k-1

    my.right.j = my.j+1
    my.right.k = 2*my.k

    left.child.name  = paste0(my.left.j, "_", my.left.k) 
    right.child.name = paste0(my.right.j, "_", my.right.k)

    ## If left child doesn't exit
    if(!left.child.name %in% nodedat[,"name"]){
        cat("my node is", mynode, "adding left child", left.child.name, fill= T)
       node.count = node.count+1
        ## Add to node list and edge list as ZERO
       empty.edgedat[node.count,"from"] = paste0( my.j, "_", my.k)
       empty.edgedat[node.count,"to"]   = left.child.name
       empty.nodedat[node.count,"name"] = left.child.name
       empty.nodedat[node.count,"value"] = paste0("b0[", my.left.j, ",", my.left.k, "]=", -1)
    }

    ## If right child doesn't exit
    if(!right.child.name %in% nodedat[,"name"]){
       cat("my node is", mynode, "adding right child", right.child.name, fill= T)
       print(right.child.name)
       node.count = node.count+1
        ## Add to node list and edge list as ZERO
       empty.edgedat[node.count,"from"] = paste0(my.j, "_", my.k)
       empty.edgedat[node.count,"to"]   = right.child.name
       empty.nodedat[node.count,"name"] = right.child.name
       empty.nodedat[node.count,"value"] = paste0("b0[", my.right.j, ",", my.right.k, "]=", -1)
    }
}

nodedat = rbind(nodedat, empty.nodedat)
edgedat = rbind(edgedat, empty.edgedat)


################
## Plot Tree ###
################
w = 10; h = 7 
mar = c(1,1,1,1)
# create graph
g <- graph.data.frame(edgedat[-1,], directed = TRUE, vertices = nodedat)
V(g)$label.cex = 1
vv = V(g)$value 
vv = substr(vv, nchar(vv)-1, nchar(vv))
V(g)$color = sapply(stringdist(vv, "-1")==0, function(mydist) (if(mydist==0) "skyblue" else "pink"))

pdf(file = file.path(outputdir, "example-tree.pdf"), width = w, height = h)
par(mar=mar)
plot(g,  layout = layout.reingold.tilford,
     edge.label="", vertex.label = V(g)$value,
     vertex.size = 30, vertex.frame.color= "white", edge.arrow.size = .5)
graphics.off()


