
levs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)[8:9]
results.by.lev = list()
mc.cores=8
nsim=3000
nsims=c(seq(from=3000,to=1000,length=5), round(seq(from=600, to=300, length=4) ))[8:9]
visc=visc.fourjump = unlist(lapply(c(1,2,3,4)*(n/5), function(cp)cp+c(-1,0,1)))
n=200

## Compare p-values from three methods
dosim_flrandplus <- function(type=c("wbs","fl.nonrand","fl.rand","fl.rand.plus","sbs.rand",
                                 "sbs.nonrand", "wbs.rand", "wbs.nonrand",
                                 "cbs.rand", "cbs.nonrand"),
                          n, lev, numIntervals=n, sigma.add=0.2, numIS=100, meanfun=onejump, visc=NULL, numSteps=1, bits=50){

    type = match.arg(type)
    if(is.null(visc))visc=1:n
    sigma = 1
    mn = meanfun(lev=lev,n=n)
    y = mn + rnorm(n, 0, sigma)
    cumsum.y = cumsum(y)
    inference.type = "pre-multiply"
    improve.nomass.problem = TRUE
    if(visc==NULL) retain=1:n ## default value.


    if(type=="fl.rand.plus"){


        ## Draw new noise
        new.noise = rnorm(n,0,sigma.add)

        ## 2. Get stopped FL p-values
        numSteps = 10
        consec=2
        g.fudged = dualpathSvd2(y + new.noise, maxsteps=numSteps, D = makeDmat(n,ord=0))
        ic_obj = get_ic(g.fudged$cp, g.fudged$y, consec=consec, sigma=sigma, type="bic")
        stoptime = ic_obj$stoptime

        ## Get stopped polyhedron
        if(ic_obj$flag=="normal" ){

            ## Get model selection event polyhedron
            Gobj.fudged = genlassoinf::getGammat.naive(obj=g.fudged, y=y, condition.step=stoptime+consec)
            poly.fudged = polyhedra(obj=Gobj.fudged$G, u=Gobj.fudged$u)

            ## Get ic-stopped model selection polyhedron
            ic_poly = ic_obj$poly

            ## Combine them
            combined.poly = polyhedra(obj = rbind(poly.fudged$gamma, ic_poly$gamma),
                                      u = c(poly.fudged$u, ic_poly$u))

            ## Postprocess and retain vicinity contrasts
            cp = g.fudged$cp[1:stoptime]
            cp.sign = g.fudged$cp.sign[1:stoptime]
            cp = declutter(coords=cp, coords.sign=cp.sign, how.close=3)
            cp.sign = sign(cp)
            cp = abs(cp)
            vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
            retain = which(abs(as.numeric(names(vlist))) %in% visc)
            vlist = vlist[retain]

            ## Collect the p-values
            pvs = sapply(vlist, function(v){
                pv = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                        sigma.add=sigma.add,
                                        orig.fudged.poly= combined.poly)
            })
            locs = as.numeric(names(vlist))
            return(data.frame(pvs=pvs, locs=locs))
        } else {
            return(data.frame(pvs=NA, locs=NA))
        }
    }
