trials <- 1000
contrast <- c(rep(1/5,5), rep(-1/5,5))
doMC::registerDoMC(cores = 3)

func <- function(i){
  if(i %% floor(trials/10) == 0) cat('*')
  set.seed(10*i)
  y <- rnorm(10)
  noise <- rnorm(10)
  sampler_kevin(y, noise, 1, 1, 100, contrast)
}

vec <- foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                          func(i))
vec <- unlist(vec)

# plot(sort(vec), seq(0,1,length.out = length(vec)))
# lines(c(0,1), c(0,1), col = "red", lwd = 2)

