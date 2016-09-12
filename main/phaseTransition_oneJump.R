rm(list=ls())
source("phaseTransition_source.R")
oneJumpRuleFl <- .oneJumpRuleClosure(200, fusedlasso1d)
firstJumpCriterion <- .extractJumpClosure(1)
trials <- 100
cores <- 20

jump.height <- seq(0, 5, length.out = 11)[-1]
jump.loc <- seq(0, 1, length.out = 11)[2:10]
paramMat <- as.matrix(expand.grid(0, jump.height, jump.loc))

resFl1Jump <- simulationGenerator(oneJumpRuleFl, paramMat, firstJumpCriterion,
  trials, cores)

save.image(file = paste0("res/oneJump-", Sys.Date(), ".RData"))
