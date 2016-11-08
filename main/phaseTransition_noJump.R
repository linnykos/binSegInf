source("phaseTransition_source.R")

noJumpRuleFl <- .noJumpRuleClosure(200, fusedlasso1d)
noJumpRuleBs <- .noJumpRuleClosure(200, sbs)
firstJumpCriterion <- .extractJumpClosure(1)
paramMat <- matrix(0, 1, 1)
trials <- 500
cores <- 20

resFl <- simulationGenerator(noJumpRuleFl, paramMat, firstJumpCriterion, trials, cores)
resBs <- simulationGenerator(noJumpRuleBs, paramMat, firstJumpCriterion, trials, cores)

save.image(file = paste0("res/noJump-", Sys.Date(), ".RData"))
quit()
