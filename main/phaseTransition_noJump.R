source("phaseTransition_source.R")

noJumpRuleFl <- .noJumpRuleClosure(100, fusedlasso1d)
noJumpRuleBs <- .noJumpRuleClosure(100, sbs)
noJumpCriterion <- .extractJumpClosure(1)
paramMat <- matrix(0, 1, 1)
trials <- 10

resFl <- simulationGenerator(noJumpRuleFl, paramMat, noJumpCriterion, trials)
resBs <- simulationGenerator(noJumpRuleBs, paramMat, noJumpCriterion, trials)

save.image(file = paste0("res/noJump-", Sys.Date(), ".RData"))
