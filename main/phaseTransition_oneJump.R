source("phaseTransition_source.R")

noJumpRuleFl <- .noJumpRuleClosure(100, fusedlasso1d)
noJumpRuleBs <- .noJumpRuleClosure(100, sbs)
noJumpCriterion <- .extractJumpClosure(1)
paramMat <- matrix(0, 1, 1)
trials <- 10

resFl <- simulationGenerator(noJumpRuleFl, paramMat, noJumpCriterion, trials)
resBs <- simulationGenerator(noJumpRuleBs, paramMat, noJumpCriterion, trials)



oneJumpRuleFl <- .oneJumpRuleClosure(100, fusedlasso1d)
paramMat <- matrix(c(0,1,.5,0,1,.2), ncol = 3, nrow = 2, byrow = T)

resFl1Jump <- simulationGenerator(oneJumpRuleFl, paramMat, noJumpCriterion,
  trials)


