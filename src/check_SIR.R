library(rcolgem)

births <- rbind(c('parms$beta * I * S'))
deaths <- c(I='parms$gamma * I')
rownames(births) <- colnames(births) <- c('I')

nonDemeDynamics <- c(S='-parms$beta * I * S', R='parms$gamma * I')

params <- list(gamma=0.001, beta=0.0001)

pop.list <- c(S=9999, I=1, R=0)

start.time <- 0
end.time <- 1000

tfgy <- make.fgy(start.time, end.time, births, deaths, nonDemeDynamics, pop.list, parms=params, fgyResolution=1000, integrationMethod = "adam")

show.demographic.process <- function(tfgy, ...) {
#    tfgy <- demo.model(theta, x0, t0, t1, res = 1000, integrationMethod = integrationMethod)
    o <- tfgy[[5]]
    t <- o[, 1]
    matplot(t, o[, 2:ncol(o)], type = "l", xlab = "Time", ylab = "", ...)
    legend("bottomright", inset = 0.05, legend = colnames(o)[2:ncol(o)], pch = 1, horiz = TRUE, ...)
#    tfgy
}

show.demographic.process(tfgy, col=c('red', 'green', 'blue'))

sampleStates <- t(t(rep(0, 100)))
colnames(sampleStates) <- c('I')
sampleTimes <- rep(1000, 100)

debug(simulate.binary.dated.tree.fgy)

trees <- simulate.binary.dated.tree.fgy(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod="lsoda", n.reps=100)

LTT <- function(tree) {

}