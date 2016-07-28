library(rcolgem)
library(phytools)

demes <- c('I.1', 'I.2')
births <- rbind(c('S / parms$N * parms$beta.1 * I.1', '0'), c('S / parms$N * parms$beta.2 * I.2', '0'))
migrations <- rbind(c('0', 'parms$gamma.1 * I.1'), c('0', '0'))
deaths <- c(I.1='0', I.2='parms$gamma.2 * I.2')
rownames(births) <- colnames(births) <- demes
rownames(migrations) <- colnames(migrations) <- demes

nonDemeDynamics <- c(S='-S / parms$N * (parms$beta.1 * I.1 + parms$beta.2 * I.2)')

R.0 <- 20
q <- 0.42
N <- 1000

params <- list(gamma.1=1/365,  gamma.2=1/3650, beta.1=R.0*q/(1+q)/3650, beta.2=R.0/(1+q)/3650, N=N)

pop.list <- c(S=N - 2, I.1=1, I.2=1)

start.time <- 0
end.time <- 10000

cat("making fgy...")

tfgy <- make.fgy(start.time, end.time, births, deaths, nonDemeDynamics, pop.list, parms=params, migrations=migrations, fgyResolution=10000, integrationMethod = "lsoda")

show.demographic.process <- function(tfgy, ...) {
#    tfgy <- demo.model(theta, x0, t0, t1, res = 1000, integrationMethod = integrationMethod)
    o <- tfgy[[5]]
    t <- o[, 1]
    matplot(t, o[, 2:ncol(o)], type = "l", xlab = "Time", ylab = "", ...)
    legend("bottomright", inset = 0.05, legend = colnames(o)[2:ncol(o)], pch = 1, horiz = TRUE, ...)
#    tfgy
}

par(mfcol=c(1, 2))

show.demographic.process(tfgy, col=c('red', 'green', 'blue'))

Y <- tfgy[[5]]
sample.index <- which(Y[, "S"] < 1)[1]

sampleStates <- t(matrix(c(rep(c(1, 0), floor(Y[sample.index, "I.1"])), rep(c(0, 1), floor(Y[sample.index, "I.2"]))), nrow=2))
colnames(sampleStates) <- c('I.1', 'I.2')
sampleTimes <- rep(Y[sample.index, "time"], nrow(sampleStates))
names(sampleTimes) <- NULL

cat("making trees...")

trees <- simulate.binary.dated.tree.fgy(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod="lsoda", n.reps=100)

cat(length(trees))

ltts <- lapply(trees, function(tree) {l <- ltt(tree, plot=F, log.lineages=F); l$times <- sampleTimes[1] - l$times[length(l$times)] + l$times; l})

suppress <- lapply(ltts, function(l) {lines(l$times, l$ltt, col=rgb(0, 0, 0, .2))})

detach("package:rcolgem")

library(phydynR)

dm <- build.demographic.process(births=births, deaths=deaths, migrations=migrations, nonDemeDynamics=nonDemeDynamics, parameterNames=names(params), rcpp=F, sde=F)

cat("making fgy...")

tfgy.phy <- dm(params, pop.list, start.time, end.time, res = 10000, integrationMethod = "lsoda")

show.demographic.process(tfgy.phy, col=c('red', 'green', 'blue'))

Y <- tfgy.phy[[5]]
sample.index <- which(Y[, "S"] < 1)[1]

sampleStates <- t(matrix(c(rep(c(1, 0), floor(Y[sample.index, "I.1"])), rep(c(0, 1), floor(Y[sample.index, "I.2"]))), nrow=2))
colnames(sampleStates) <- c('I.1', 'I.2')
sampleTimes <- rep(Y[sample.index, "time"], nrow(sampleStates))
names(sampleTimes) <- NULL

cat("making trees...")

trees.phy <- lapply(1:100, function(x) {cat(" "); cat(x); sim.co.tree.fgy(tfgy.phy, sampleTimes, sampleStates)})

ltts.phy <- lapply(trees.phy, function(tree) {l <- ltt(tree, plot=F, log.lineages=F); l$times <- sampleTimes[1] - l$times[length(l$times)] + l$times; l})

suppress <- lapply(ltts.phy, function(l) {lines(l$times, l$ltt, col=rgb(0, 0, 0, .2))})