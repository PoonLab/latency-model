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
phi <- 0.1

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

pdf("../data/VolzHIV/VolzHIV.pdf", height=12, width=12)

show.demographic.process(tfgy, col=c('red', 'green', 'blue'))

dev.off()

Y <- tfgy[[5]]
sample.index <- which(Y[, "S"] < 1)[1]

sampleStates <- t(matrix(c(rep(c(1, 0), floor(Y[sample.index, "I.1"] * phi)), rep(c(0, 1), floor(Y[sample.index, "I.2"] * phi))), nrow=2))
colnames(sampleStates) <- c('I.1', 'I.2')
sampleTimes <- rep(Y[sample.index, "time"], nrow(sampleStates))
names(sampleTimes) <- NULL

cat("making trees...")

trees <- simulate.binary.dated.tree.fgy(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod="lsoda", n.reps=50)

cat(length(trees))

trees <- lapply(trees, function(tree) {tree$tip.label <- paste(tree$tip.label, apply(tree$sampleStates, 1, function(x) c("I.1", "I.2")[which(x == 1)]), sep="."); tree})

suppress <- lapply(1:length(trees), function(i) write.tree(trees[[i]], file=sprintf("../data/VolzHIV/tree.%03d.tre", i)))

