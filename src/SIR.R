library(rcolgem)
library(phytools)

demes <- c('I')
births <- rbind(c('parms$beta*S*I / (S+I)'))
migrations <- rbind(c('0'))
deaths <- c(I='(parms$mu+parms$gamma)*I')
rownames(births) <- colnames(births) <- demes
rownames(migrations) <- colnames(migrations) <- demes

nonDemeDynamics <- c(S='-parms$mu*S + parms$mu*S + (parms$mu+parms$gamma)*I-S*(parms$beta*I) / (S+I)')

params <- list(beta=0.01, gamma=1/520, mu=1/3640)

pop.list <- c(S=1000 - 1, I=1)

start.time <- 0
end.time <- 30. * 52

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

pdf("../data/SIR/SIR.pdf", height=12, width=12)

show.demographic.process(tfgy, col=c('red', 'green', 'blue'))

suppress <- dev.off()

sampleStates <- as.matrix(data.frame(I=rep(1, floor(tfgy[[5]][length(tfgy[[5]][,"I"]), "I"] * 0.1))))
#names(sampleStates) <- c('I')
sampleTimes <- rep(1000, nrow(sampleStates))
names(sampleTimes) <- NULL

cat("making trees...")

trees <- simulate.binary.dated.tree.fgy(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod="lsoda", n.reps=50)

cat(length(trees))

#trees <- lapply(trees, function(tree) {tree$tip.label <- paste(tree$tip.label, apply(tree$sampleStates, 1, function(x) c("I")[which(x == 1)]), sep="."); tree})

suppress <- lapply(1:length(trees), function(i) write.tree(trees[[i]], file=sprintf("../data/SIR/tree.%03d.tre", i)))

