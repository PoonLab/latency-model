#!/usr/bin/Rscript -f
#PLoS Comput Biol. 2009 Oct;5(10):e1000533. doi: 10.1371/journal.pcbi.1000533. Epub 2009 Oct 16.

library(rcolgem)
#library(phydynR)

args <- commandArgs(trailing=T)

#params.file = NA

if (length(args > 0)) {
	params.file <- args[1]
}

show.demographic.process <- function(demo.model, theta, x0, t0, t1, res = 1000, integrationMethod = "lsoda", ...) 
{
    tfgy <- demo.model(theta, x0, t0, t1, res = 1000, integrationMethod = integrationMethod)
    o <- tfgy[[5]]
    t <- o[, 1]
    matplot(t, log10(o[, 2:ncol(o)]), type = "l", xlab = "Time", ylab = "", ...)
    legend("bottomright", inset = 0.05, legend = colnames(o)[2:ncol(o)], pch = 1, horiz = TRUE, ...)
    tfgy
}

make.mod <- function(treament.start, treatment.end, alpha, trans=10) {
	function(tp) {
		if (tp > treatment.start && tp < treament.start + trans) {
			log(alpha)/trans
		} else if (tp > treatment.end && tp < treatment.end + trans) {
			-log(alpha)/trans
		} else {
			0
		}
	}
}

# init params
lambda <- 1e4
d.T <- 0.01
eta <- 0.001
d.0 <- 0.001
a.L <- 0.1
delta <- 1
N <- 2000
c <- 23
r.L <- 0

V.0 <- 452860
L.0 <- 52
Ts.0 <- 5208

T <- 479178
E <- 2.4e-8
tp <- 0

ntimes <- 1
nsamples <- 5
end.time <- 2000
treatment.start <- 10000
treatment.end <- 10000
epsilon <- 0.85
do.mod <- make.mod(treatment.start, treatment.end, 1-epsilon)

if (!is.null(params.file) && !is.na(params.file)) {
	params <- read.table(params.file, header=F, sep="\t")
	for (i in 1:(length(params$V1))) {
		eval(parse(text=paste(params$V1[i], params$V2[i], sep="=")))
	}
}

demes <- c("V", "L", "Ts")

births <- t(rbind(c('0', '0', 'parms$N * parms$delta * Ts'), c('parms$eta * E * T * V', 'parms$r.L * L', '0'), c('(1 - parms$eta) * E * T * V', 'parms$a.L * L', '0')))
rownames(births) <- colnames(births) <- demes

migrations <- t(rbind(c('0', '0', '0'), c('0', '0', '0'), c('0', '0', '0')))
rownames(migrations) <- colnames(migrations) <- demes

deaths <- c('parms$c * V', 'parms$d.0 * L + parms$a.L * L', 'parms$delta * Ts')
names(deaths) <- demes

nonDemeDynamics <- c('parms$lambda - parms$d.T * T - E * V * T', 'do.mod(tp) * E', '1')
names(nonDemeDynamics) <- c('T', 'E', 'tp')

params.list <- list(lambda=lambda, d.T=d.T, eta=eta, d.0=d.0, a.L=a.L, delta=delta, N=N, c=c, r.L=r.L)
params.names <- names(params.list)
pop.list <- c(V=V.0, L=L.0, Ts=Ts.0, T=T, E=E, tp=tp)

dm <- build.demographic.process(births=births, deaths=deaths, nonDemeDynamics=nonDemeDynamics, parameterNames=params.names, rcpp=F, sde=F)

par(mfcol=c(1, 2))

tfgy <- show.demographic.process(dm, theta=params.list, x0=pop.list, t0=0, t1=end.time, col=c("red", "green", "blue", "yellow", "cyan", "magenta"))

cells <- tfgy[[5]][, "L"] + tfgy[[5]][, "Ts"]
sampleStates <-t(matrix(c(1, 0, 0), nrow=3, ncol=3 * ntimes * nsamples))

for (i in 0:(ntimes-1)) {
	for (j in 0:(nsamples-1)) {		
		sampleStates[i * nsamples * 3 + j * 3 + 1:2,] <- t(rmultinom(n=2, size=1, prob=c(0, tfgy[[5]][1000 / ntimes * (i + 1), "L"]/cells[1000 / ntimes * (i + 1)], tfgy[[5]][1000 / ntimes * (i + 1), "Ts"]/cells[1000 / ntimes * (i + 1)])))
	}
}

colnames(sampleStates) <- demes
sampleTimes <- unlist(lapply((1:ntimes)*end.time/ntimes, rep, 3 * nsamples))

tree <- sim.co.tree(params.list, dm, x0=pop.list, t0=0, sampleTimes=sampleTimes, sampleStates=sampleStates, res=2000)

tree$tip.label <- apply(tree$sortedSampleStates, 1, function(x) c("V", "L", "Ts")[which(x == 1)])

cols <- unlist(apply(tree$sortedSampleStates, 1, function(x) c("red", "green", "blue")[which(x == 1)]))

plot.phylo(tree, cex=.75, tip.color=cols)
