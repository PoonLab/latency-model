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

make.mod <- function(treatment.start, treatment.end, alpha, trans=10) {
	function(tp) {
		if (tp > treatment.start && tp < treatment.start + trans) {
			log(alpha)/trans
		} else if (tp > treatment.end && tp < treatment.end + trans) {
			-log(alpha)/trans
		} else {
			0
		}
	}
}

N <- 2000
delta <- 1
c <- 23
lambda <- 1e4
d.T <- 0.01

V.0 <- 453000
Ts.0 <- 5210

T <- 479000
E <- 2.4e-8

ntimes <- 10
nsamples <- 5
start.time <- 0
end.time <- 1000
f <- make.mod(10000, 10000, 1 - .85)

if (!is.null(params.file) && !is.na(params.file)) {
	params <- read.table(params.file, header=F, sep="\t")
	for (i in 1:(length(params$V1))) {
		eval(parse(text=paste(params$V1[i], params$V2[i], sep="=")))
	}
}

demes <- c("V", "Ts")

births <- t(rbind(c('0', 'parms$N * parms$delta * Ts'), c('E * T * V', '0')))
rownames(births) <- colnames(births) <- demes

migrations <- t(rbind(c('0', '0'), c('0', '0')))
rownames(migrations) <- colnames(migrations) <- demes

deaths <- c('parms$c * V', 'parms$delta * Ts')
names(deaths) <- demes

nonDemeDynamics <- c('parms$lambda - parms$d.T * T - E * V * T', 'f(t) * E')
names(nonDemeDynamics) <- c('T', 'E')

params.list <- list(N=N, delta=delta, c=c, lambda=lambda, d.T=d.T)
params.names <- names(params.list)
pop.list <- c(V=V.0, Ts=Ts.0, T=T, E=E)

if (F) {

dm <- build.demographic.process(births=births, deaths=deaths, migrations=migrations, nonDemeDynamics=nonDemeDynamics, parameterNames=params.names, rcpp=F, sde=F)

par(mfcol=c(1, 2))

tfgy <- show.demographic.process(dm, theta=params.list, x0=pop.list, t0=start.time, t1=end.time, res=2000, col=c("red", "green", "blue", "yellow"))
}

if (F) {
sampleTimes <- unlist(lapply((1:ntimes)*end.time/ntimes, rep, nsamples * 2))
sampleStates <- t(matrix(rbind(c(1, 0), c(0, 1)), nrow=2, ncol=nsamples * ntimes * 2))
colnames(sampleStates) <- demes
}

sampleTimes <- c(rep(0, 4), rep(1064, 16), rep(1309, 19), rep(1706, 24), rep(1887, 21))
sampleStates <- t(matrix(c(1, 0), nrow=2, ncol=4+16+19+24+21))
colnames(sampleStates) <- demes

tree <- sim.co.tree(params.list, dm, x0=pop.list, t0=start.time, sampleTimes=sampleTimes, sampleStates=sampleStates, res=2000)

if (F) {
tree$tip.label <- apply(tree$sortedSampleStates, 1, function(x) c("V", "T")[which(x == 1)])

cols <- unlist(apply(tree$sortedSampleStates, 1, function(x) c("red", "green")[which(x == 1)]))
}

#plot.phylo(tree, show.tip.label=F)
