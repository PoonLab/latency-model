library(rcolgem)
#library(phydynR)

args <- commandArgs(trailing=T)

#params.file = NA

if (length(args > 0)) {
	params.file <- args[1]
}

show.demographic.process <- function (demo.model, theta, x0, t0, t1, res=1000, integrationMethod="lsoda", ...) {
    tfgy <- demo.model(theta, x0, t0, t1, res = res, integrationMethod = integrationMethod)
    o <- tfgy[[5]]
    t <- o[, 1]
    if (((ncol(o) - 1) == 2) & tail(colnames(o), 1) == "V2") {
        plot(t, o[, 2], type = "l", xlab = "Time", ylab = colnames(o)[2], 
            ...)
    }
    else {
        matplot(t, log(o[, 2:ncol(o)]) / log(10), type = "l", xlab = "Time", 
            ylab = "", ...)
             pch = 1, horiz = TRUE, ...)
    }
    
    tfgy
}

do.mod <- function(t.0) {
	if (t.0 > 249 && t.0 < 255) {
		log(1-.7)/10
	} else if (t.0 > 445 && t.0 < 455) {
		-log(1-.7)/10
	} else {
		0
	}
}

V.0 <- 7.5e5
L.0 <- 10 * 1e4
delta <- 1
c.0 <- 23
f <- 3e-6
N <- 2e4
E <- (1 - 0) * 2.248e-6
Ts.0 <- c.0 * V.0 / (N * delta)
T = 595
a <- 0.0058
am <- 0
omega <- 0.00939
r <- -0.00171
delta.l <- -r
p.bs <- r + delta.l
end.time <- 1000
mod <- 0
t.0 <- 0
ntimes <- 5
nsamples <- 5
lambda <- 1
d <- 1
p <- 1
Tmax <- 1000 

if (!is.na(params.file)) {
	params <- read.table(params.file, header=F, sep="\t")
	for (i in 1:(length(params$V1))) {
		eval(parse(text=paste(params$V1[i], params$V2[i], sep="=")))
	}
}

demes <- c("V", "L", "Ts")

births <- t(rbind(c('0', '0', 'parms$N * parms$delta * Ts'), c('parms$f * E * parms$T * V', 'parms$p.bs * L', '0'), c('(1 - parms$f) * E * parms$T * V', 'a * L', '0')))
rownames(births) <- colnames(births) <- demes

migrations <- t(rbind(c('0', '0', '0'), c('0', '0', '0'), c('0', '0', '0')))
rownames(migrations) <- colnames(migrations) <- demes

deaths <- c('parms$c.0  * V', 'parms$delta.l * L + a * L', 'parms$delta * Ts')
names(deaths) <- demes

nonDemeDynamics <- c('parms$omega * (parms$am - a)', 'parms$mod * do.mod(t.0) * E', '1') #, 'parms$lambda - parms$d * T + parms$p * T * (1 - T/parms$Tmax) - parms$E*V*T')
names(nonDemeDynamics) <- c('a', 'E', 't.0')#, 'T')
#nonDemeDynamics <- c('log(1/0.01)/100 * E') 
#names(nonDemeDynamics) <- c('E')

params.list <- list(delta=delta, N=N, f=f, c.0=c.0, delta.l=delta.l, p.bs=p.bs, omega=omega, am=am, mod=mod, lambda=lambda, d=d, p=p, Tmax=Tmax, T=T)
params.names <- c('delta', 'N', 'f', 'c.0', 'delta.l', 'p.bs', 'omega', 'am', 'mod', 'lambda', 'd', 'p', 'Tmax', 'T')
pop.list <- c(V=V.0, L=L.0, Ts=Ts.0, a=a, E=E, t.0=t.0)#, T=T)

dm <- build.demographic.process(births=births, deaths=deaths, nonDemeDynamics=nonDemeDynamics, parameterNames=params.names, rcpp=F, sde=F)

par(mfcol=c(1, 2))

tfgy <- show.demographic.process(dm, theta=params.list, x0=pop.list, t0=0, t1=end.time, col=seq(2, 5))

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
