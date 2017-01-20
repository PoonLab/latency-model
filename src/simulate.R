#library(rcolgem)
setwd('~/git/latency-model/src')
source('rcolgem-bd.R')
source('rong-perelson.R')


require(ade4)

# output destination folder
folder <- "../data/test/"




start.time <- 0
end.time <- 700  # time elapsed in units of days
  # note for longer simulation times it is probably better to use Brad's
  # steady state solution - the oscillations damped out rapidly and it is not 
  # worth carrying out detailed numerical solution of the ODE


#integ.method <- 'lsoda'
integ.method <- 'rk4'
fgy.resol <- 1e4  # this needs to be pretty high for longer simulation time frames (end.time)

# set initial conditions
x0 <- c(V=50, T=600, Ts=0.3, L=2)

tfgy <- make.fgy(start.time, end.time, births, deaths, nonDemeDynamics, x0, migrations=migrations, parms=params, fgyResolution=fgy.resol, integrationMethod=integ.method)

## It's a good idea to visually inspect the numerical solution of ODE first
# plot(tfgy[[5]])


ntimes <- 2  # number of time points we sample HIV RNA and DNA
nsamples.V <- 100  # number of HIV RNA (free virus) samples PER TIME POINT
nsamples.T <- 100  # number of cellular HIV DNA samples PER TIME POINT
nsamples <- nsamples.V + nsamples.T

# time series of number of infected cells
cells <- tfgy[[5]][, "L"] + tfgy[[5]][, "Ts"]

# initialize sample state matrix
sampleStates <- matrix(0, ncol=3, nrow=ntimes * nsamples)
colnames(sampleStates) <- demes

# assume uniform sampling over times
time.points <- seq(fgy.resol, 0, -ceiling(fgy.resol/ntimes))[1:ntimes]
#time.points <- (1:ntimes)*end.time/ntimes
sampleTimes <- tfgy[[5]][rep(time.points, each=nsamples), 1]
#sampleTimes <- rev(unlist(lapply((1:ntimes)*end.time/ntimes, rep, nsamples.V + nsamples.T)))

for (i in 1:length(time.points)) {
	tp <- time.points[i]
	row <- as.list(tfgy[[5]][tp,])  # extract time slice
	nsamples.L <- rbinom(1, nsamples.T, row$L/(row$L+row$Ts))  # number of latent cells in sample as binomial outcome
	nsamples.Ts <- nsamples.T - nsamples.L
	to.col <- c(rep(1, nsamples.V), rep(2, nsamples.L), rep(3, nsamples.Ts))
	
	# use (to.col) to assign 1's to random permutation of rows within time point block
	permut <- sample(1:nsamples, nsamples)
	for (j in 1:nsamples) {
		sampleStates[permut[j]+(i-1)*nsamples, to.col[j]] <- 1
	}
}



#cat(" building tree...")
n.reps <- 1  # number of trees to simulate

trees <- simulate.binary.dated.tree.fgy.wrapper(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, n.reps=n.reps, method='rk4')

# label tips -- by default they are only numbered
trees <- lapply(trees, function(tree) {
	tree$tip.label <- paste(tree$tip.label, apply(tree$sampleStates, 1, function(x) demes[which(x == 1)]), sep="."); 
	tree  # return
})
'multiPhylo' -> class(trees)

# write trees out to file
outfile <- paste0(folder, "aL-", params["a.L"], ".nwk")
write.tree(trees, file=outfile, append=FALSE)


