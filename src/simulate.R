#library(rcolgem)
setwd('~/git/latency-model/src')
source('rcolgem-bd.R')

require(ade4)

# output destination folder
folder <- "../data/test/"




cat("initializing...\n")


###################################

## Define the model


## virus lineages can be sampled from:
# V = free virus
# L = latently infected cells
# Ts = active infected cells
demes <- c("V", "L", "Ts")


births <- t(rbind(
	# infected cells burst at rate (delta) and produce (N) virions
	c('0', '0', 'parms$N * parms$delta * Ts'), 
	
	# (k) is infection rate of susceptible cells (unsampled compartment T)
	# goes latent with probability (eta)
	c('parms$eta * parms$k * T * V', '0', '0'), 
	
	# otherwise infection remains active
	c('(1 - parms$eta) * parms$k * T * V', '0', '0'))
)
rownames(births) <- colnames(births) <- demes

# latent cells reactive at rate (a.L)
migrations <- t(rbind(c('0', '0', '0'), c('0', '0', '0'), c('0', 'parms$a.L * L', '0')))
rownames(migrations) <- colnames(migrations) <- demes

# free viruses removed at rate (c)
# latently-infected cells removed at rate (d.0)
# active infected cells removed at rate by bursting (delta)
deaths <- c('parms$c * V', 'parms$d.0 * L', 'parms$delta * Ts')
names(deaths) <- demes

# susceptible cells grow at constant rate (lambda), die at rate (d.T) 
#  and are depleted by infection
nonDemeDynamics <- c('parms$lambda - parms$d.T * T - parms$k * V * T')
names(nonDemeDynamics) <- c('T')


###################################


# default model parameters, taken from Rong and Perelson (2009; PLOS Comput Biol e1000533)
params <- list()
params["lambda"] <- 1e4  # growth rate of uninfected cells (per mL per day)
params["d.T"] <- 0.01    # death rate of uninfected cells (per day)
params["k"] <- 2.4e-8    # rate of infection (mL/day)
params["eta"] <- 0.01   # probability of entering latent state
params["d.0"] <- 0.001   # death rate of latently-infected cells
params["a.L"] <- 0.2     # rate of transition from latently to productively infected cells
params["delta"] <- 1.0   # death rate of productively infected cells (per day)
params["N"] <- 2000      # number of virions produced by cell death
params["c"] <- 23.        # clearance rate of free virus (per day)


# # get.steady.state <- function(params) {
	## Equation (2) from paper (epsilon = 0)
	# V.0 <- with(params, N * lambda / c * (1 - d.0 / (d.0 + a.L) * eta) - d.T / k)
	# T.0 <- with(params, lambda / (d.T + k * V.0))
	# L.0 <- with(params, eta * k * V.0 * T.0 / (d.0 + a.L))
	# Ts.0 <- with(params, c * V.0 / (N * delta))
	
	# c(V=V.0, T=T.0, L=L.0, Ts=Ts.0)
# }


start.time <- 0
end.time <- 700  # time elapsed in units of days


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
n.reps <- 20  # number of trees to simulate

trees <- simulate.binary.dated.tree.fgy.wrapper(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, n.reps=n.reps, method='rk4')

# label tips -- by default they are only numbered
trees <- lapply(trees, function(tree) {
	tree$tip.label <- paste(tree$tip.label, apply(tree$sampleStates, 1, function(x) demes[which(x == 1)]), sep="."); 
	tree  # return
})
'multiPhylo' -> class(trees)

# write trees out to file
outfile <- paste0(folder, "aL-0.2.nwk")
write.tree(trees, file=outfile, append=FALSE)


