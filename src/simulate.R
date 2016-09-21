#library(rcolgem)
setwd('~/git/latency-model/src')
source('rcolgem-bd.R')

library(ade4)

folder <- "../data/pca.extreme/"




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


# default model parameters
params <- list()
params["lambda"] <- 1e4  # growth rate of uninfected cells
params["d.T"] <- 1e-2    # death rate of uninfected cells
params["k"] <- 2.4e-8    # rate of infection
params["eta"] <- 0.01   # probability of entering latent state
params["d.0"] <- 0.001   # death rate of latently-infected cells
params["a.L"] <- 0.01     # reactivation rate of latently-infected cells
params["delta"] <- 1.0   # death rate of active infected cells
params["N"] <- 2e3       # number of virions produced by cell death
params["c"] <- 5.4       # removal rate of free virus 


get.steady.state <- function(params) {
	V.0 <- with(params, N * lambda / c * (1 - d.0 / (d.0 + a.L) * eta) - d.T / k)
	T.0 <- with(params, lambda / (d.T + k * V.0))
	L.0 <- with(params, eta * k * V.0 * T.0 / (d.0 + a.L))
	Ts.0 <- with(params, c * V.0 / (N * delta))
	
	c(V=V.0, T=T.0, L=L.0, Ts=Ts.0)
}

n.trials <- 1
n.reps <- 3  # number of trees to simulate
start.time <- 0
end.time <- 50  # time elapsed in units of ???
ntimes <- 3  # number of time points we sample HIV RNA and DNA
nsamples.V <- 50  # number of HIV RNA (free virus) samples
nsamples.T <- 50  # number of cellular HIV DNA samples

#integ.method <- 'lsoda'
integ.method <- 'rk4'
fgy.resol <- 500

x0 <- c(V=10, T=1e5, Ts=0, L=0)



# initial population vector
#pop.list <- get.steady.state(p)
pop.list <- x0

tfgy <- make.fgy(start.time, end.time, births, deaths, nonDemeDynamics, pop.list, migrations=migrations, parms=params, fgyResolution=fgy.resol, integrationMethod=integ.method)

cells <- tfgy[[5]][, "L"] + tfgy[[5]][, "Ts"]
sampleStates <- t(matrix(c(1, 0, 0), nrow=3, ncol=ntimes * (nsamples.V + nsamples.T)))
sampleTimes <- rev(unlist(lapply((1:ntimes)*end.time/ntimes, rep, nsamples.V + nsamples.T)))

# assign deme sample states at random
for (j in 0:(ntimes-1)) {
	sampleStates[j * (nsamples.T + nsamples.V) + 1:nsamples.T,] <- t(rmultinom(n=nsamples.T, size=1, prob=c(0, tfgy[[5]][fgy.resol / ntimes * (j + 1), "L"]/cells[fgy.resol / ntimes * (j + 1)], tfgy[[5]][fgy.resol / ntimes * (j + 1), "Ts"]/cells[fgy.resol / ntimes * (j + 1)])))
}

colnames(sampleStates) <- demes

#cat(" building tree...")

trees <- simulate.binary.dated.tree.fgy.wrapper(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, n.reps=n.reps, method='rk4')

# label tips
trees <- lapply(trees, function(tree) {tree$tip.label <- paste(tree$tip.label, apply(tree$sampleStates, 1, function(x) c("V", "L", "T")[which(x == 1)]), sep="."); tree})

# write to files
lapply(1:n.reps, function(j) write.tree(trees[[j]], sprintf("%stree.%02d.tre", folder, j)))
#cat("\n")

