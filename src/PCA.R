#library(rcolgem)
setwd('~/git/latency-model/src')
source('rcolgem-bd.R')

library(ade4)

folder <- "../data/pca.extreme/"




cat("initializing...\n")


###################################

## Define the model

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

deaths <- c('parms$c * V', 'parms$d.0 * L', 'parms$delta * Ts')
names(deaths) <- demes

# susceptible cells grow at constant rate (lambda), die at rate (d.T) 
#  and are depleted by infection
nonDemeDynamics <- c('parms$lambda - parms$d.T * T - parms$k * V * T')
names(nonDemeDynamics) <- c('T')


###################################

# read parameter limits from file
params.table <- read.table(paste0(folder, "params.txt"), col.names=c("names", "init", "less", "greater"), row.names=1)

params.names <- row.names(params.table)

# convert data frame into associative lists
params.list <- params.table$init
names(params.list) <- params.names
params.list <- lapply(params.list, function(x) x)

params.list.less <- params.table$less
names(params.list.less) <- params.names
params.list.less <- lapply(params.list.less, function(x) x)

params.list.greater <- params.table$greater
names(params.list.greater) <- params.names
params.list.greater <- lapply(params.list.greater, function(x) x)


n.reps <- 5
ntimes <- 10
nsamples.V <- 40
nsamples.T <- 60
start.time <- 0
end.time <- 1000

n <- length(unlist(params.list))
n.trials <- 1 + 2 * n  # every parameter is modified twice from baseline
#n.trials <- 1 + 2 * n + 4 * n * (n - 1) / 2

var.tab <- as.data.frame(t(matrix(unlist(params.list), nrow=n, ncol=n.trials)))
names(var.tab) <- params.names

for (i in 2:n.trials) {
	if (i <= n + 1)
		var.tab[i, i - 1] <- params.list.less[i - 1]
	else if (i <= 2 * n + 1)
		var.tab[i, i - n -  1] <- params.list.greater[i - n - 1]
}


# for (i in 2:n.trials) {
	# if (i <= n + 1)
		# var.tab[i, i-1] <- var.tab[i, i - 1] * 1.1
	# else if (i <= 2 * n + 1)
		# var.tab[i, i - n -  1] <- var.tab[i, i - n - 1] * 0.9
	# else if (i <= 3 * n + 1)
		# var.tab[i, i - 2 * n -  1] <- var.tab[i, i - 2 * n - 1] * 1.25
	# else
		# var.tab[i, i - 3 * n -  1] <- var.tab[i, i - 3 * n - 1] * 0.75
# }



# combs <- combn(seq(1, n), 2)

# for (i in 2:n.trials) {
	# if (i %% 1000 == 0)
		# cat(paste0(i, " "))
		
	# if (i <= n + 1)
		# var.tab[i, i - 1] <- var.tab[i, i - 1] * 1.1
	# else if (i <= 2*n + 1)
		# var.tab[i, i - n -  1] <- var.tab[i, i - n - 1] * 0.9
	# else {
		# r <- i - 2 * n - 2
	
		# s.j <- r %/% (n * (n - 1))
		# s.k <- r %% (n * (n - 1)) %/% (n * (n-1) / 2)
		# v <- combs[, r %% (n * (n - 1) / 2) + 1]
				
		# var.tab[i, v[1]] <- var.tab[i, v[1]] * (1.1 - .2 * s.j)
		# var.tab[i, v[2]] <- var.tab[i, v[2]] * (1.1 - .2 * s.k)
	# }
# }


cat("\n")


get.steady.state <- function(params) {
	V.0 <- with(params, N * lambda / c * (1 - d.0 / (d.0 + a.L) * eta) - d.T / k)
	T.0 <- with(params, lambda / (d.T + k * V.0))
	L.0 <- with(params, eta * k * V.0 * T.0 / (d.0 + a.L))
	Ts.0 <- with(params, c * V.0 / (N * delta))
	
	c(V=V.0, T=T.0, L=L.0, Ts=Ts.0)
}


suppress <- lapply(1:n.trials, function(i) {
		cat(i)
		p <- var.tab[i,]

		pop.list <- get.steady.state(p)

		tfgy <- make.fgy(start.time, end.time, births, deaths, nonDemeDynamics, pop.list, migrations=migrations, parms=p, fgyResolution=1000, integrationMethod="lsoda")

		cells <- tfgy[[5]][, "L"] + tfgy[[5]][, "Ts"]
		sampleStates <- t(matrix(c(1, 0, 0), nrow=3, ncol=ntimes * (nsamples.V + nsamples.T)))
		sampleTimes <- rev(unlist(lapply((1:ntimes)*end.time/ntimes, rep, nsamples.V + nsamples.T)))

		for (j in 0:(ntimes-1)) {
			sampleStates[j * (nsamples.T + nsamples.V) + 1:nsamples.T,] <- t(rmultinom(n=nsamples.T, size=1, prob=c(0, tfgy[[5]][1000 / ntimes * (j + 1), "L"]/cells[1000 / ntimes * (j + 1)], tfgy[[5]][1000 / ntimes * (j + 1), "Ts"]/cells[1000 / ntimes * (j + 1)])))
		}

		colnames(sampleStates) <- demes

		cat(" building tree...")

		trees <- simulate.binary.dated.tree.fgy.wrapper(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, n.reps=n.reps)

		trees <- lapply(trees, function(tree) {tree$tip.label <- paste(tree$tip.label, apply(tree$sampleStates, 1, function(x) c("V", "L", "T")[which(x == 1)]), sep="."); tree})

		lapply(1:n.reps, function(j) write.tree(trees[[j]], sprintf("%stree.%04d.%02d.tre", folder, i, j)))
		cat("\n")
	})
