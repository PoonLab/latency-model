###################################

## Define the model


## virus lineages can be sampled from:
# V = free virus
# L = latently infected cells
# Ts = productively infected cells
demes <- c("V", "L", "Ts")

# TODO: clonal expansion of T cells
births <- rbind(
	# productively infected cells burst at rate (delta) and produce (N) virions
	c('0', '0', '0'),
	c('0', '0', '0'), 
	c('parms$N*parms$delta*Ts', '0', '0')  # Ts->V
)
rownames(births) <- colnames(births) <- demes


# (k) is infection rate of susceptible cells
# (eta) is probability of latent infection (V->L)
# (a.L) reactivation rate of latently-infected cells (L->Ts)
migrations <- rbind(
	c('0', 'parms$eta*parms$k*T*V', '(1-parms$eta)*parms$k*T*V'), 
	c('0', '0', 'parms$a.L * L'), 
	c('0', '0', '0')
)
rownames(migrations) <- colnames(migrations) <- demes


deaths <- c(
	'parms$c * V',      # free viruses removed at rate (c)
	'parms$d.0 * L',    # latently-infected cells removed at rate (d.0)
	'parms$delta * Ts'  # active infected cells removed at rate by bursting (delta)
)
names(deaths) <- demes


# susceptible cells grow at constant rate (lambda), die at rate (d.T) 
#  and are depleted by infection
nonDemeDynamics <- c('parms$lambda - parms$d.T * T - parms$k * V * T')
names(nonDemeDynamics) <- c('T')


## default model parameters, taken from Rong and Perelson (2009; PLOS Comput Biol e1000533)
params <- list()
params["lambda"] <- 1e4  # growth rate of uninfected cells (per mL per day)
params["d.T"] <- 0.01    # death rate of uninfected cells (per day)
params["k"] <- 2.4e-8    # rate of infection (mL/day)
params["eta"] <- 0.01    # probability of entering latent state
params["d.0"] <- 0.001   # death rate of latently-infected cells
params["a.L"] <- 0.05    # rate of transition from latently to productively infected cells
params["delta"] <- 1.0   # death rate of productively infected cells (per day)
params["N"] <- 2000      # number of virions produced by cell death
params["c"] <- 23.       # clearance rate of free virus (per day)


get.steady.state <- function(params) {
	# Equation (2) from paper (epsilon = 0)
	V.0 <- with(params, N * lambda / c * (1 - d.0 / (d.0 + a.L) * eta) - d.T / k)
	T.0 <- with(params, lambda / (d.T + k * V.0))
	L.0 <- with(params, eta * k * V.0 * T.0 / (d.0 + a.L))
	Ts.0 <- with(params, c * V.0 / (N * delta))
	
	c(V=V.0, T=T.0, L=L.0, Ts=Ts.0)
}

