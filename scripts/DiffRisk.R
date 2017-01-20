library(rcolgem)
library(kernlab)

n.reps <- 10
start.time <- 0
end.time <- 60.*52

N = 3000  # total population size
n.tips <- 1000
p = 0.5  # frequency of risk group 1

# model parameters
beta = 0.01
gamma = 1/520.
mu = 1/3640.
c1 = 2^(-3:3)
c2 = 1.0
rho = 0.9


# initial population frequencies
S1 = p*N-1
S2 = (1-p)*N
I1 = 1
I2 = 0
pop.list <- c(I1=I1, I2=I2, S1=S1, S2=S2)

# define ODE system
demes <- c('I1', 'I2')  # two risk groups

p11 <- '(parms$rho + (1-parms$rho) * parms$c1*(S1+I1) / (parms$c1*(S1+I1) + parms$c2*(S2+I2)))'
p12 <- '(1-parms$rho) * parms$c2*(S2+I2) / (parms$c1*(S1+I1) + parms$c2*(S2+I2))'
p21 <- '(1-parms$rho) * parms$c1*(S1+I1) / (parms$c1*(S1+I1) + parms$c2*(S2+I2))'
p22 <- '(parms$rho + (1-parms$rho) * parms$c2*(S2+I2) / (parms$c1*(S1+I1) + parms$c2*(S2+I2)))'

births <- rbind(c(paste(sep='*', 'parms$beta*parms$c1', p11, 'I1/(S1+I1)*S1'), 
				paste(sep='*', 'parms$beta*parms$c2', p21, 'I1/(S1+I1)*S2')),
				c(paste(sep='*', 'parms$beta*parms$c1', p12, 'I2/(S2+I2)*S1'),
				paste(sep='*', 'parms$beta*parms$c2', p22, 'I2/(S2+I2)*S2')))

rownames(births)=colnames(births) <- demes

migrations <- rbind(c('0', '0'), c('0', '0'))
rownames(migrations)=colnames(migrations) <- demes

deaths <- c('(parms$mu+parms$gamma)*I1', '(parms$mu+parms$gamma)*I2')
names(deaths) <- demes

# dynamics for susceptible classes (S)
nonDemeDynamics <- c(paste(sep='', '-parms$mu*S1 + parms$mu*S1 + (parms$mu+parms$gamma)*I1', paste(sep='*', '-S1*(parms$beta*parms$c1', p11, 'I1/(S1+I1) + parms$beta*parms$c1', p12, 'I2/(S2+I2))')),
	paste(sep='', '-parms$mu*S2 + parms$mu*S2 + (parms$mu+parms$gamma)*I2', paste(sep='*', '-S2*(parms$beta*parms$c2', p21, 'I1/(S1+I1) + parms$beta*parms$c2', p22, 'I2/(S2+I2))')))

names(nonDemeDynamics) <- c('S1', 'S2')


supress <- lapply(0:6, function(i) {
cat(paste0(i, ": making fgy...\n"))
params <- list(beta=beta, gamma=gamma, mu=mu, c1=c1[i+1], c2=c2, rho=rho)

tfgy <- make.fgy(start.time, end.time, births, deaths, nonDemeDynamics, pop.list, parms=params, migrations=migrations, fgyResolution=10000, integrationMethod = "lsoda")

show.demographic.process <- function(tfgy, ...) {
#    tfgy <- demo.model(theta, x0, t0, t1, res = 1000, integrationMethod = integrationMethod)
    o <- tfgy[[5]]
    t <- o[, 1]
    matplot(t, o[, 2:ncol(o)], type = "l", xlab = "Time", ylab = "", ...)
    legend("bottomright", inset = 0.05, legend = colnames(o)[2:ncol(o)], pch = 1, horiz = TRUE, ...)
#    tfgy
}

pdf(sprintf("../data/DiffRisk/plot.%02d.pdf", i), height=8, width=8)

show.demographic.process(tfgy, col=c('red', 'green', 'cyan', 'magenta'))

dev.off()

demes.t.end <- tfgy[[4]][[1]]

demes.sample <- sample(rep(1:length(demes), times=round(demes.t.end)), size=n.tips)

sampleStates <- matrix((1:2 == demes.sample) * 1, nrow=n.tips, ncol=length(demes))
colnames(sampleStates) <- demes
for (x in 1:n.tips) {
	sampleStates[x, demes.sample[x]] <- 1
}
rownames(sampleStates) <- paste(1:n.tips, demes.sample, sep='_')

sampleTimes <- rep(end.time, nrow(sampleStates))
names(sampleTimes) <- NULL

cat("making trees...\n")

trees <- simulate.binary.dated.tree.fgy(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod="lsoda", n.reps=n.reps)

#trees <- lapply(trees, function(tree) {tree$tip.label <- paste(tree$tip.label, apply(tree$sampleStates, 1, function(x) c("I")[which(x == 1)]), sep="."); tree})

lapply(1:length(trees), function(j) write.tree(trees[[j]], file=sprintf("../data/DiffRisk/tree.%04d.tre", i * n.reps + j)))
})

cat("computing kernel...")
system("python prepareKernel.py ../data/DiffRisk/ \"\"")

cat("doing PCA...")
tab <- read.csv("../data/DiffRisk/kernels.u.csv", header=T)
mat.tab <- as.matrix(tab)
diag(mat.tab) <- diag(mat.tab) / 2
mat.tab <- mat.tab + t(mat.tab)

kpc <- kpca(as.kernelMatrix(mat.tab), features=2) 

pdf("../data/DiffRisk/PCA.pdf", height=8, width=8)
cols=unlist(as.list(t(matrix(rep(seq(1:7), 10), nrow=7))))
plot(rotated(kpc), col=cols)
legend("topright", legend=c1, col=1:7, pch=16)
suppress <- dev.off()