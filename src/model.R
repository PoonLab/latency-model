library(phydynR)

show.demographic.process <- function (demo.model, theta, x0, t0, t1, res=1000, integrationMethod="lsoda", ...) {
    tfgy <- demo.model(theta, x0, t0, t1, res = res, integrationMethod = integrationMethod)
    o <- tfgy[[5]]
    t <- o[, 1]
    if (((ncol(o) - 1) == 2) & tail(colnames(o), 1) == "V2") {
        plot(t, o[, 2], type = "l", xlab = "Time", ylab = colnames(o)[2], 
            ...)
    }
    else {
        matplot(t, o[, 2:ncol(o)], type = "l", xlab = "Time", 
            ylab = "", ...)
        legend("bottomright", inset = 0.05, legend = colnames(o)[2:ncol(o)], 
            pch = 1, horiz = TRUE, ...)
    }
}

demes <- c("V", "L", "T")

births <- rbind(c('0', '0', 'parms$N * parms$delta * T'), c('parms$f * parms$E * V', 'parms$p.bs * L', '0'), c('(1 - parms$f) * parms$E * V', 'parms$a * L', '0'))
rownames(births) <- demes

migrations <- rbind(c('0', '0', '0'), c('0', '0', '0'), c('0', '0', '0'))
rownames(migrations) <- demes

deaths <- c('parms$c.0  * V', 'parms$delta.l * L -parms$a * L', 'parms$delta * T')

nonDemeDynamics <- c(a='omega * (a - am)')

V <- 7.5e5
L <- 10e5
delta <- 1
c.0 <- 23
f <- 3e-6
N <- 2e4
E <- (1-.133) * 2.248e-6 * 595
T <- c.0 * V / (N * delta)
a <- c.0 * V / (N * L) - (1 - f) * E * V / L
am <- 0
omega <- 0.00939
r <- -0.00171 #a - f * E * V / L
delta.l <- delta
p.bs <- r + 1
params.list <- list(delta=delta, N=N, f=f, E=E, c.0=c.0, delta.l=delta.l, p.bs=p.bs, omega=omega, am=am)
pop.list <- c(V=V, L=L, T=T, a=a)

sampleStates <- matrix(1/3, nrow=40, ncol=3)
rownames(sampleStates) <- rep(seq(1, 4), 10) * seq(0, .9, length.out=40)

dm <- build.demographic.process(births=births, deaths=deaths, nonDemeDynamics=nonDemeDynamics, parameterNames=c('delta', 'N', 'f', 'E', 'c.0', 'delta.l', 'p.bs', 'omega', 'am'), rcpp=F, sde=T)

show.demographic.process(dm, theta=params.list, x0=pop.list, t0=0, t1=2 * 365)

tree <- sim.co.tree(params.list, dm, x0=pop.list, t0=0, sampleTimes=rep(seq(1, 4), 10) * 365, sampleStates=sampleStates, res=1000)

plot(tree)