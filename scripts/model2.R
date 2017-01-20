#!/usr/bin/Rscript -f
#PLoS Comput Biol. 2009 Oct;5(10):e1000533. doi: 10.1371/journal.pcbi.1000533. Epub 2009 Oct 16.

library(rcolgem)
#library(phydynR)

args <- commandArgs(trailing=T)

#params.file = NA

params.file <- if (length(args > 0)) {
	args[1]
} else {
	NA
}

#show.demographic.process <- function(demo.model, theta, x0, t0, t1, res = 1000, integrationMethod = "lsoda", ...) {
show.demographic.process <- function(tfgy, ...) {
#    tfgy <- demo.model(theta, x0, t0, t1, res = 1000, integrationMethod = integrationMethod)
    o <- tfgy[[5]]
    t <- o[, 1]
    matplot(t, log10(o[, 2:ncol(o)]), type = "l", xlab = "Time", ylab = "", ...)
    legend("bottomright", inset = 0.05, legend = colnames(o)[2:ncol(o)], pch = 1, horiz = TRUE, ...)
    tfgy
}

simulate.binary.dated.tree.fgy <- function (times, births, migrations, demeSizes, sampleTimes, 
    sampleStates, integrationMethod = "rk4", n.reps = 1, cluster = NULL) {
    s <- round(coef(lm(x ~ y, data.frame(x = 1:length(times), 
        y = sort(times))))[1], digits = 9)
    if (s != 1) 
        warning("Tree simulator assumes times given in equal increments")
    n <- length(sampleTimes)
    m <- ncol(sampleStates)
    if (m == 1) 
        return(simulate.bdt.onedeme(times, births, demeSizes, 
            sampleTimes, sampleStates, integrationMethod, n.reps, 
            cluster))
    maxSampleTime <- max(sampleTimes)
    if (length(names(sampleTimes)) == 0) {
        sampleNames <- as.character(1:length(sampleTimes))
        names(sampleTimes) <- sampleNames
        rownames(sampleStates) <- sampleNames
    }
    if (length(rownames(sampleStates)) == 0) 
        warning("simulate.binaryDatedTree.fgy: sampleStates should have row names")
    sampleHeights <- maxSampleTime - sampleTimes
    ix <- sort(sampleHeights, index.return = TRUE)$ix
    sortedSampleHeights <- sampleHeights[ix]
    sortedSampleStates <- as.matrix(sampleStates[ix, ])
    uniqueSortedSampleHeights <- unique(sortedSampleHeights)
    maxtime <- max(times)
    mintime <- min(times)
    maxHeight <- maxSampleTime - mintime
    times_ix <- sort(times, index.return = TRUE, decreasing = TRUE)$ix
    fgyParms <- list()
    fgyParms$FGY_RESOLUTION <- length(times)
    fgyParms$F_DISCRETE <- lapply(times_ix, function(k) births[[k]])
    fgyParms$G_DISCRETE <- lapply(times_ix, function(k) migrations[[k]])
    fgyParms$Y_DISCRETE <- lapply(times_ix, function(k) demeSizes[[k]])
    fgyParms$hoffset = hoffset <- maxtime - maxSampleTime
    if (hoffset < 0) 
        stop("Time axis does not cover the last sample time")
    fgyParms$get.index <- function(h) min(1 + fgyParms$FGY_RESOLUTION * 
        (h + hoffset)/(maxtime - mintime), fgyParms$FGY_RESOLUTION)
    fgyParms$F. <- function(h) {
        fgyParms$F_DISCRETE[[fgyParms$get.index(h)]]
    }
    fgyParms$G. <- function(h) {
        fgyParms$G_DISCRETE[[fgyParms$get.index(h)]]
    }
    fgyParms$Y. <- function(h) {
        fgyParms$Y_DISCRETE[[fgyParms$get.index(h)]]
    }
    fgyParms$FGY_H_BOUNDARIES <- sort(maxSampleTime - times)
    F. <- fgyParms$F.
    G. <- fgyParms$G.
    Y. <- fgyParms$Y.
    FGY_RESOLUTION <- fgyParms$FGY_RESOLUTION
    get.fgy <- function(h) {
        .Y <- fgyParms$Y.(h)
        .F <- fgyParms$F.(h)
        .G <- fgyParms$G.(h)
        list(.F = .F, .G = .G, .Y = .Y)
    }
    heights <- sort(fgyParms$FGY_H_BOUNDARIES)
    heights <- heights[heights <= maxHeight]
    heights <- heights[heights >= 0]
    fgyParms$heights <- heights
    fgymat <- t(sapply(fgyParms$heights, function(h) with(get.fgy(h), 
        {
            c(as.vector(.F), as.vector(.G), as.vector(.Y))
        })))
    fgymat <- pmax(fgymat, 0)
    .solve.Q.A.L <- function(h0, h1, A0, L0) {
        Q0 <- diag(m)
        parameters <- c(m, maxHeight, length(fgyParms$heights), 
            sum(A0), as.vector(fgymat))
        y0 <- c(as.vector(Q0), A0, L0)
        tryCatch({
            o <- ode(y = y0, c(h0, h1), func = "dQAL", parms = parameters, 
                dllname = "rcolgem", initfunc = "initfunc", method = integrationMethod)
        }, error = function(e) {
            return(list())
        })
        Q1 <- t(matrix(abs(o[nrow(o), 2:(1 + m^2)]), nrow = m))
        A1 <- o[nrow(o), (1 + m^2 + 1):(1 + m^2 + m)]
        L1 <- o[nrow(o), ncol(o)]
        return(list(unname(Q1), unname(A1), unname(L1)))
    }
    run1 <- function(repl) {
    	cat(paste0(repl, ", "))
    
        cumSortedSampleStates <- sapply(1:m, function(k) cumsum(sortedSampleStates[, 
            k]))
        cumSortedNotSampledStates <- t(cumSortedSampleStates[n, 
            ] - t(cumSortedSampleStates))
        nsy.index <- approxfun(sort(jitter(sortedSampleHeights, 
            factor = max(1e-06, sortedSampleHeights[length(sortedSampleHeights)]/1e+06))), 
            1:n, method = "constant", rule = 2)
        not.sampled.yet <- function(h) {
            cumSortedNotSampledStates[nsy.index(h), ]
        }
        dA <- function(h, A, parms, ...) {                
            nsy <- not.sampled.yet(h)
            with(get.fgy(h), {
                A_Y <- (A - nsy)/.Y
                A_Y[is.nan(A_Y)] <- 0
                csFpG <- colSums(.F + .G)
                list(setNames(as.vector(.G %*% A_Y - csFpG * 
                  A_Y + (.F %*% A_Y) * pmax(1 - A_Y, 0)), names(A)))
            })
        }
        AIntervals <- c(sortedSampleHeights[2:length(uniqueSortedSampleHeights)], 
            maxHeight)
        h0 <- 0
        sampled.at.h <- function(h) which(sortedSampleHeights == 
            h)
        haxis <- seq(0, maxHeight, length.out = fgyParms$FGY_RESOLUTION)
        odA <- ode(y = colSums(sortedSampleStates), times = haxis, 
            func = dA, parms = NA, method = integrationMethod)
        haxis <- odA[, 1]
        AplusNotSampled <- odA[, 2:(m + 1)]
        AplusNotSampled <- as.matrix(AplusNotSampled, nrows = length(haxis))
        Amono <- rowSums(AplusNotSampled)
        Amono[is.na(Amono)] <- min(Amono[!is.na(Amono)])
        Amono <- Amono - min(Amono)
        Amono <- (max(Amono) - Amono)/max(Amono)
        nodeHeights <- sort(approx(Amono, haxis, xout = runif(n - 
            1, 0, 1))$y)
        uniqueSampleHeights <- unique(sampleHeights)
        eventTimes <- c(uniqueSampleHeights, nodeHeights)
        isSampleEvent <- c(rep(TRUE, length(uniqueSampleHeights)), 
            rep(FALSE, length(nodeHeights)))
        ix <- sort(eventTimes, index.return = TRUE)$ix
        eventTimes <- eventTimes[ix]
        isSampleEvent <- isSampleEvent[ix]
        get.A.index <- function(h) {
            min(1 + floor(fgyParms$FGY_RESOLUTION * (h)/(maxHeight)), 
                fgyParms$FGY_RESOLUTION)
        }
        get.A <- function(h) {
            i <- get.A.index(h)
            AplusNotSampled[i, ] - not.sampled.yet(h)
        }
        S <- 1
        L <- 0
        Nnode <- n - 1
        edge.length <- rep(-1, Nnode + n - 1)
        edge <- matrix(-1, (Nnode + n - 1), 2)
        if (length(names(sortedSampleHeights)) == 0) {
            tip.label <- as.character(1:n)
        }
        else {
            tip.label <- names(sortedSampleHeights)
        }
        maxSampleTime <- max(sampleTimes)
        heights <- rep(0, (Nnode + n))
        parentheights <- rep(-1, (Nnode + n))
        heights[1:n] <- sortedSampleHeights
        inEdgeMap <- rep(-1, Nnode + n)
        outEdgeMap <- matrix(-1, (Nnode + n), 2)
        parent <- 1:(Nnode + n)
        daughters <- matrix(-1, (Nnode + n), 2)
        lstates <- matrix(-1, (Nnode + n), m)
        mstates <- matrix(-1, (Nnode + n), m)
        ustates <- matrix(-1, (Nnode + n), m)
        ssm <- matrix(0, nrow = n, ncol = m)
        lstates[1:n, ] <- sortedSampleStates
        mstates[1:n, ] <- lstates[1:n, ]
        isExtant <- rep(FALSE, Nnode + n)
        isExtant[sampled.at.h(h0)] <- TRUE
        extantLines <- which(isExtant)
        nExtant <- sum(isExtant)
        if (length(extantLines) > 1) {
            A0 <- colSums(as.matrix(sortedSampleStates[extantLines, 
                ], nrow = length(extantLines)))
        }
        else {
            A0 <- sortedSampleStates[extantLines, ]
        }
        lineageCounter <- n + 1
        for (ih in 1:(length(eventTimes) - 1)) {
            h0 <- eventTimes[ih]
            h1 <- eventTimes[ih + 1]
            fgy <- get.fgy(h1)
            nExtant <- sum(isExtant)
            A0 <- get.A(h0)
            out <- .solve.Q.A.L(h0, h1, A0, L)
            Q <- out[[1]]
            A <- out[[2]]
            L <- out[[3]]
            
            if (is.nan(L)) {
                L <- Inf
            }
            if (sum(is.nan(Q)) > 0) 
                Q <- diag(length(A))
            if (sum(is.nan(A)) > 0) 
                A <- A0
            if (nExtant > 1) {
                mstates[isExtant, ] <- t(t(Q) %*% t(mstates[isExtant, 
                  ]))
                mstates[isExtant, ] <- abs(mstates[isExtant, 
                  ])/rowSums(as.matrix(abs(mstates[isExtant, 
                  ]), nrow = length(isExtant)))
                A <- colSums(as.matrix(mstates[isExtant, ], nrow = length(isExtant)))
            }
            else {
                mstates[isExtant, ] <- t(t(Q) %*% mstates[isExtant, 
                  ])
                mstates[isExtant, ] <- abs(mstates[isExtant, 
                  ])/sum(abs(mstates[isExtant, ]))
                A <- mstates[isExtant, ]
            }
            if (isSampleEvent[ih + 1]) {
                sat_h1 <- sampled.at.h(h1)
                isExtant[sat_h1] <- TRUE
                heights[sat_h1] <- h1
            }
            else {
                .F <- fgy$.F
                .G <- fgy$.G
                .Y <- fgy$.Y
                if (nExtant > 1 && is.element(0, .Y)) {
                  return(NA)
                }
                a <- A/.Y
                extantLines <- which(isExtant)
                tryCatch({
                  .lambdamat <- (t(t(a)) %*% a) * .F
                  kl <- sample.int(m^2, size = 1, prob = as.vector(.lambdamat))
                  k <- 1 + ((kl - 1)%%m)
                  l <- 1 + floor((kl - 1)/m)
                  probstates <- as.matrix(mstates[extantLines, 
                    ], nrow = length(extantLines))
                  u_i <- sample.int(nExtant, size = 1, prob = probstates[, 
                    k])
                  probstates[u_i, ] <- 0
                  u <- extantLines[u_i]
                  v <- sample(extantLines, size = 1, prob = probstates[, 
                    l])
                }, error = function(e) browser())
                ustates[u, ] <- mstates[u, ]
                ustates[v, ] <- mstates[v, ]
                a_u <- pmin(1, mstates[u, ]/.Y)
                a_v <- pmin(1, mstates[v, ]/.Y)
                lambda_uv <- ((a_u) %*% t(a_v)) * .F + ((a_v) %*% 
                  t(a_u)) * .F
                palpha <- rowSums(lambda_uv)/sum(lambda_uv)
                alpha <- lineageCounter
                lineageCounter <- lineageCounter + 1
                isExtant[alpha] <- TRUE
                isExtant[u] <- FALSE
                isExtant[v] <- FALSE
                lstates[alpha, ] = mstates[alpha, ] <- palpha
                heights[alpha] <- h1
                uv <- c(u, v)
                inEdgeMap[uv] <- alpha
                outEdgeMap[alpha, ] <- uv
                parent[uv] <- alpha
                parentheights[uv] <- h1
                daughters[alpha, ] <- uv
                edge[u, ] <- c(alpha, u)
                edge.length[u] <- h1 - heights[u]
                edge[v, ] <- c(alpha, v)
                edge.length[v] <- h1 - heights[v]
            }
        }
        tryCatch({
            self <- list(edge = edge, edge.length = edge.length, 
                Nnode = Nnode, tip.label = tip.label, heights = heights, 
                parentheights = parentheights, parent = parent, 
                daughters = daughters, lstates = lstates, mstates = mstates, 
                ustates = ustates, m = m, sampleTimes = sampleTimes, 
                sampleStates = sampleStates, maxSampleTime = maxSampleTime, 
                inEdgeMap = inEdgeMap, outEdgeMap = outEdgeMap)
            class(self) <- c("binaryDatedTree", "phylo")
            
            order.nodes <- c(1:length(self$tip.label), (length(self$tip.label) + self$Nnode):(length(self$tip.label) + 1))
			self$edge <- t(apply(self$edge, 1, function(e) order.nodes[e]))
			self$lstates <- self$lstates[order.nodes,]
			self$ustates <- self$ustates[order.nodes,]
			self$mstates <- self$mstates[order.nodes,]
#            sampleTimes2 <- sampleTimes[names(sortedSampleHeights)]
#            sampleStates2 <- as.matrix(lstates[1:n, ], nrow = n)
#            rownames(sampleStates2) <- tip.label
#            phylo <- read.tree(text = write.tree(self))
#            sampleTimes2 <- sampleTimes2[phylo$tip.label]
#            sampleStates2 <- as.matrix(sampleStates2[phylo$tip.label, 
#                ], nrow = length(phylo$tip.label))
#            bdt <- binaryDatedTree(phylo, sampleTimes2, sampleStates = sampleStates2)
#            bdt$lstates <- lstates
#            bdt$mstates <- mstates
#            bdt$ustates <- ustates
        }, error = function(e) browser())
#        return(bdt)
		return(self)
    }
    if (any(is.null(cluster))) {
        result <- lapply(1:n.reps, run1)
    }
    else {
        result <- parLapply(cluster, 1:n.reps, run1)
    }
    return(result[!is.na(result)])
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
eta <- 0.0001
d.0 <- 0.001
a.L <- 0.1
delta <- 1
N <- 2000
c <- 23
r.L <- 0

V.0 <- 453000
L.0 <- 4.48
Ts.0 <- 5210

T <- 479000
E <- 2.4e-8
tp <- 0

ntimes <- 10
nsamples <- 20
start.time <- 0
end.time <- 1000
treatment.start <- 10000000
treatment.end <- 100000000
epsilon <- 0.7
do.mod <- make.mod(treatment.start, treatment.end, 1-epsilon)

if (!is.null(params.file) && !is.na(params.file)) {
	params <- read.table(params.file, header=F, sep="\t")
	for (i in 1:(length(params$V1))) {
		eval(parse(text=paste(params$V1[i], params$V2[i], sep="=")))
	}
}

demes <- c("V", "L", "Ts")

births <- t(rbind(c('0', '0', 'parms$N * parms$delta * Ts'), c('parms$eta * E * T * V', 'parms$r.L * L', '0'), c('(1 - parms$eta) * E * T * V', '0', '0')))
rownames(births) <- colnames(births) <- demes

migrations <- t(rbind(c('0', '0', '0'), c('0', '0', '0'), c('0', 'parms$a.L * L', '0')))
rownames(migrations) <- colnames(migrations) <- demes

deaths <- c('parms$c * V', 'parms$d.0 * L', 'parms$delta * Ts')
names(deaths) <- demes

nonDemeDynamics <- c('parms$lambda - parms$d.T * T - E * V * T', 'do.mod(t) * E')
names(nonDemeDynamics) <- c('T', 'E')

params.list <- list(lambda=lambda, d.T=d.T, eta=eta, d.0=d.0, a.L=a.L, delta=delta, N=N, c=c, r.L=r.L)
params.names <- names(params.list)
pop.list <- c(V=V.0, L=L.0, Ts=Ts.0, T=T, E=E)

tfgy <- make.fgy(start.time, end.time, births, deaths, nonDemeDynamics, pop.list, migrations=migrations, parms=params.list, fgyResolution=1000, integrationMethod = "lsoda")
#dm <- build.demographic.process(births=births, deaths=deaths, migrations=migrations, nonDemeDynamics=nonDemeDynamics, parameterNames=params.names, rcpp=F, sde=F)

par(mfcol=c(1, 2))

show.demographic.process(tfgy, col=c("red", "green", "blue", "yellow", "cyan"))
#tfgy <- show.demographic.process(dm, theta=params.list, x0=pop.list, t0=0, t1=end.time, col=c("red", "green", "blue", "yellow", "cyan"), res=2000)

cells <- tfgy[[5]][, "L"] + tfgy[[5]][, "Ts"]
if (T) {
sampleStates <- t(matrix(c(1, 0, 0), nrow=3, ncol=ntimes * nsamples * 3))
sampleTimes <- rev(unlist(lapply((1:ntimes)*end.time/ntimes, rep, 3 * nsamples)))

if (T) {
for (i in 0:(ntimes-1)) {
	for (j in 0:(nsamples-1)) {
		sampleStates[i * nsamples * 3 + j * 3 + 1:2,] <- t(rmultinom(n=2, size=1, prob=c(0, tfgy[[5]][1000 / ntimes * (i + 1), "L"]/cells[1000 / ntimes * (i + 1)], tfgy[[5]][1000 / ntimes * (i + 1), "Ts"]/cells[1000 / ntimes * (i + 1)])))
	}
}
}

colnames(sampleStates) <- demes
} else {
sampleTimes <- c(rep(14, 7), rep(1197, 20), rep(1454, 18), rep(1734, 34), rep(2560, 26))
sampleStates <- t(matrix(c(1, 0, 0), nrow=3, ncol=105))
colnames(sampleStates) <- demes

sampleStates[1:7,] <- t(rmultinom(n=7, size=1, prob=c(0, tfgy[[5]][1000 / 2560 * 14, "L"]/cells[1000 / 2560 * 14], tfgy[[5]][1000 / 2560 * 14, "Ts"]/cells[1000 / 2560 * 14])))
sampleStates[8:27,] <- t(rmultinom(n=20, size=1, prob=c(0, tfgy[[5]][1000 / 2560 * 1197, "L"]/cells[1000 / 2560 * 1197], tfgy[[5]][1000 / 2560 * 1197, "Ts"]/cells[1000 / 2560 * 1197])))
sampleStates[28:45,] <- t(rmultinom(n=18, size=1, prob=c(0, tfgy[[5]][1000 / 2560 * 1454, "L"]/cells[1000 / 2560 * 1454], tfgy[[5]][1000 / 2560 * 1454, "Ts"]/cells[1000 / 2560 * 1454])))
sampleStates[46:71,] <- t(rmultinom(n=26, size=1, prob=c(0, tfgy[[5]][1000 / 2560 * 1734, "L"]/cells[1000 / 2560 * 1734], tfgy[[5]][1000 / 2560 * 1734, "Ts"]/cells[1000 / 2560 * 1734])))
sampleStates[80:97,] <- t(rmultinom(n=18, size=1, prob=c(0, tfgy[[5]][1000 / 2560 * 1734, "L"]/cells[1000 / 2560 * 1734], tfgy[[5]][1000 / 2560 * 1734, "Ts"]/cells[1000 / 2560 * 1734])))
}

trees <- simulate.binary.dated.tree.fgy(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod="lsoda", n.reps=100)
tree <- trees[[1]]
#tree <- sim.co.tree.fgy(tfgy, sampleTimes, sampleStates)

tree$tip.label <- apply(tree$sampleStates, 1, function(x) c("V", "L", "Ts")[which(x == 1)])

cols <- unlist(apply(tree$sampleStates, 1, function(x) c("red", "green", "blue")[which(x == 1)]))

edge.cols <- unlist(apply(tree$edge, 1, function(x) {
		co <- tree$lstates[x[2], ]
		rgb(co[1], co[2], co[3])
	}))
	
tree$node.label <- unlist(apply(tree$lstates[(1:tree$Nnode) + length(tree$tip.label),], 1, function(co) rgb(co[1], co[2], co[3])))

plot.phylo(tree, cex=.75, tip.color=cols, edge.color=edge.cols, show.node.label=F)