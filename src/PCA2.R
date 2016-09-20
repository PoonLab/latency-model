library(kernlab)

folder <- "../data/pca.extreme/"

cat("initializing...\n")

params.table <- read.table(paste0(folder, "params.txt"), col.names=c("names", "init", "less", "greater"))

params.names <- unlist(params.table$names)
params.list <- params.table$init
names(params.list) <- params.names
params.list <- lapply(params.list, function(x) x)

params.list.less <- params.table$less
names(params.list.less) <- params.names
params.list.less <- lapply(params.list.less, function(x) x)

params.list.greater <- params.table$greater
names(params.list.greater) <- params.names
params.list.greater <- lapply(params.list.greater, function(x) x)

n.reps <- 10
n <- length(unlist(params.list))
n.trials <- 1 + 2 * n

var.tab <- as.data.frame(t(matrix(unlist(params.list), nrow=n, ncol=n.trials)))
names(var.tab) <- params.names

var.tab <- as.data.frame(t(matrix(unlist(params.list), nrow=n, ncol=n.trials)))
names(var.tab) <- params.names

if (T) {
for (i in 2:n.trials) {
	if (i <= n + 1)
		var.tab[i, i - 1] <- params.list.less[i - 1]
	else if (i <= 2 * n + 1)
		var.tab[i, i - n -  1] <- params.list.greater[i - n - 1]
}
}

if (F) {
for (i in 2:n.trials) {
#	if (i <= n + 1)
#		var.tab[i, i-1] <- var.tab[i, i - 1] * 1.1
#	else if (i <= 2 * n + 1)
#		var.tab[i, i - n -  1] <- var.tab[i, i - n - 1] * 0.9
	if (i <= n + 1)
		var.tab[i, i - 1] <- var.tab[i, i - 1] * 1.25
	else
		var.tab[i, i - n - 1] <- var.tab[i, i - n - 1] * 0.75
}
}

if (F) {
n <- length(unlist(params.list))
n.trials <- 1 + 2 * n + 4 * n * (n - 1) / 2

var.tab <- as.data.frame(t(matrix(unlist(params.list), nrow=n, ncol=n.trials)))
names(var.tab) <- params.names

combs <- combn(seq(1, n), 2)

for (i in 2:n.trials) {
	if (i %% 1000 == 0)
		cat(paste0(i, " "))
		
	if (i <= n + 1)
		var.tab[i, i - 1] <- var.tab[i, i - 1] * 1.1
	else if (i <= 2*n + 1)
		var.tab[i, i - n -  1] <- var.tab[i, i - n - 1] * 0.9
	else {
		r <- i - 2 * n - 2
	
		s.j <- r %/% (n * (n - 1))
		s.k <- r %% (n * (n - 1)) %/% (n * (n-1) / 2)
		v <- combs[, r %% (n * (n - 1) / 2) + 1]
				
		var.tab[i, v[1]] <- var.tab[i, v[1]] * (1.1 - .2 * s.j)
		var.tab[i, v[2]] <- var.tab[i, v[2]] * (1.1 - .2 * s.k)
	}
}
}

cat("reading kernels...\n")

kernels <- as.matrix(read.csv(paste0(folder, "kernels.csv"), header=T))

# convert into symmetric matrix
diag(kernels) <- diag(kernels) / 2
kernels <- kernels + t(kernels)

# normalize matrix
norm.km <- matrix(0, nrow=nrow(kernels), ncol=ncol(kernels))
for (i in 1:nrow(kernels)) {
	for (j in 1:ncol(kernels)) {
		norm.km[i,j] <- kernels[i,j] / sqrt(kernels[i,i]*kernels[j,j])
	}
}

kpc <- kpca(as.kernelMatrix(norm.km))

var.tab <- var.tab[rep(seq(1, nrow(var.tab)), rep(n.reps, nrow(var.tab))),]

cols <- function(r) {
	unlist(lapply(r, function(x) c("#ff000080", "#00ff0080", "#0000ff80")[abs(x / c(r[1], r[1] * 1.25, r[1] * 0.75) - 1) < 1e-8]))
}

ex.cols <- c("#ffffff",
	"#ff0000", "#aa5500", "#55aa00",
	"#00ff00", "#00aa55", "#0055aa",
	"#0000ff", "#5500aa", "#aa0055",
	"#800000", "#552b00", "#2b5500",
	"#008000", "#00552b", "#002b55",
	"#000080", "#2b0055", "#55002b")

kernels.u <- as.matrix(read.csv(paste0(folder, "kernels.u.csv"), header=T))

diag(kernels.u) <- diag(kernels.u) / 2

kernels.u <- kernels.u  + t(kernels.u)

kernels.u <- as.kernelMatrix(kernels.u)

kpc.u <- kpca(kernels.u)

cat("doing PCA...\n")
k.pca <- kpca(kernels, 2)

if (F) {
score1 <- unlist(kernels[1,])
score2 <- unlist(kernels[2,])
score2[1] <- score1[2]

score1.u <- unlist(kernels2[1,])
score2.u <- unlist(kernels2[2,])
score2.u[1] <- score1.u[2]

var.tab <- var.tab[rep(seq(1, nrow(var.tab)), rep(n.reps, nrow(var.tab))),]
var.tab <- data.frame(var.tab, kernel1=score1, kernel2=score2, u.kernel1=score1.u, u.kernel2=score2.u)

dudi <- dudi.pca(var.tab, scannf=F, nf=2)
dudi2 <- dudi.pca(var.tab[3:nrow(var.tab),], scannf=F, nf=2)
}