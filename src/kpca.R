require(kernlab)

# levels are 0.01, 0.02, 0.05, 0.1, 0.2

# set up plot region
par(mfrow=c(1,2), mar=c(5,5,2,0))

# this does not use labels
m <- read.csv('~/git/latency-model/data/test/test.kmat.csv', header=F)
km <- as.kernelMatrix(as.matrix(m))
kp <- kpca(km)
eig(kp) / sum(eig(kp))  # 51.5% in first two components

plot(rotated(kp), bg=rep(rainbow(5, v=0.8, alpha=0.5), each=20), cex=1.5, pch=rep(21:25, each=20), xlab='First component', ylab='Second component', cex.lab=1.2)
title('Unlabeled tree kernel', adj=0, line=0.5)

# m2 <- read.csv('~/git/latency-model/src/test.lab-kmat.csv', header=F)
# km.2 <- as.kernelMatrix(as.matrix(m2))
# kp.2 <- kpca(km.2)
# plot(rotated(kp.2), col=rep(rainbow(5, v=0.8), each=20), pch=20, cex=2)


m3 <- read.csv('~/git/latency-model/data/test/test.kmat0.csv', header=F)
km.3 <- as.kernelMatrix(as.matrix(m3))
kp.3 <- kpca(km.3)

par(mar=c(5,0,2,5))
plot(rotated(kp.3), bg=rep(rainbow(5, v=0.8, alpha=0.5), each=20), cex=1.5, pch=rep(21:25, each=20), xlab='First component', ylab=NA, yaxt='n', cex.lab=1.2)
title('Labeled tree kernel', adj=0, line=0.5)
par(xpd=NA); legend(x=1.5, y=0.5, legend=c(0.01, 0.02, 0.05, 0.1, 0.2), pch=21:25, pt.cex=1.5, pt.bg=rainbow(5, v=0.8, alpha=0.5), bty='n'); par(xpd=FALSE)

eig(kp.3) / sum(eig(kp.3))  # 57.6% in first two components


# evaluate classifier performance by cross-validation
rsqr <- function(kmat) {
	labels <- rep(c(0.01, 0.02, 0.05, 0.1, 0.2), each=20)
	res <- {}
	for (rep in 1:100) {
		#train <- sample(1:nrow(kmat), nrow(kmat)/2)
		train <- as.vector(sapply(1:5, function(i) sample(1:20,10)+(i-1)*20))
		fit <- ksvm(kmat[train,train], labels[train], type='nu-svr', kernel='matrix')
		test.km <- as.kernelMatrix(kmat[-train,-train][,SVindex(fit), drop=FALSE])
		pred <- predict(fit, test.km)
		res <- c(res, summary(lm(pred~labels[-train]))$r.squared)
	}
	return(res)
}

r1 <- rsqr(km)
r2 <- rsqr(km.3)

boxplot(r1, r2)


