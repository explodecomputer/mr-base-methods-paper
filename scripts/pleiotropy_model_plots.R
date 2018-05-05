# ivw

set.seed(10)

pdf(file="../images/pleiotropy_models.pdf")
par(mfrow=c(2,2))

slope <- 1
interecept <- 0

bx <- runif(10)
by <- bx + rnorm(10, sd=bx * 0.1)
plot(by ~ bx, 
	xlim=c(0, 1), ylim=c(0, 1),
	xlab='SNP effect on exposure', ylab='SNP effect on outcome',
	pch=16, col="black", cex=1,
	main="IVW\nBalanced horizontal pleiotropy"
)
abline(a=0,b=1)


# egger

set.seed(10)

slope <- 1
interecept <- 0.3

bx <- runif(10)
by <- interecept + bx + rnorm(10, sd=bx * 0.1)
plot(by ~ bx, 
	xlim=c(0, 1), ylim=c(0, 1),
	xlab='', ylab='',
	pch=16, col="black", cex=1,
	main="Egger regression\nDirectional horizontal pleiotropy"
)
abline(lm(by ~ bx))
abline(lm(by ~ bx+0), col="grey")


# median

set.seed(10)

bx1 <- runif(10)
by1 <- bx1 + rnorm(10, sd=bx1 * 0.1)

bx2 <- runif(7)
by2 <- bx2 * 0.3 + rnorm(7, sd=bx2 * 0.1)

bx <- c(bx1, bx2)
by <- c(by1, by2)

plot(by1 ~ bx1, 
	xlim=c(0, 1), ylim=c(0, 1),
	xlab='', ylab='',
	pch=16, col="black", cex=1,
	main="Median-based estimator\nMinority horizontal pleiotropy"
)
points(by2 ~ bx2, col="red", pch=16)
abline(lm(by1 ~ bx1 + 0))
abline(lm(by ~ bx + 0), col="grey")


# mode

set.seed(10)

bx1 <- runif(10)
by1 <- bx1 + rnorm(10, sd=bx1 * 0.1)

bx2 <- runif(8)
by2 <- bx2 * 0.3 + rnorm(8, sd=bx2 * 0.01)

bx3 <- runif(8)
by3 <- bx3 * 3 + rnorm(8, sd=bx3 * 0.01)

bx4 <- runif(8)
by4 <- bx4 * 5 + rnorm(8, sd=bx4 * 0.01)

bx <- c(bx1, bx2, bx3, bx4)
by <- c(by1, by2, by3, by4)

plot(by1 ~ bx1, 
	xlim=c(0, 1), ylim=c(0, 1),
	xlab='', ylab='',
	pch=16, col="black", cex=1,
	main="Mode-based estimator\nMajority horizontal pleiotropy"
)
points(by2 ~ bx2, col="red", pch=16)
points(by3 ~ bx3, col="red", pch=16)
points(by4 ~ bx4, col="red", pch=16)
abline(b=1, a=0)
abline(b=0.3, a=0, col="grey")
abline(b=3, a=0, col="grey")
abline(b=5, a=0, col="grey")

dev.off()
