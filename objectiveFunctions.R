#! /usr/bin/env Rscript

cat("start ...\n")
library(fastICA)
library(clusterICA)
library(jvcoords)

# rotate the points anti-clockwise, so that the scatter plot looks nice
rotate.deg = 10

cat("loading points ...\n")
x <- read.csv("points.csv", header = FALSE)
x <- data.matrix(x)
n <- nrow(x)

######################################################################
# create a scatter plot together with the marginals

#cat("Calculating objective function over [-pi, pi] ...\n")
cat("setting up fastICA objective functions ...\n")

# find the best direction for fastICA
fi <- fastICA(x, n.comp = 2, alg.typ = "deflation")
stopifnot(isTRUE(all.equal(cov(fi$S) * (n - 1) / n, diag(2))))
stopifnot(isTRUE(all.equal(x %*% fi$K %*% fi$W, fi$S)))
stopifnot(isTRUE(all.equal(t(fi$W) %*% fi$W, diag(2))))
w <- fi$K %*% fi$W %*% c(1, 0) * sqrt((n - 1) / n)
stopifnot(abs(sum(w^2) - 1) < 1e-6)
xf <- as.vector(x %*% w)
stopifnot(isTRUE(all.equal(xf * sqrt(n / (n - 1)), fi$S[, 1])))
if (w[1] > 0) {
    w <- -w
    xf <- -xf
}
df <- density(xf, bw = "SJ")

# find the best direction for clusterICA
ci <- clusterICA(x)
stopifnot(isTRUE(all.equal(cov(ci$y), diag(2), check.attributes = FALSE)))
v <- fromCoords(ci, c(1, 0), apply.shift = FALSE)
stopifnot(abs(sum(v^2) - 1) < 1e-6)
xc <- as.vector(x %*% v)
if (v[1] < 0) {
    v <- -v
    xc <- -xc
}
dc <- density(xc, bw = "SJ")

##################################
# helper functions

G1 <- function(x) log(cosh(x))
G2 <- function(x) -exp(-x^2/2)
G1.base <- mean(G1(rnorm(1e7)))
G2.base <- -1/sqrt(2)

# systematically try a grid of directions to project onto
M <- 100
theta <- seq(0, pi, length.out = M + 1)[1:M]
W <- rbind(cos(theta), sin(theta))

cat("plotting objective functions ...\n")

entropyGaussian <- 0.5 * (1 + log(2 * pi))

ficaScale <- function(negFica, negEnt) {
                maxFica <- max(negFica)
                maxNegentropy <- max(negEnt)
                scaleConst <- maxNegentropy / maxFica
                res <- scaleConst * negFica
                res
}

pdf("objectivesNegentropy.pdf", 5, 2.4, family = "serif", pointsize = 10)
par(mai = c(.6, .3, .1, .1))

all.proj <- apply(x %*% W, 2, sort)
m <- floor(sqrt(n))

ent <- mSpacingEntropy(x=t(all.proj), m=m)
negEnt <- entropyGaussian - ent
fica1 <- -(colMeans(G1(all.proj)) - G1.base)^2
ficaScaled <- ficaScale(negFica=-fica1, negEnt=negEnt)

plot(theta, negEnt, type = "l", ylim = c(0, 0.45),
     xaxt = "n", 
    yaxt = "n", 
    mgp = c(2.5, 0, 0),
     xlab = expression(theta ~~ "(projection angle)"), ylab = NA)
title(ylab = "non-Gaussianity", mgp = c(0.7, 0, 0))
axis(1, at = (0:4)*pi/4, padj = c(-.8, -.4, -.4, -.4, -1.2), cex.axis = .8,
     labels = c(expression(0),
                expression(1/4 ~ pi),
                expression(1/2 ~ pi),
                expression(3/4 ~ pi),
                expression(pi)))
axis(2, at = 0, padj = 0.6, cex.axis = 0.8, labels = expression(0)) 
abline(v = atan2(v[2], v[1]) %% pi)

lines(theta, ficaScaled, lty = 3)
abline(v = atan2(w[2], w[1]) %% pi, lty = 3)

invisible(dev.off())

cat("done\n")
