#! /usr/bin/env Rscript

cat("start ...\n")
library(fastICA)
library(clusterICA)
library(jvcoords)

n <- 3000
n1 <- 2000

# good values 36, 42, 45, 75
seed <- 1000

set.seed(seed)

#cat("creating random points ...\n")
x <- cbind(c(rnorm(n1, 1, 0.4), rnorm(n-n1, -1, 0.4)), 
           c(rnorm(n1, 0, 0.4), rnorm(n-n1, 0, 0.4)))
x <- whiten(x)$y

pdf("pointsSeed1000Fast100.pdf", 5, 6, 
    family = "serif", pointsize = 10)
    
par(mai = c(0, 0, 0, 0))

plot(x, xlim = c(-6, 4), ylim = c(-6, 6), asp = 1, pch = 20, cex = .6,
     col = rgb(0, 0, 0, .5),
     axes = FALSE, xlab = "", ylab = "")

d.dist <- 4

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

Q <- matrix(c(v[1], -v[2], v[2], v[1]), 2, 2)
lines(c(-4 * v[1], 4 * v[1]), c(-4 * v[2], 4 * v[2]))
for (pos in seq(-pi, pi, l = 7)) {
    a <- c(pos, 0) %*% Q
    b <- c(pos, d.dist) %*% Q
    arrows(a[1], a[2], b[1], b[2], length = .08, lwd = .6)
}
z <- cbind(dc$x, 4 * dc$y + d.dist)
A <- z %*% Q
polygon(A, border = NA, col = "grey")
z <- rbind(c(-4, d.dist), c(4, d.dist))
lines(z %*% Q)
lines(A)

# find the best direction for fastICA
# multiple fastICA runs
cat("Running fastICA 100 times...\n")
for (run in 1:100) {
        fi <- fastICA(x, n.comp = 2, alg.typ = "deflation")
        w <- fi$K %*% fi$W %*% c(1, 0) * sqrt((n - 1) / n)

        # find the best direction for clusterICA
        ci <- clusterICA(x)
        v <- fromCoords(ci, c(1, 0), apply.shift = FALSE)

        cosTheta <- t(w) %*% v
        sinTheta <- sqrt(1 - cosTheta)

        if (sinTheta > 1e-1 & sinTheta < sqrt(2) - 1e-1) {
                cat(seed, sinTheta, "\n")
        } else {
                cat(".")
        }

        # add projection to plot
        Q <- matrix(c(w[1], -w[2], w[2], w[1]), 2, 2)
        lines(c(-4 * w[1], 4 * w[1]), c(-4 * w[2], 4 * w[2]), lty = 3)
}

invisible(dev.off())

cat("done...\n")
