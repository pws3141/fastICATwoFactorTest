#! /usr/bin/env Rscript

cat("start ...\n")
library(fastICA)
library(clusterICA)
library(jvcoords)

n <- 3000
n1 <- 2000

# good values 36, 42, 45, 75
for (seed in 1:1000) {
        set.seed(seed)

        #cat("creating random points ...\n")
        x <- cbind(c(rnorm(n1, 1, 0.4), rnorm(n-n1, -1, 0.4)), 
                   c(rnorm(n1, 0, 0.4), rnorm(n-n1, 0, 0.4)))
        x <- whiten(x)$y

        # find the best direction for fastICA
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
}

cat("done\n")


