balanceHD
=========

Estimation of average treatment effects in high dimensions via
approximate residual balancing, as proposed by Athey et al. (2016). Fork
of Stefan Wager’s original package with some convenience functions +
continued development using alternative nuisance function estimators.

To install this package in R, run the following commands:

    library(devtools)
    install_github("apoorvalal/balanceHD")

This package currently works with three optimizers: `mosek`, `pogs`, and
`quadprog`. Mosek is a commercial interior point solver, pogs is a
first-order optimizer, based on ADMM, while quadprog is a standard `R`
optimization library. In general, we achieved best performance with
mosek, and recommend trying optimizers in the order listed above. We
found pogs to be somewhat slower than mosek on the problems we tried.
(Note that we offer two solution strategies based on pogs: `pogs` and
`pogs.dual`. We usually recommend the former, except when p is much
larger than n.) Finally, quadprog performors well on small problems, but
can be much slower for larger problems.

In terms of availability, the optimizer quadprog is easiest to access,
and is available directly from CRAN. Pogs needs to be installed
separately, but is still free. To install pogs, simply install one of
the pre-compiled binaries available from the project
[repository](https://github.com/foges/pogs); see [this
page](https://github.com/foges/pogs/blob/master/src/interface_r/README.md)
for furhter instructions. Finally, mosek is a commercial solver;
however, academic
[licenses](https://www.mosek.com/resources/academic-license) are
available for free. One mosek has been installed, we call into it using
the R package Rmosek.

Example usage:

    set.seed(42)
    library(balanceHD)

    n = 400; p = 1000
    tau = 7
    nclust = 10
    beta = 2 / (1:p) / sqrt(sum(1 / (1:p)^2))
    # cluster assingment
    clust.ptreat = rep(c(0.1, 0.9), nclust / 2)
    cluster.center = 0.5 * matrix(rnorm(nclust * p), nclust, p)
    cluster = sample.int(nclust, n, replace = TRUE)
    X = cluster.center[cluster, ] + matrix(rnorm(n * p), n, p)

    W = rbinom(n, 1, clust.ptreat[cluster])
    Y = X %*% beta + rnorm(n, 0, 1) + tau * W

    # bad
    naive.ate(X, Y, W)

    ## [1] 6.339

main function
-------------

    library(tictoc)

    ## 
    ## Attaching package: 'tictoc'

    ## The following object is masked from 'package:data.table':
    ## 
    ##     shift

    tic()
    tau.hat = residualBalance.ate(X, Y, W, estimate.se = TRUE, optimizer = "mosek")
    toc()

    ## 1.355 sec elapsed

    print(paste("true tau:", tau))

    ## [1] "true tau: 7"

    print(paste("point estimate:", round(tau.hat[1], 2)))

    ## [1] "point estimate: 6.84"

    print(paste0(
      "95% CI for tau: (", round(tau.hat[1] - 1.96 * tau.hat[2], 2),
      ", ", round(tau.hat[1] + 1.96 * tau.hat[2], 2), ")"
    ))

    ## [1] "95% CI for tau: (6.49, 7.2)"

alternate optimizers
--------------------

    tic()
    tau.hat = residualBalance.ate(X, Y, W, estimate.se = TRUE, optimizer = "quadprog")

    ## Warning in approx.balance(XW, balance.target, zeta = zeta,
    ## allow.negative.weights = allow.negative.weights, : bound.gamma = TRUE not
    ## implemented for this optimizer

    ## Warning in approx.balance(XW, balance.target, zeta = zeta,
    ## allow.negative.weights = allow.negative.weights, : bound.gamma = TRUE not
    ## implemented for this optimizer

    toc()

    ## 0.499 sec elapsed

    print(paste("true tau:", tau))

    ## [1] "true tau: 7"

    print(paste("point estimate:", round(tau.hat[1], 2)))

    ## [1] "point estimate: 6.85"

    print(paste0(
      "95% CI for tau: (", round(tau.hat[1] - 1.96 * tau.hat[2], 2),
      ", ", round(tau.hat[1] + 1.96 * tau.hat[2], 2), ")"
    ))

    ## [1] "95% CI for tau: (6.5, 7.2)"

POGS is good for big problems but takes some setup on linux.

    tic()
    tau.hat = residualBalance.ate(X, Y, W, estimate.se = TRUE, optimizer = "pogs")
    toc()

    ## 0.517 sec elapsed

    print(paste("true tau:", tau))

    ## [1] "true tau: 7"

    print(paste("point estimate:", round(tau.hat[1], 2)))

    ## [1] "point estimate: 6.85"

    print(paste0(
      "95% CI for tau: (", round(tau.hat[1] - 1.96 * tau.hat[2], 2),
      ", ", round(tau.hat[1] + 1.96 * tau.hat[2], 2), ")"
    ))

    ## [1] "95% CI for tau: (6.5, 7.2)"

Standard methods
----------------

### AIPW

Default elastic net for both outcome model and pscore.

    ipw.ate(X, Y, W, estimate.se = TRUE)

    ## [1] 6.8227 0.1432

use grf for pscore

    ipw.ate(X, Y, W, prop.method = "randomforest", estimate.se = TRUE)

    ## [1] 6.9013 0.1283

### IPW

    ipw.ate(X, Y, W, fit.method = "none", estimate.se = TRUE)

    ## [1] 6.3583 0.2557

### OM

    elnet.ate(X, Y, W, estimate.se = TRUE)

    ## [1] 6.895 0.134

### double lasso

    twostep.lasso.ate(X, Y, W, estimate.se = TRUE)

    ##  [1]   1   2   3   4   5   9  10  16 108 242 284 319 389 509 742 767 808 865 925
    ## [20]  42  82 106 131 149 160 187 231 250 262 278 286 309 316 353 357 359 375 383
    ## [39] 394 395 404 428 456 458 493 496 500 542 545 572 665 705 729 734 739 785 799
    ## [58] 807 814 818 820 892 895 904 911 937 948 958 961 981
    ##  [1]   1   2   3   4   5   9  10  16 108 242 284 319 389 509 742 767 808 865 925
    ## [20]  42  82 106 131 149 160 187 231 250 262 278 286 309 316 353 357 359 375 383
    ## [39] 394 395 404 428 456 458 493 496 500 542 545 572 665 705 729 734 739 785 799
    ## [58] 807 814 818 820 892 895 904 911 937 948 958 961 981

    ##      1        
    ## 6.6556 0.1765

#### References

Susan Athey, Guido Imbens, and Stefan Wager. <b>Approximate Residual
Balancing: De-Biased Inference of Average Treatment Effects in High
Dimensions.</b> 2016.
\[<a href="http://arxiv.org/pdf/1604.07125.pdf">arxiv</a>\]
