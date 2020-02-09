#' Beta MIxture Quantile dilation (BMIQ)
#'
#' The BMIQ was developed by Andrew Teschendorff in 2013 and this implementation is the based on the version modified by Steve Horvath. I also removed and modified some parameters to simplify the function.
#' @param beta.v A list of DNA methylation level
#' @param design.v A list of probe design('I' or 'II')
#' @param nfit number of probes of a given design type to use for the fitting. Default is 50000. Smaller values (~10000) will make BMIQ run faster at the expense of a small loss in accuracy. For most applications, 5000 or 10000 is ok.
#' @param th1.v thresholds used for the initialisation of the EM-algorithm, they should represent buest guesses for calling probes hemi-methylated and methylated, and will be refined by the EM algorithm. Default values work well in most cases.
#' @param niter maximum number of EM iterations to do. This number should be large enough to yield good fits to the type1 distribution. By default 5.
#' @param tol tolerance convergence threshold for EM algorithm. By default 0.001.
#' @return A data frame contains normalized m and um, p value, and DNA methylation level
BMIQ <- function(beta.v, design.v, nfit = 50000, th1.v = c(0.2,0.75),
                 niter = 5, tol = 0.001){

  type1.idx <- which(design.v == 'I')
  type2.idx <- which(design.v == 'II')
  beta1.v <- beta.v[type1.idx]
  beta2.v <- beta.v[type2.idx]

  ### estimate initial weight matrix from type1 distribution
  w0.m <- matrix(0, nrow = length(beta1.v), ncol = 3)
  w0.m[which(beta1.v < th1.v[1]), 1] <- 1
  w0.m[intersect(which(beta1.v >= th1.v[1]), which(beta1.v <= th1.v[2])), 2] <- 1
  w0.m[which(beta1.v > th1.v[2]), 3] <- 1

  ### fit type1 probe
  set.seed(1)
  hypo.m = which(beta1.v < 0.2)
  hyper.m = which(beta1.v > 0.75)
  median.m = which(beta1.v >= 0.2 & beta1.v <= 0.75)
  rand.idx = c(sample(hypo.m, 10000, replace = T),
               sample(hyper.m, 10000, replace = T),
               sample(median.m, 10000, replace = T))
  #rand.idx <- sample(1:length(beta1.v), min(c(nfit, length(beta1.v))), replace = FALSE)
  em1.o <- .blc2(Y = matrix(beta1.v[rand.idx], ncol=1),
                w = w0.m[rand.idx, ], maxiter = niter, tol = tol)
  subsetclass1.v <- apply(em1.o$w, 1, which.max)
  subsetth1.v <- c(mean(c(max(beta1.v[rand.idx[subsetclass1.v == 1]]),
                        min(beta1.v[rand.idx[subsetclass1.v == 2]]))),
                   mean(c(max(beta1.v[rand.idx[subsetclass1.v==2]]),
                        min(beta1.v[rand.idx[subsetclass1.v==3]]))))
  class1.v <- rep(2, length(beta1.v))
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3
  nth1.v <- subsetth1.v
  print(subsetth1.v)
  ### Estimate Modes
  if (sum(class1.v == 1) == 1){ mod1U <- beta1.v[class1.v == 1] }
  if (sum(class1.v == 3) == 1){ mod1M <- beta1.v[class1.v == 3] }
  if (sum(class1.v == 1) > 1){
    d1U.o <- density(beta1.v[class1.v == 1])
    mod1U <- d1U.o$x[which.max(d1U.o$y)]
  }
  if (sum(class1.v == 3) > 1){
    d1M.o <- density(beta1.v[class1.v == 3])
    mod1M <- d1M.o$x[which.max(d1M.o$y)]
  }

  d2U.o <- density(beta2.v[which(beta2.v < 0.4)])
  d2M.o <- density(beta2.v[which(beta2.v > 0.6)])
  mod2U <- d2U.o$x[which.max(d2U.o$y)]
  mod2M <- d2M.o$x[which.max(d2M.o$y)]

  ### now deal with type2 fit
  th2.v <- vector()
  th2.v[1] <- nth1.v[1] + (mod2U - mod1U)
  th2.v[2] <- nth1.v[2] + (mod2M - mod1M)

  ### estimate initial weight matrix
  w0.m <- matrix(0, nrow = length(beta2.v), ncol = 3)
  w0.m[which(beta2.v <= th2.v[1]), 1] <- 1
  w0.m[intersect(which(beta2.v > th2.v[1]), which(beta2.v <= th2.v[2])), 2] <- 1
  w0.m[which(beta2.v > th2.v[2]), 3] <- 1

  hypo.m = which(beta2.v < th2.v[1])
  hyper.m = which(beta2.v > th2.v[2])
  median.m = which(beta2.v >= th2.v[1] & beta2.v <= th2.v[2])
  rand.idx = c(sample(hypo.m, 10000, replace = T),
               sample(hyper.m, 10000, replace = T),
               sample(median.m, 10000, replace = T))

  #rand.idx <- sample(1:length(beta2.v), min(c(nfit, length(beta2.v))), replace = FALSE)
  em2.o <- .blc2(Y = matrix(beta2.v[rand.idx], ncol=1),
                w = w0.m[rand.idx, ], maxiter = niter, tol = tol)

  ### for type II probes assign to state (unmethylated, hemi or full methylation)
  subsetclass2.v <- apply(em2.o$w, 1, which.max)

  if (sum(subsetclass2.v == 2) > 0){
    subsetth2.v <- c(mean(c(max(beta2.v[rand.idx[subsetclass2.v == 1]]),
                          min(beta2.v[rand.idx[subsetclass2.v == 2]]))),
                    mean(c(max(beta2.v[rand.idx[subsetclass2.v == 2]]),
                         min(beta2.v[rand.idx[subsetclass2.v == 3]]))))
  }
  if (sum(subsetclass2.v == 2) == 0){
    subsetth2.v <- c(1 / 2 * max(beta2.v[rand.idx[subsetclass2.v == 1]]) +
                       1 / 2 * mean(beta2.v[rand.idx[subsetclass2.v == 3]]),
                     1 / 3 * max(beta2.v[rand.idx[subsetclass2.v == 1]]) +
                       2 / 3 * mean(beta2.v[rand.idx[subsetclass2.v == 3]]))
  }

  class2.v <- rep(2, length(beta2.v))
  class2.v[which(beta2.v <= subsetth2.v[1])] <- 1
  class2.v[which(beta2.v >= subsetth2.v[2])] <- 3

  classAV1.v <- vector()
  classAV2.v <- vector()
  for(l in 1:3){
    classAV1.v[l] <- em1.o$mu[l, 1]
    classAV2.v[l] <- em2.o$mu[l, 1]
  }

  ### start normalising input probes
  nbeta2.v <- beta2.v
  ### select U probes
  lt <- 1
  selU.idx <- which(class2.v == lt)
  selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])]
  selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])]
  ### find prob according to typeII distribution
  p.v <- pbeta(beta2.v[selUR.idx], em2.o$a[lt, 1],
               em2.o$b[lt, 1], lower.tail = FALSE) ######
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v, em1.o$a[lt, 1],
               em1.o$b[lt, 1], lower.tail = FALSE)
  nbeta2.v[selUR.idx] <- q.v

  p.v <- pbeta(beta2.v[selUL.idx], em2.o$a[lt, 1],
               em2.o$b[lt, 1], lower.tail=TRUE)
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v, em1.o$a[lt, 1],
               em1.o$b[lt, 1], lower.tail=TRUE)
  nbeta2.v[selUL.idx] <- q.v

  ### select M probes
  lt <- 3
  selM.idx <- which(class2.v == lt)
  selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])]
  selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])]
  ### find prob according to typeII distribution
  p.v <- pbeta(beta2.v[selMR.idx], em2.o$a[lt, 1],
               em2.o$b[lt, 1], lower.tail = FALSE)
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v, em1.o$a[lt, 1],
               em1.o$b[lt, 1], lower.tail = FALSE)
  nbeta2.v[selMR.idx] <- q.v


  lt <- 2
  selH.idx <- c(which(class2.v == lt), selML.idx)
  minH <- min(beta2.v[selH.idx])
  maxH <- max(beta2.v[selH.idx])
  deltaH <- maxH - minH
  deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
  deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])

  nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM
  ## new minimum of H probes should be
  nminH <- max(nbeta2.v[selU.idx]) + deltaUH
  ndeltaH <- nmaxH - nminH

  ### perform conformal transformation (shift+dilation)
  hf <- ndeltaH / deltaH
  nbeta2.v[selH.idx] <- nminH + hf * (beta2.v[selH.idx] - minH)


  pnbeta.v <- beta.v
  pnbeta.v[type1.idx] <- beta1.v
  pnbeta.v[type2.idx] <- nbeta2.v

  return(pnbeta.v)
}

.betaEst2 <- function (y, w, weights) {
  yobs <- !is.na(y)
  if (sum(yobs) <= 1)
    return(c(1, 1))
  y <- y[yobs]
  w <- w[yobs]
  weights <- weights[yobs]
  N <- sum(weights * w)
  p <- sum(weights * w * y) / N
  v <- sum(weights * w * y * y) / N - p * p
  logab = log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p)/v -1))
  if (sum(yobs) == 2)
    return(exp(logab))
  opt <- try(optim(logab, betaObjf, ydata = y, wdata = w, weights = weights,
                   method = "Nelder-Mead", control = list(maxit = 50)), silent = TRUE)
  if (inherits(opt, "try-error"))
    return(c(1, 1))
  exp(opt$par)
}

.blc2 <- function(Y, w, maxiter = 25, tol = 1e-06, weights = NULL, verbose = FALSE) {
  Ymn <- min(Y[Y > 0], na.rm = TRUE)
  Ymx <- max(Y[Y < 1], na.rm = TRUE)
  Y <- pmax(Y, Ymn/2)
  Y <- pmin(Y, 1 - (1 - Ymx) / 2)
  Yobs <- !is.na(Y)
  J <- dim(Y)[2]
  K <- dim(w)[2]
  n <- dim(w)[1]
  if (n != dim(Y)[1])
    stop("Dimensions of w and Y do not agree")
  if (is.null(weights))
    weights = rep(1, n)
  mu = a = b = matrix(Inf, K, J)
  crit <- Inf
  for (i in 1:maxiter) {
    warn0 <- options()$warn
    options(warn = -1)
    eta <- apply(weights * w, 2, sum)/sum(weights)
    mu0 <- mu
    for (k in 1:K) {
      for (j in 1:J) {
        ab <- .betaEst2(Y[, j], w[, k], weights)
        a[k, j] <- ab[1]
        b[k, j] <- ab[2]
        mu[k, j] <- ab[1]/sum(ab)
      }
    }
    ww <- array(0, dim = c(n, J, K))
    for (k in 1:K) {
      for (j in 1:J) {
        ww[Yobs[, j], j, k] <- dbeta(Y[Yobs[, j], j],
                                     a[k, j], b[k, j], log = TRUE)
      }
    }
    options(warn = warn0)
    w <- apply(ww, c(1, 3), sum, na.rm = TRUE)
    wmax <- apply(w, 1, max)
    for (k in 1:K) w[, k] = w[, k] - wmax
    w <- t(eta * t(exp(w)))
    like <- apply(w, 1, sum)
    w <- (1/like) * w
    llike <- weights * (log(like) + wmax)
    crit <- max(abs(mu - mu0))
    if (verbose)
      print(crit)
    if (crit < tol)
      break
  }
  return(list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike)))
}


