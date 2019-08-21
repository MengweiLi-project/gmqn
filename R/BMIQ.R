#' Beta MIxture Quantile dilation (BMIQ)
#'
#' The BMIQ was developed by Andrew Teschendorff in 2013 and this implementation is the based on the version modified by Steve Horvath. I also removed and modified some parameters to simplify the function.
#' @param beta.v A list of DNA methylation level
#' @param design.v A list of probe design('I' or 'II')
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
  rand.idx <- sample(1:length(beta1.v), min(c(nfit, length(beta1.v))), replace = FALSE)
  em1.o <- .blc2(Y = matrix(beta1.v[rand.idx], ncol=1),
                w = w0.m[rand.idx, ], maxiter = niter, tol = tol)
  subsetclass1.v <- apply(em1.o$w, 1, which.max)
  subsetth1.v <- c(mean(max(beta1.v[rand.idx[subsetclass1.v == 1]]),
                        min(beta1.v[rand.idx[subsetclass1.v == 2]])),
                   mean(max(beta1.v[rand.idx[subsetclass1.v==2]]),
                        min(beta1.v[rand.idx[subsetclass1.v==3]])))
  class1.v <- rep(2, length(beta1.v))
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3
  nth1.v <- subsetth1.v

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

  rand.idx <- sample(1:length(beta2.v), min(c(nfit, length(beta2.v))), replace = FALSE)
  em2.o <- .blc2(Y = matrix(beta2.v[rand.idx], ncol=1),
                w = w0.m[rand.idx, ], maxiter = niter, tol = tol)

  ### for type II probes assign to state (unmethylated, hemi or full methylation)
  subsetclass2.v <- apply(em2.o$w, 1, which.max)

  if (sum(subsetclass2.v == 2) > 0){
    subsetth2.v <- c(mean(max(beta2.v[rand.idx[subsetclass2.v == 1]]),
                          min(beta2.v[rand.idx[subsetclass2.v == 2]])),
                    mean(max(beta2.v[rand.idx[subsetclass2.v == 2]]),
                         min(beta2.v[rand.idx[subsetclass2.v == 3]])))
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
               em2.o$b[lt, 1], lower.tail = FALSE)
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
