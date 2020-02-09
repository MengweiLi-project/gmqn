# We slightly modified swan function in minfi package. When the signal value is less than 0, we make it equal to the non-zero minimum value not the background signal intensity.

.modified_SWAN <- function(Meth, Unmeth, probe) {
  set.seed(1)
  counts <- CpG.counts[CpG.counts$Name %in% probe, ]
  subset <- min(
    table(counts$nCpG[counts$Type == "I" & counts$nCpG %in% 1:3]),
    table(counts$nCpG[counts$Type == "II" & counts$nCpG %in% 1:3]))

  xNormSet <- lapply(c("I", "II"), function(type) {
    .getSubset(counts$nCpG[counts$Type == type], subset)
  })

  M <- .SWAN(
    x = Meth,
    xNormSet = xNormSet,
    counts = counts)
  U <- .SWAN(
    x = Unmeth,
    xNormSet = xNormSet,
    counts = counts)

  return(list(M = M, U = U))
}

.getSubset <- function(counts, subset){
  x <- integer(0)
  for (i in 1:3) {
    x <- c(x, sample(seq.int(1, length(counts))[counts == i], subset))
  }
  seq.int(1, length(counts)) %in% x
}

.SWAN = function(x, xNormSet, counts) {
  normalized_x <- matrix(NA_real_,
                         ncol = ncol(x),
                         nrow = nrow(x),
                         dimnames = dimnames(x))
  typeI_idx <- rownames(x) %in% counts$Name[counts$Type == "I"]
  typeII_idx <- rownames(x) %in% counts$Name[counts$Type == "II"]
  for (j in seq_len(ncol(x))) {
    normalized_x[, j] <- .normaliseChannel(
      intensityI = x[typeI_idx, j],
      intensityII = x[typeII_idx, j],
      xNormSet = xNormSet)
  }
  normalized_x
}

.normaliseChannel <- function(intensityI, intensityII, xNormSet) {
  xTarget <- .aveQuantile(
    list(intensityI[xNormSet[[1]]], intensityII[xNormSet[[2]]]))
  xNorm <- unlist(
    .subsetQuantileNorm(
      list(intensityI, intensityII), xNormSet, xTarget))
  names(xNorm) <- c(names(intensityI), names(intensityII))
  xNorm
}

.aveQuantile <- function(X) {
  nbrOfChannels <- length(X)
  if (nbrOfChannels == 1) {
    return(X)
  }
  nbrOfObservations <- unlist(lapply(X, FUN = length), use.names = FALSE)
  maxNbrOfObservations <- max(nbrOfObservations)
  if (maxNbrOfObservations == 1) {
    return(X)
  }
  ## nbrOfFiniteObservations <- rep(maxNbrOfObservations, times = nbrOfChannels)
  quantiles <- (0:(maxNbrOfObservations - 1))/(maxNbrOfObservations - 1)
  xTarget <- vector("double", maxNbrOfObservations)
  for (cc in 1:nbrOfChannels) {
    Xcc <- X[[cc]]
    Scc <- sort(Xcc)
    nobs <- length(Scc)
    if (nobs < maxNbrOfObservations) {
      ## tt <- !is.na(Xcc)
      bins <- (0:(nobs - 1))/(nobs - 1)
      Scc <- approx(x = bins, y = Scc, xout = quantiles,ties = "ordered")$y
    }
    xTarget <- xTarget + Scc
  }
  rm(Scc, Xcc)
  xTarget <- xTarget/nbrOfChannels
  xTarget
}


.subsetQuantileNorm <- function(x, xNormSet, xTarget) {
  for(i in 1:length(x)){
    n <- length(x[[i]])
    nTarget <- length(xTarget)
    nNormSet <- sum(xNormSet[[i]])

    if(nNormSet != nTarget){
      targetQuantiles <- (0:(nTarget - 1))/(nTarget - 1)
      r <- rank(x[xNormSet[,i], i])
      xNew <-(r - 1)/(nNormSet - 1)
      xNew <- xNew[order(xNew)]
      xNorm <- approx(x = targetQuantiles, y = xTarget, xout = xNew, ties = "ordered", rule = 2)$y
    } else {
      xNorm<-xTarget
    }

    r <- rank(x[[i]])
    xNew <-(r - 1)/(n - 1)
    quantiles <- xNew[xNormSet[[i]]]
    quantiles <- quantiles[order(quantiles)]
    xmin <- min(x[[i]][xNormSet[[i]]]) #get min value from subset
    xmax <- max(x[[i]][xNormSet[[i]]]) #get max value from subset
    kmax <- which(xNew > max(quantiles))
    kmin<- which(xNew < min(quantiles))
    offsets.max <- x[[i]][kmax]-xmax
    offsets.min <- x[[i]][kmin]-xmin
    x[[i]] <- approx(x = quantiles, y = xNorm, xout = xNew, ties = "ordered")$y #interpolate
    x[[i]][kmax]<- max(xNorm) + offsets.max
    x[[i]][kmin]<- min(xNorm) + offsets.min
    #here is what we modified
    x[[i]] = ifelse(x[[i]] <= 0, min(x[[i]][which(x[[i]] > 0)]), x[[i]])
  }
  x
}
