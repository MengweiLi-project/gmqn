#' Gaussian Mixture Quantile Normalization (GMQN)
#'
#' The implementation of Gaussian Mixture Quantile Normalization (GMQN)
#' @param m The vector of methylation signal intensities.
#' @param um The vector of umethylation signal intensities.
#' @param probe the probe id for m and um.
#' @param type The type of methylation array which can be '450k'(default) or '850k'
#' @param ref The reference used to normalize data. By default, it uses the reference set up from GSE105018 which is also used in EWAS Datahub.
#' @return A data frame contains normalized m and um, p value, and DNA methylation level
gmqn_normalize_gamma <- function(m, um, probe, type = '450k', ref) {

  if (type == '450k') {
    t1.red <- row.names(annon_450k)[which(annon_450k$color == 'Red' & annon_450k$probe_type == 'I')]
    t1.green <- row.names(annon_450k)[which(annon_450k$color == 'Grn' & annon_450k$probe_type == 'I')]
    t2 <- row.names(annon_450k)[which(annon_450k$probe_type == 'II')]
    type <- annon_450k[probe, 'probe_type']
  } else if (type == '850k') {
    t1.red <- row.names(annon_850k)[which(annon_850k$color == 'Red' & annon_850k$probe_type == 'I')]
    t1.green <- row.names(annon_850k)[which(annon_850k$color == 'Grn' & annon_850k$probe_type == 'I')]
    t2 <- row.names(annon_850k)[which(annon_850k$probe_type == 'II')]
    type <- annon_850k[probe, 'probe_type']
  }

  set.seed(1)
  print("Fitting Gamma mixture model for probe of type 1 red----------")
  t1.red.index <- match(t1.red, probe)
  t1.red.index <- t1.red.index[which(!is.na(t1.red.index))]
  t1.red.signal <- c(m[t1.red.index], um[t1.red.index])
  t1.red.model <- gammamixEM(t1.red.signal[which(t1.red.signal > 0)][1:10000], lambda = c(0.5, 0.5))
  t1.red.alpha <- as.numeric(t1.red.model$gamma.pars[1,])
  t1.red.beta <- as.numeric(t1.red.model$gamma.pars[2,])

  print("Normalizing probe of type 1 red---------------------------------")
  sub.pro = cbind(pgamma(t1.red.signal, t1.red.alpha[1], t1.red.beta[1]),
                  pgamma(t1.red.signal, t1.red.alpha[2], t1.red.beta[2]))
  subsetclass <- apply(sub.pro, 1, which.max)
  temp = t1.red.signal
  t1.red.signal[which(subsetclass == 1)] <- qgamma(pgamma(t1.red.signal[which(subsetclass == 1)], t1.red.alpha[1], 1/t1.red.beta[1]), ref$t1.red.ref.alpha[1], 1/ref$t1.red.ref.beta[1])

  t1.red.signal[which(subsetclass == 2)] <- qgamma(pgamma(t1.red.signal[which(subsetclass == 2)], t1.red.alpha[2], 1/t1.red.beta[2]), ref$t1.red.ref.alpha[2], 1/ref$t1.red.ref.beta[2])

  t1.red.signal[which(t1.red.signal == Inf)] = temp[which(t1.red.signal == Inf)]

  m[t1.red.index] <- t1.red.signal[1:(length(t1.red.signal)/2)]
  um[t1.red.index] <- t1.red.signal[(length(t1.red.signal)/2+1):length(t1.red.signal)]

  print("Fitting Gamma mixture model for probe of type 1 green--------")
  t1.green.index <- match(t1.green, probe)
  t1.green.index <- t1.green.index[which(!is.na(t1.green.index))]
  t1.green.signal <- c(m[t1.green.index], um[t1.green.index])
  t1.green.model <- gammamixEM(t1.green.signal[which(t1.green.signal > 0)][1:10000], lambda = c(0.5, 0.5))
  t1.green.alpha <- as.numeric(t1.green.model$gamma.pars[1,])
  t1.green.beta <- as.numeric(t1.green.model$gamma.pars[2,])

  print("Normalizing probe of type 1 green-------------------------------")
  sub.pro = cbind(pgamma(t1.green.signal, t1.green.alpha[1], t1.green.beta[1]),
                  pgamma(t1.green.signal, t1.green.alpha[2], t1.green.beta[2]))
  subsetclass <- apply(sub.pro, 1, which.max)
  temp = t1.green.signal
  t1.green.signal[which(subsetclass == 1)] <- qgamma(pgamma(t1.green.signal[which(subsetclass == 1)], t1.green.alpha[1], 1/t1.green.beta[1]), ref$t1.green.ref.alpha[1], 1/ref$t1.green.ref.beta[1])

  t1.green.signal[which(subsetclass == 2)] <- qgamma(pgamma(t1.green.signal[which(subsetclass == 2)], t1.green.alpha[2], 1/t1.green.beta[2]), ref$t1.green.ref.alpha[2], 1/ref$t1.green.ref.beta[2])

  t1.green.signal[which(t1.green.signal == Inf)] = temp[which(t1.green.signal == Inf)]

  m[t1.green.index] <- t1.green.signal[1:(length(t1.green.signal)/2)]
  um[t1.green.index] <- t1.green.signal[(length(t1.green.signal)/2+1):length(t1.green.signal)]

  print("Detecting p value-----------------------------------------------")
  t2.index <- match(t2, probe)
  t2.index <- t2.index[which(!is.na(t2.index))]
  pIR <- apply(cbind(1 - pgamma(m[t1.red.index], t1.red.alpha[1], 1/t1.red.beta[1]),
                     1 - pgamma(um[t1.red.index], t1.red.alpha[1], 1/t1.red.beta[1])),
               1, min)
  pIG <- apply(cbind(1 - pgamma(m[t1.green.index], t1.green.alpha[1], 1/t1.green.beta[1]),
                     1 - pgamma(um[t1.green.index], t1.green.alpha[1], 1/t1.green.beta[1])),
               1, min)
  pII <- apply(cbind(1 - pgamma(um[t2.index], t1.green.alpha[1], 1/t1.green.beta[1]),
                     1 - pgamma(m[t2.index], t1.red.alpha[1], 1/t1.red.beta[1])),
               1, min)
  p <- cbind(c(pIR, pIG , pII),
             c(probe[t1.red.index], probe[t1.green.index], probe[t2.index]))
  row.names(p) <- c(probe[t1.red.index], probe[t1.green.index], probe[t2.index])
  m[which(m <= 0)] <- min(m[which(m > 0)])
  um[which(um <= 0)] <- min(um[which(um > 0)])
  normalized.signal <- data.frame(cbind(round(as.numeric(m), 0),
                                        round(as.numeric(um), 0),
                                        as.numeric(p[probe, 1])))
  row.names(normalized.signal) <- probe
  names(normalized.signal) <- c("m", "um", "p")

  print("Normalizing between probe type----------------------------------")
  normalized.signal[which(normalized.signal$m == Inf), "m"] <- max(normalized.signal$m[normalized.signal$m != Inf])
  normalized.signal[which(normalized.signal$m == 0), "m"] <- min(normalized.signal$m[normalized.signal$m > 0])
  normalized.signal[which(normalized.signal$um == Inf), "um"] <- max(normalized.signal$um[normalized.signal$um != Inf])
  normalized.signal[which(normalized.signal$um == 0), "um"] <- min(normalized.signal$um[normalized.signal$um > 0])
  beta <- normalized.signal$m / (normalized.signal$m + normalized.signal$um)
  beta.adjust <- BMIQ(beta, type)
  beta.adjust <- beta.adjust
  normalized.signal$beta <- beta.adjust
  return(normalized.signal)
  print("Finished--------------------------------------------------------")
}
