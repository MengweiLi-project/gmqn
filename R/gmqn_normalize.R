#' Gaussian Mixture Quantile Normalization (GMQN)
#'
#' The implementation of Gaussian Mixture Quantile Normalization (GMQN)
#' @param m The vector of methylation signal intensities.
#' @param um The vector of umethylation signal intensities.
#' @param probe the probe id for m and um.
#' @param type The type of methylation array which can be '450k'(default) or '850k'
#' @param ref The reference used to normalize data. By default, it uses the reference set up from GSE105018 which is also used in EWAS Datahub.
#' @return A data frame contains normalized m and um, p value
gmqn_normalize <- function(m, um, probe, type = '450k', ref = 'default', verbose = TRUE) {

  if (ref == 'default') {
    ref = list(t1.green.ref.mean=c(268.9645,6884.7993),
               t1.green.ref.sd=c(203.8351,3682.3554),
               t1.red.ref.mean=c(519.537,8662.520),
               t1.red.ref.sd=c(380.9977,5064.9526))
  }

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
  if(verbose) message("Fitting Gaussian mixture model for probe of type 1 red")
  t1.red.index <- match(t1.red, probe)
  t1.red.index <- t1.red.index[which(!is.na(t1.red.index))]
  t1.red.signal <- c(m[t1.red.index], um[t1.red.index])
  t1.red.model <- Mclust(t1.red.signal, G=2, verbose = F)
  t1.red.mean <- t1.red.model$parameters$mean
  t1.red.sd <- sqrt(t1.red.model$parameters$variance$sigmasq)

  pIR <- apply(cbind(1 - pnorm(m[t1.red.index], t1.red.mean[1], t1.red.sd[1]),
                     1 - pnorm(um[t1.red.index], t1.red.mean[1], t1.red.sd[1])),
               1, min)

  if(verbose) message("Normalizing probe of type 1 red")

  temp = t1.red.signal
  #t1.red.signal[which(t1.red.model$classification == 1)] <- qnorm(pnorm(t1.red.signal[which(t1.red.model$classification == 1)], t1.red.mean[1], t1.red.sd[1]), ref$t1.red.ref.mean[1], ref$t1.red.ref.sd[1])

  t1.red.signal[which(t1.red.model$classification == 1 & pnorm(t1.red.signal, t1.red.mean[1], t1.red.sd[1]) < 0.9 & pnorm(t1.red.signal, t1.red.mean[1], t1.red.sd[1]) > 0.1)] <- qnorm(pnorm(t1.red.signal[which(t1.red.model$classification == 1 & pnorm(t1.red.signal, t1.red.mean[1], t1.red.sd[1]) < 0.9 & pnorm(t1.red.signal, t1.red.mean[1], t1.red.sd[1]) > 0.1)], t1.red.mean[1], t1.red.sd[1]), ref$t1.red.ref.mean[1], ref$t1.red.ref.sd[1])

  # t1.red.signal[which(t1.red.model$classification == 2)] <- qnorm(pnorm(t1.red.signal[which(t1.red.model$classification == 2)], t1.red.mean[2], t1.red.sd[2]), ref$t1.red.ref.mean[2], ref$t1.red.ref.sd[2])
  t1.red.signal[which(t1.red.model$classification == 2 & pnorm(t1.red.signal, t1.red.mean[2], t1.red.sd[2]) > 0.1 & pnorm(t1.red.signal, t1.red.mean[2], t1.red.sd[2]) < 0.9)] <- qnorm(pnorm(t1.red.signal[which(t1.red.model$classification == 2 & pnorm(t1.red.signal, t1.red.mean[2], t1.red.sd[2]) > 0.1 & pnorm(t1.red.signal, t1.red.mean[2], t1.red.sd[2]) < 0.9)], t1.red.mean[2], t1.red.sd[2]), ref$t1.red.ref.mean[2], ref$t1.red.ref.sd[2])


  t1.red.signal[which(t1.red.signal == Inf)] = temp[which(t1.red.signal == Inf)]

  m[t1.red.index] <- t1.red.signal[1:(length(t1.red.signal)/2)]
  um[t1.red.index] <- t1.red.signal[(length(t1.red.signal)/2+1):length(t1.red.signal)]

  if(verbose) message("Fitting Gaussian mixture model for probe of type 1 green")
  t1.green.index <- match(t1.green,probe)
  t1.green.index <- t1.green.index[which(!is.na(t1.green.index))]
  t1.green.signal <- c(m[t1.green.index],um[t1.green.index])
  t1.green.model <- Mclust(t1.green.signal, G=2, verbose = F)
  t1.green.mean <- t1.green.model$parameters$mean
  t1.green.sd <- sqrt(t1.green.model$parameters$variance$sigmasq)

  pIG <- apply(cbind(1 - pnorm(m[t1.green.index], t1.green.mean[1], t1.green.sd[1]),
                     1 - pnorm(um[t1.green.index], t1.green.mean[1], t1.green.sd[1])),
               1, min)

  if(verbose) message("Normalizing probe of type 1 green")
  temp = t1.green.signal
  # t1.green.signal[which(t1.green.model$classification == 1)] <- qnorm(pnorm(t1.green.signal[which(t1.green.model$classification == 1)],t1.green.mean[1], t1.green.sd[1]), ref$t1.green.ref.mean[1], ref$t1.green.ref.sd[1])
  t1.green.signal[which(t1.green.model$classification == 1 & pnorm(t1.green.signal, t1.green.mean[1], t1.green.sd[1]) < 0.9 & pnorm(t1.green.signal, t1.green.mean[1], t1.green.sd[1]) > 0.1)] <- qnorm(pnorm(t1.green.signal[which(t1.green.model$classification == 1 & pnorm(t1.green.signal, t1.green.mean[1], t1.green.sd[1]) < 0.9 & pnorm(t1.green.signal, t1.green.mean[1], t1.green.sd[1]) > 0.1)],t1.green.mean[1], t1.green.sd[1]), ref$t1.green.ref.mean[1], ref$t1.green.ref.sd[1])
  t1.green.signal[which(t1.green.model$classification == 2 & pnorm(t1.green.signal, t1.green.mean[2], t1.green.sd[2]) < 0.9 & pnorm(t1.green.signal, t1.green.mean[2], t1.green.sd[2]) > 0.1)] <- qnorm(pnorm(t1.green.signal[which(t1.green.model$classification == 2 & pnorm(t1.green.signal, t1.green.mean[2], t1.green.sd[2]) < 0.9 & pnorm(t1.green.signal, t1.green.mean[2], t1.green.sd[2]) > 0.1)],t1.green.mean[2], t1.green.sd[2]), ref$t1.green.ref.mean[2], ref$t1.green.ref.sd[2])

  # t1.green.signal[which(t1.green.model$classification == 2)] <- qnorm(pnorm(t1.green.signal[which(t1.green.model$classification == 2)],t1.green.mean[2], t1.green.sd[2]), ref$t1.green.ref.mean[2], ref$t1.green.ref.sd[2])

  t1.green.signal[which(t1.green.signal == Inf)] = temp[which(t1.green.signal == Inf)]

  m[t1.green.index] <- t1.green.signal[1:(length(t1.green.signal)/2)]
  um[t1.green.index] <- t1.green.signal[(length(t1.green.signal)/2+1):length(t1.green.signal)]

  if(verbose) message("Detecting p value")
  t2.index <- match(t2, probe)
  t2.index <- t2.index[which(!is.na(t2.index))]
  pII <- apply(cbind(1 - pnorm(um[t2.index], t1.red.mean[1], t1.red.sd[1]),
                     1 - pnorm(m[t2.index], t1.green.mean[1], t1.green.sd[1])),
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
  return(normalized.signal)
}
