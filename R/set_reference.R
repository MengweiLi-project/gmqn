#' Set reference for GMQN using user defined data
#'
#' @param m The data frame of methylation signal intensities.
#' @param um The data frame of umethylation signal intensities.
#' @param type The type of methylation array which can be '450k'(default) or '850k'
#' @return A list contains referenced parameters
set_reference <- function(m, um, type = '450k') {

  set.seed(1)

  if (sum(row.names(m) == row.names(um)) != dim(m)[1]) {
    stop('The rownames of m and um are different!')
  }
  print('Processing the input data---------------------------------------')
  probe = row.names(m)
  m = as.numeric(apply(m, 1, mean))
  um = as.numeric(apply(um, 1, mean))

  if (type == '450k') {
    t1.red <- row.names(annon_450k)[which(annon_450k$color == 'Red' & annon_450k$probe_type == 'I')]
    t1.green <- row.names(annon_450k)[which(annon_450k$color == 'Grn' & annon_450k$probe_type == 'I')]
    type <- annon_450k[probe, 'probe_type']
  } else if (type == '850k') {
    t1.red <- row.names(annon_850k)[which(annon_850k$color == 'Red' & annon_850k$probe_type == 'I')]
    t1.green <- row.names(annon_850k)[which(annon_850k$color == 'Grn' & annon_850k$probe_type == 'I')]
    type <- annon_850k[probe, 'probe_type']
  }

  t1.red.index <- match(t1.red, probe)
  t1.red.signal <- c(m[t1.red.index], um[t1.red.index])
  t1.red.model <- Mclust(t1.red.signal, G=2, verbose = F)
  t1.red.ref.mean <- t1.red.model$parameters$mean
  t1.red.ref.sd <- sqrt(t1.red.model$parameters$variance$sigmasq)

  t1.green.index <- match(t1.green,probe)
  t1.green.signal <- c(m[t1.green.index],um[t1.green.index])
  t1.green.model <- Mclust(t1.green.signal, G=2, verbose = F)
  t1.green.ref.mean <- t1.green.model$parameters$mean
  t1.green.ref.sd <- sqrt(t1.green.model$parameters$variance$sigmasq)

  return(list(t1.green.ref.mean = t1.green.ref.mean,
              t1.green.ref.sd = t1.green.ref.sd,
              t1.red.ref.mean = t1.red.ref.mean,
              t1.red.ref.sd = t1.red.ref.sd))

}
