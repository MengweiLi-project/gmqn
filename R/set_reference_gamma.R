#' Set reference for GMQN using user defined data
#'
#' @param m The data frame of methylation signal intensities.
#' @param um The data frame of umethylation signal intensities.
#' @param type The type of methylation array which can be '450k'(default) or '850k'
#' @return A list contains referenced parameters
set_reference_gamma <- function(m, um, type = '450k') {

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
  t1.red.model <- gammamixEM(t1.red.signal[which(t1.red.signal > 0)][1:10000], lambda = c(0.5, 0.5))
  t1.red.ref.alpha <- as.numeric(t1.red.model$gamma.pars[1,])
  t1.red.ref.beta <- as.numeric(t1.red.model$gamma.pars[2,])

  t1.green.index <- match(t1.green,probe)
  t1.green.signal <- c(m[t1.green.index],um[t1.green.index])
  t1.green.model <- gammamixEM(t1.green.signal[which(t1.green.signal > 0)][1:10000], lambda = c(0.5, 0.5))
  t1.green.ref.alpha <- as.numeric(t1.green.model$gamma.pars[1,])
  t1.green.ref.beta <- as.numeric(t1.green.model$gamma.pars[2,])

  return(list(t1.red.ref.alpha = t1.red.ref.alpha,
              t1.red.ref.beta = t1.red.ref.beta,
              t1.green.ref.alpha = t1.green.ref.alpha,
              t1.green.ref.beta = t1.green.ref.beta))

}
