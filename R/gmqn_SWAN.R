#' GMQN + SWAN
#'
#' Using SWAN to conduct between array normalization and using BMIQ to conduct within array normalization. We slightly modified swan function in minfi package. When the signal value is less than 0, we make it equal to the non-zero minimum value not the background signal intensity.
#' @param m The vector of methylation signal intensities.
#' @param um The vector of umethylation signal intensities.
#' @param probe the probe id for m and um.
#' @param type The type of methylation array which can be '450k'(default) or '850k'
#' @param ref The reference used to normalize data. By default, it uses the reference set up from GSE105018 which is also used in EWAS Datahub.
#' @return A data frame contains normalized m and um, p value, and DNA methylation level
gmqn_swan <- function(m, um, probe, type = '450k', ref = 'default', verbose = TRUE) {

  normalized.signal = gmqn_normalize(m, um, probe, type = '450k', ref = 'default', verbose = TRUE)

  if(verbose) message("Start SWAN normalizing")
  m = data.frame(m)
  swan.processed = .modified_SWAN(data.frame(m, row.names = probe), data.frame(um, row.names = probe), probe)
  normalized.signal$m = swan.processed$M
  normalized.signal$um = swan.processed$U
  beta <- swan.processed$M / (swan.processed$M + swan.processed$U)
  beta[which(normalized.signal$p >= 0.01)] = NA
  normalized.signal$beta = normalized.signal$beta
  return(normalized.signal)
}
