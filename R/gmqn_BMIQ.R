#' GMQN + BMIQ
#'
#' Using GMQN to conduct between array normalization and using BMIQ to conduct within array normalization. The BMIQ was developed by Andrew Teschendorff in 2013 and this implementation is the based on the version modified by Steve Horvath. I also removed and modified some parameters to simplify the function.
#' @param m The vector of methylation signal intensities.
#' @param um The vector of umethylation signal intensities.
#' @param probe the probe id for m and um.
#' @param type The type of methylation array which can be '450k'(default) or '850k'
#' @param ref The reference used to normalize data. By default, it uses the reference set up from GSE105018 which is also used in EWAS Datahub.
#' @return A data frame contains normalized m and um, p value, and DNA methylation level
gmqn_bmiq <- function(m, um, probe, type = '450k', ref = 'default', verbose = TRUE) {

  normalized.signal = gmqn_normalize(m, um, probe, type = type, ref = ref, verbose = TRUE)

  if(verbose) message("Start BMIQ normalizing")
  beta <- normalized.signal$m / (normalized.signal$m + normalized.signal$um)
  beta[which(normalized.signal$p < 0.01)] <- .BMIQ(beta[which(normalized.signal$p < 0.01)], type[which(p < 0.01)])
  beta[which(normalized.signal$p >= 0.01)] = NA
  normalized.signal$beta <- beta
  return(normalized.signal)
}
