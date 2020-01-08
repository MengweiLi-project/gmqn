#' GMQN + SWAN in parallel
#'
#' Using SWAN to conduct between array normalization and using BMIQ to conduct within array normalization. We slightly modified swan function in minfi package. When the signal value is less than 0, we make it equal to the non-zero minimum value not the background signal intensity.
#' @param m The vector of methylation signal intensities.
#' @param um The vector of umethylation signal intensities.
#' @param probe the probe id for m and um.
#' @param type The type of methylation array which can be '450k'(default) or '850k'
#' @param ref The reference used to normalize data. By default, it uses the reference set up from GSE105018 which is also used in EWAS Datahub.
#' @return A data frame contains normalized m and um, p value, and DNA methylation level
gmqn_swan_parallel <- function(m, um, type = '450k', ref = 'default', ncpu = 4, verbose = TRUE) {

  registerDoParallel(max(ncol(m), ncpu))
  beta.GMQN.swan = foreach (i=1:dim(m)[2], .combine=cbind) %dopar% {
    res = gmqn::gmqn_swan(m[, i], um[, i], row.names(m), type = type, ref = ref)
    res$beta
  }
  beta.GMQN.swan = data.frame(beta.GMQN.swan)
  names(beta.GMQN.swan) = names(m)
  row.names(beta.GMQN.swan) = row.names(m)
  return(res)
}
