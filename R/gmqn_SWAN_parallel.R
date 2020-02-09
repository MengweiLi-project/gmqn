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
  registerDoParallel(min(ncol(m), ncpu))


  if (type == '450k') {
    m.tem = matrix(nrow = 485512, ncol = ncol(m))
    m.tem = data.frame(m.tem)
    row.names(m.tem) = row.names(annon_450k)
  } else if (type == '850k') {
    m.tem = matrix(nrow = 866836, ncol = ncol(m))
    m.tem = data.frame(m.tem)
    row.names(m.tem) = row.names(annon_850k)
  } else {
    stop("Incorrect platform type!")
  }

  names(m.tem) = names(m)
  um.tem = m.tem
  m.tem[intersect(row.names(m), row.names(m.tem)),] = m[intersect(row.names(m), row.names(m.tem)),]
  um.tem[intersect(row.names(um), row.names(m.tem)),] = um[intersect(row.names(um), row.names(um.tem)),]
  m = m.tem
  um = um.tem

  beta.GMQN.swan = foreach (i=1:ncol(m), .combine=cbind) %dopar% {
    beta = rep(NA, nrow(m))
    no.na.index = which(is.na(m[,i]) == F & is.na(um[,i]) == F)
    res = gmqn::gmqn_swan(m[no.na.index,i], um[no.na.index,i], row.names(m)[no.na.index], type = type, ref = ref)
    beta[no.na.index] = res$beta
    beta
  }
  beta.GMQN.swan = data.frame(beta.GMQN.swan)
  names(beta.GMQN.swan) = names(m)
  row.names(beta.GMQN.swan) = row.names(m)
  return(beta.GMQN.swan)
}
