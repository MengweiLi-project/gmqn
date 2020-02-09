#' BMIQ in parallel
#'
#' Using BMIQ to conduct within array normalization.
#' @param m The vector of methylation signal intensities.
#' @param um The vector of umethylation signal intensities.
#' @param probe the probe id for m and um.
#' @param type The type of methylation array which can be '450k'(default) or '850k'
#' @param ref The reference used to normalize data. By default, it uses the reference set up from GSE105018 which is also used in EWAS Datahub.
#' @return A data frame contains normalized m and um, p value, and DNA methylation level
bmiq_parallel <- function(m, um, type = '450k', ref = 'default', ncpu = 4, verbose = TRUE) {

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

  beta.BMIQ = foreach (i=1:ncol(m), .combine=cbind) %dopar% {
    beta = rep(NA, nrow(m))
    no.na.index = which(is.na(m[,i]) == F & is.na(um[,i]) == F)
    beta.tem = m[no.na.index,i]/(m[no.na.index,i]+um[no.na.index,i])
    design = CpG.counts[row.names(m)[no.na.index],"Type"]
    res = gmqn::BMIQ(beta.tem, design)
    beta[no.na.index] = res
    beta
  }
  beta.BMIQ = data.frame(beta.BMIQ)
  names(beta.BMIQ) = names(m)
  row.names(beta.BMIQ) = row.names(m)
  return(beta.BMIQ)
}
