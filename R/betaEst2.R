# Internal functions -----------------------------------------------------------

.betaEst2 <- function (y, w, weights) {
  yobs <- !is.na(y)
  if (sum(yobs) <= 1)
    return(c(1, 1))
  y <- y[yobs]
  w <- w[yobs]
  weights <- weights[yobs]
  N <- sum(weights * w)
  p <- sum(weights * w * y) / N
  v <- sum(weights * w * y * y) / N - p * p
  logab = log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p)/v -1))
  if (sum(yobs) == 2)
    return(exp(logab))
  opt <- try(optim(logab, betaObjf, ydata = y, wdata = w, weights = weights,
                   method = "Nelder-Mead", control = list(maxit = 50)), silent = TRUE)
  if (inherits(opt, "try-error"))
    return(c(1, 1))
  exp(opt$par)
}
