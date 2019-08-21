# Internal functions -----------------------------------------------------------

.blc2 <- function(Y, w, maxiter = 25, tol = 1e-06, weights = NULL, verbose = FALSE) {
  Ymn <- min(Y[Y > 0], na.rm = TRUE)
  Ymx <- max(Y[Y < 1], na.rm = TRUE)
  Y <- pmax(Y, Ymn/2)
  Y <- pmin(Y, 1 - (1 - Ymx) / 2)
  Yobs <- !is.na(Y)
  J <- dim(Y)[2]
  K <- dim(w)[2]
  n <- dim(w)[1]
  if (n != dim(Y)[1])
    stop("Dimensions of w and Y do not agree")
  if (is.null(weights))
    weights = rep(1, n)
  mu = a = b = matrix(Inf, K, J)
  crit <- Inf
  for (i in 1:maxiter) {
    warn0 <- options()$warn
    options(warn = -1)
    eta <- apply(weights * w, 2, sum)/sum(weights)
    mu0 <- mu
    for (k in 1:K) {
      for (j in 1:J) {
        ab <- .betaEst2(Y[, j], w[, k], weights)
        a[k, j] <- ab[1]
        b[k, j] <- ab[2]
        mu[k, j] <- ab[1]/sum(ab)
      }
    }
    ww <- array(0, dim = c(n, J, K))
    for (k in 1:K) {
      for (j in 1:J) {
        ww[Yobs[, j], j, k] <- dbeta(Y[Yobs[, j], j],
                                    a[k, j], b[k, j], log = TRUE)
      }
    }
    options(warn = warn0)
    w <- apply(ww, c(1, 3), sum, na.rm = TRUE)
    wmax <- apply(w, 1, max)
    for (k in 1:K) w[, k] = w[, k] - wmax
    w <- t(eta * t(exp(w)))
    like <- apply(w, 1, sum)
    w <- (1/like) * w
    llike <- weights * (log(like) + wmax)
    crit <- max(abs(mu - mu0))
    if (verbose)
      print(crit)
    if (crit < tol)
      break
  }
  return(list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike)))
}

