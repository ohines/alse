# Modified version of stats::density.default
# Returns the density at each point of x instead of for a grid of points
# This function is not used by the estimation algorithms
# but is instead used to approximate the weights
density_point <- function(
    x, bw = "nrd0", kernel = c(
      "gaussian",
      "epanechnikov", "rectangular", "triangular", "biweight",
      "cosine", "optcosine"
    )) {
  kernel <- match.arg(kernel)
  nx <- length(x)
  weights <- rep.int(1 / nx, nx)
  bw <- switch(tolower(bw),
    nrd0 = bw.nrd0(x),
    nrd = bw.nrd(x),
    ucv = bw.ucv(x),
    bcv = bw.bcv(x),
    sj = ,
    `sj-ste` = bw.SJ(x,
      method = "ste"
    ),
    `sj-dpi` = bw.SJ(x, method = "dpi"),
    stop("unknown bandwidth rule")
  )
  n <- 2^ceiling(log2(nx))
  lo <- min(x) - 7 * bw
  up <- max(x) + 7 * bw
  y <- .Call(stats:::C_BinDist, x, weights, lo, up, n)
  kords <- seq.int(0, 2 * (up - lo), length.out = 2L * n)
  kords[(n + 2):(2 * n)] <- -kords[n:2]
  kords <- switch(kernel,
    gaussian = dnorm(kords, sd = bw),
    rectangular = {
      a <- bw * sqrt(3)
      ifelse(abs(kords) < a, 0.5 / a, 0)
    },
    triangular = {
      a <- bw * sqrt(6)
      ax <- abs(kords)
      ifelse(ax < a, (1 - ax / a) / a, 0)
    },
    epanechnikov = {
      a <- bw * sqrt(5)
      ax <- abs(kords)
      ifelse(ax < a, 3 / 4 * (1 - (ax / a)^2) / a, 0)
    },
    biweight = {
      a <- bw * sqrt(7)
      ax <- abs(kords)
      ifelse(ax < a, 15 / 16 * (1 - (ax / a)^2)^2 / a, 0)
    },
    cosine = {
      a <- bw / sqrt(1 / 3 - 2 / pi^2)
      ifelse(abs(kords) < a, (1 + cos(pi * kords / a)) / (2 * a), 0)
    },
    optcosine = {
      a <- bw / sqrt(1 - 8 / pi^2)
      ifelse(abs(kords) < a, pi / 4 * cos(pi * kords / (2 * a)) / a, 0)
    }
  )
  kords <- fft(fft(y) * Conj(fft(kords)), inverse = TRUE)
  kords <- pmax.int(0, Re(kords)[1L:n] / length(y))
  xords <- seq.int(lo, up, length.out = n)
  approx(xords, kords, x)$y
}


approximate_weights <- function(z, ...) {
  f <- density_point(z, ...)
  ord <- order(z)
  inverse_premutation <- order(ord)
  z_sorted <- z[ord]
  f <- f[ord]
  w <- - cumsum(z_sorted) / length(z) / f
  w[inverse_premutation]
}
