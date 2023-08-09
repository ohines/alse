require(dplyr)
require(parallel)


alse_estimands <- function(alse_functions) {
  with(alse_functions, {
    n <- length(y)
    y_res <- y - mu
    a_res <- a - pi
    ya_res <- y_res * a_res
    aa_res <- a_res^2

    eta <- sum(aa_res) / n
    psi_1 <- sum(ya_res) / n / eta
    infl_1 <- (ya_res - psi_1 * aa_res) / eta
    std_err_1 <- sqrt(sum(infl_1^2) / n / n)

    infl_r_uncentred <- (ya_res - lambda_r * aa_res) * beta_inv_r + lambda_r
    psi_r <- sum(infl_r_uncentred) / n
    std_err_r <- sqrt((sum(infl_r_uncentred^2) / n - psi_r^2) / n)

    infl_d_uncentred <- (ya_res - lambda_d * aa_res) * beta_inv_d + lambda_d
    psi_d <- sum(infl_d_uncentred) / n
    std_err_d <- sqrt((sum(infl_d_uncentred^2) / n - psi_d^2) / n)

    # summary stats for debugging
    min_r <- min(beta_inv_r)
    max_r <- max(beta_inv_r)
    min_d <- min(beta_inv_d)
    max_d <- max(beta_inv_d)

    list(
      psi_0 = psi_1,
      std_err_0 = std_err_1,
      psi_r = psi_r,
      std_err_r = std_err_r,
      psi_d = psi_d,
      std_err_d = std_err_d,
      min_r = min_r,
      max_r = max_r,
      min_d = min_d,
      max_d = max_d
    )
  })
}


get_alse_functions <- function(y, a, x, y_new, a_new, x_new, fitfunc) {
  n <- length(y)
  n_new <- length(y_new)
  x_all <- rbind(x, x_new, deparse.level = 0)

  # fit the standard regression functions
  pi_hat <- fitfunc(a, x, x_all)
  mu_hat <- fitfunc(y, x, x_all)
  ya_hat <- fitfunc(y * a, x, x_new)
  aa_hat <- fitfunc(a^2, x, x_new)

  # fit the R-learner type learners
  a_res <- a - pi_hat$pred[1:n]
  v <- (y - mu_hat$pred[1:n]) / a_res
  w <- a_res^2
  lambda_hat <- fitfunc(v, x, x_new, weights = w)
  beta_inv_hat <- fitfunc(1 / w, x, x_new, weights = w)

  tibble(
    y = y_new,
    a = a_new,
    mu = mu_hat$pred[(n + 1):(n + n_new)],
    pi = pi_hat$pred[(n + 1):(n + n_new)],
    beta_inv_d = 1 / (aa_hat$pred - pi^2),
    lambda_d = (ya_hat$pred - (mu * pi)) * beta_inv_d,
    beta_inv_r = beta_inv_hat$pred,
    lambda_r = lambda_hat$pred,
  )
}


alse_analysis <- function(
    y,
    a,
    x,
    folds,
    fitfunc,
    parallel = FALSE) {
  alg_1 <- get_alse_functions(y, a, x, y, a, x, fitfunc)

  fn <- function(fold) {
    in_train <- folds != fold
    in_test <- !in_train
    get_alse_functions(
      y[in_train],
      a[in_train],
      x[in_train, ],
      y[in_test],
      a[in_test],
      x[in_test, ],
      fitfunc
    )
  }
  if (parallel) {
    alg_2 <- mclapply(
      unique(folds), fn,
      mc.cores = getOption("mc.cores", 10)
    ) %>% bind_rows()
  } else {
    alg_2 <- sapply(unique(folds), fn, simplify = FALSE) %>% bind_rows()
  }

  ests_1 <- alse_estimands(alg_1)
  ests_2 <- alse_estimands(alg_2)

  names(ests_1) <- paste0(names(ests_1), "_1")
  names(ests_2) <- paste0(names(ests_2), "_2")

  as.list(unlist(c(ests_1, ests_2)))
}
