source("R/simulation/utils.R")
source("R/algorithms.R")
source("R/modelling_utilities.R")
require(tidyr)
require(dplyr)
require(ggplot2)
require(glue)
require(cowplot)
require(latex2exp)

RESULTS_DIR <- "Output/"
SIM_NAME <- "experiment1"


generate_data_simple <- function(n) {
  tibble(
    X1 = runif(n, -1, 1),
    X2 = runif(n, -1, 1),
    X3 = runif(n, -1, 1),
    A = rnorm(n, mean = X1 + 0.5 * X1^3 - 2 * X2^2 + X1^2 * X2, sd = 1 + X1^2),
    Y = rnorm(n, mean = A * (1 + X1 - X1^2 - 0.5 * X2^2) + X2 * (X3 - X1^2))
  )
}


true_vals <- tibble(
  estimand = c("Psi", "psi"),
  true_value = c(107 / 294, 1 / 2)
)


get_estimates_simple <- function(df, k_folds) {
  y <- df$Y
  a <- df$A
  x <- df[c("X1", "X2", "X3")]
  n <- nrow(df)

  # no need to randomize folds, since data already randomized
  folds <- rep(1:k_folds, length.out = n)

  res_g <- alse_analysis(y, a, x, folds, fitfunc_gam)
  res_r <- alse_analysis(y, a, x, folds, fitfunc_ranger)

  c(list(n = n), res_g, res_r)
}

# names for values returned by 'get_estimates_simple'
estimate_names <- c("psi", "std_err_psi", "psi_r", "std_err_r", "psi_d", "std_err_d", "psi_t", "std_err_t", "min_r", "max_r", "min_d", "max_d")
estimate_names <- c(glue("{estimate_names}_1"), glue("{estimate_names}_2"))
estimate_names <- c("n", glue("{estimate_names}g"), glue("{estimate_names}r"))

# The meat of this script
res <- run_simulation(
  n_datasets = 1, # 1000 in final run
  sample_sizes = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000),
  generate_data = generate_data_simple,
  get_estimates = get_estimates_simple,
  estimate_names = estimate_names,
  sim_name = SIM_NAME,
  results_directory = RESULTS_DIR,
  k_folds = 5,
  append = TRUE
)
