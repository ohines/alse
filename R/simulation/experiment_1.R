source("R/simulation/utils.R")
source("R/algorithms.R")
source("R/modelling_utilities.R")
require(tidyr)
require(dplyr)
require(glue)

RESULTS_DIR <- "Output"
SIM_NAME <- "experiment1"

# We parallelize the simulation already so we set this to avoid creating too many threads
# the default for ranger is to use 2 threads.
Sys.setenv(R_RANGER_NUM_THREADS = 1)


generate_data_simple <- function(n) {
    tibble(
        X1 = runif(n, -1, 1),
        X2 = runif(n, -1, 1),
        X3 = runif(n, -1, 1),
        A = rnorm(
            n,
            mean = X1 + 0.5 * X1^3 - 2 * X2^2 + X1^2 * X2,
            sd = 1 + X1^2
        ),
        Y = rnorm(n, mean = A * (1 + X1 - X1^2 - 0.5 * X2^2) + X2 * (X3 - X1^2))
    )
}


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
estimate_names <- c(
    "psi_0",
    "std_err_0",
    "psi_r",
    "std_err_r",
    "psi_d",
    "std_err_d",
    "min_r",
    "max_r",
    "min_d",
    "max_d"
)
estimate_names <- c(glue("{estimate_names}_1"), glue("{estimate_names}_2"))
estimate_names <- c("n", glue("{estimate_names}g"), glue("{estimate_names}r"))

# The meat of this script
res <- run_simulation(
    n_datasets = 1000,
    sample_sizes = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000),
    generate_data = generate_data_simple,
    get_estimates = get_estimates_simple,
    estimate_names = estimate_names,
    sim_name = SIM_NAME,
    results_directory = RESULTS_DIR,
    k_folds = 5,
    append = FALSE,
    mc_cores = 10
)
