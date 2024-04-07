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
PLOT_OUTPUT <- paste0(RESULTS_DIR, "sim_plot.pdf")


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
estimate_names <- c("psi_0", "std_err_0", "psi_r", "std_err_r", "psi_d", "std_err_d", "min_r", "max_r", "min_d", "max_d")
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
  append = FALSE
)

# process and tabulate results
df_r <- res %>%
  pivot_longer(
    cols = c(starts_with("psi"), starts_with("max"), starts_with("min"), starts_with("std_err")),
    names_to = c(".value", "method", "cv", "learner"),
    names_pattern = "(.{3,7})_(.{1})_(.{1})(.{1})$"
  ) %>%
  mutate(
    estimand = ifelse(method == "0", "Psi", "psi"),
    learner = factor(
      learner,
      levels = c("g", "r"),
      labels = c("gam", "ranger")
    ),
    method = factor(
      method,
      levels = c("0", "r", "d"),
      labels = c("Psi", "r-learner", "direct")
    ),
    cv = ifelse(cv == "2", "SS", "noSS"),
  ) %>%
  left_join(true_vals) %>%
  mutate(
    bias = psi - true_value,
    scaled_bias = bias * sqrt(n),
    coverage_yn = abs(bias) <= 1.959964 * std_err
  )

summary_stats <- df_r %>%
  group_by(n, method, learner, estimand, cv) %>%
  summarise(
    n_samples = length(bias),
    coverage = mean(coverage_yn),
    root_n_bias = mean(scaled_bias),
    n_var = mean(scaled_bias^2) - root_n_bias^2,
    root_n_bias_std_err = sqrt(n_var / n_samples),
    bias_min = root_n_bias - 1.96 * root_n_bias_std_err,
    bias_max = root_n_bias + 1.96 * root_n_bias_std_err
  )

# pretty plots
plots_est <- function(summary_stats, est) {
  df_plot <- summary_stats %>%
    rename(
      Learner = learner,
      Algorithm = cv,
      Coverage = coverage
    ) %>%
    filter(method == est)

  p1 <- df_plot %>% ggplot(
    aes(
      x = n,
      y = root_n_bias,
      color = Learner,
      shape = Algorithm,
      ymin = bias_min,
      ymax = bias_max
    )
  ) +
    ylab(TeX("$\\sqrt{n}$(Bias)")) +
    geom_hline(aes(yintercept = 0)) +
    geom_point() +
    geom_errorbar(width = 100) +
    theme_bw()

  p2 <- df_plot %>% ggplot(
    aes(x = n, y = Coverage, color = Learner, shape = Algorithm)
  ) +
    geom_hline(aes(yintercept = 0.95)) +
    geom_point() +
    ylim(0.62, 1) +
    theme_bw()


  p3 <- df_plot %>% ggplot(
    aes(x = n, y = n_var, color = Learner, shape = Algorithm)
  ) +
    ylab(TeX("$n$(Emp. variance)")) +
    geom_hline(aes(yintercept = 0)) +
    geom_point() +
    theme_bw()

  list(p1, p3, p2)
}


psi_r <- plots_est(summary_stats, "r-learner")
psi_d <- plots_est(summary_stats, "direct")
psi <- plots_est(summary_stats, "Psi")


nice_display <- function(row1, row2, row3, labels = "AUTO") {
  # https://github.com/wilkelab/cowplot/blob/master/vignettes/shared_legends.Rmd
  legend_b <- cowplot::get_legend(
    row1[[2]] +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  prow <- cowplot::plot_grid(
    row1[[1]] + theme(legend.position = "none"),
    row1[[2]] + theme(legend.position = "none"),
    row1[[3]] + theme(legend.position = "none"),
    row2[[1]] + theme(legend.position = "none"),
    row2[[2]] + theme(legend.position = "none"),
    row2[[3]] + theme(legend.position = "none"),
    row3[[1]] + theme(legend.position = "none"),
    row3[[2]] + theme(legend.position = "none"),
    row3[[3]] + theme(legend.position = "none"),
    ncol = 3,
    align = "hv",
    axis = "b",
    labels = labels
  )
  cowplot::plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
}

plt <- nice_display(psi, psi_d, psi_r)

pdf(file = PLOT_OUTPUT, height = 9, width = 8.5)
print(plt)
dev.off()
