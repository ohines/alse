source("R/weight_approximation.R")
require(ggplot2)
require(dplyr)
require(tidyr)
require(glue)

DATA_DIR <- "Data/"

get_weights_df <- function(path, cross_validated, discrete) {
  arrow::read_parquet(path) %>%
    mutate(
      inv_sd_d = sqrt(pmax(beta_inv_d, 0)),
      inv_sd_r = sqrt(pmax(beta_inv_r, 0)),
      z_d = (a - pi) * inv_sd_d,
      z_r = (a - pi) * inv_sd_r,
      f_d = density_point(z_d),
      f_r = density_point(z_r),
      w_d = approximate_weights(z_d, f_d),
      w_r = approximate_weights(z_r, f_r),
      cross_validated = cross_validated,
      discrete = discrete,
    )
}

get_combined_weights_df <- function(discrete) {
  root <- ""
  if (discrete) root <- "non"
  path1 <- glue("{DATA_DIR}{root}discrete_raw_1.parquet")
  path2 <- glue("{DATA_DIR}{root}discrete_raw_2.parquet")
  df1 <- get_weights_df(path1, cross_validated = "noSS", discrete = discrete)
  df2 <- get_weights_df(path2, cross_validated = "SS", discrete = discrete)

  bind_rows(df1, df2) %>%
    select(w_r, w_d, cross_validated, discrete) %>%
    pivot_longer(
      cols = starts_with("w_"),
      names_to = "method",
      names_prefix = "w_",
      values_to = "weight"
    ) %>%
    mutate(method = factor(method, levels = c("d", "r"), labels = c("A", "B")))
}

get_facet_histogram <- function(df) {
  ggplot(df, aes(x = weight)) +
    geom_histogram(binwidth = 0.2) +
    facet_grid(rows = vars(cross_validated), cols = vars(method)) +
    xlim(0, 5) +
    theme_bw()
}

plt1 <- get_facet_histogram(get_combined_weights_df(discrete = TRUE))
plt2 <- get_facet_histogram(get_combined_weights_df(discrete = FALSE))

arrow::read_parquet("Output/data_illustration.parquet")
