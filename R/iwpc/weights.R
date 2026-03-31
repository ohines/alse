source("R/weight_approximation.R")
require(ggplot2)
require(dplyr)
require(tidyr)
require(glue)

DATA_DIR <- "Data"
RESULTS_DIR <- "Output"
PLOT_DISCRETE <- file.path(RESULTS_DIR, "weights_plot_discrete.pdf")
PLOT_NONDISCRETE <- file.path(RESULTS_DIR, "weights_plot_nondiscrete.pdf")

get_weights_df <- function(path, cross_validated, discrete) {
    epsilon <- 1e-8
    arrow::read_parquet(path) |>
        as_tibble() |>
        mutate(
            beta_inv_d_trimmed = pmax(beta_inv_d, epsilon),
            beta_inv_r_trimmed = pmax(beta_inv_r, epsilon),
            z_d = (a - pi) * sqrt(beta_inv_d_trimmed),
            z_r = (a - pi) * sqrt(beta_inv_r_trimmed),
            w_d0 = approximate_weights(z_d),
            w_r0 = approximate_weights(z_r),
            var_w_d = (1 / beta_inv_d_trimmed) *
                mean((a - pi)^2) /
                mean(1 / beta_inv_d_trimmed),
            var_w_r = (1 / beta_inv_r_trimmed) *
                mean((a - pi)^2) /
                mean(1 / beta_inv_r_trimmed),
            w_d = w_d0 * var_w_d,
            w_r = w_r0 * var_w_r,
            cross_validated = cross_validated,
            discrete = discrete,
        )
}

get_combined_weights_df <- function(discrete) {
    root <- ""
    if (discrete) {
        root <- "non"
    }
    path1 <- file.path(DATA_DIR, glue("{root}discrete_raw_1.parquet"))
    path2 <- file.path(DATA_DIR, glue("{root}discrete_raw_2.parquet"))
    df1 <- get_weights_df(path1, cross_validated = "noSS", discrete = discrete)
    df2 <- get_weights_df(path2, cross_validated = "SS", discrete = discrete)

    bind_rows(df1, df2) |>
        select(w_r, w_d, cross_validated, discrete) |>
        pivot_longer(
            cols = starts_with("w_"),
            names_to = "method",
            names_prefix = "w_",
            values_to = "weight"
        ) |>
        mutate(
            method = factor(method, levels = c("d", "r"), labels = c("A", "B"))
        )
}

get_facet_histogram <- function(df) {
    ggplot(df, aes(x = weight)) +
        geom_histogram(binwidth = 0.2) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        facet_grid(rows = vars(cross_validated), cols = vars(method)) +
        xlim(0, 5) +
        xlab("Weight") +
        ylab("Count") +
        theme_bw()
}

df_discrete <- get_combined_weights_df(discrete = TRUE)
df_nondiscrete <- get_combined_weights_df(discrete = FALSE)

ggsave(PLOT_DISCRETE, get_facet_histogram(df_discrete), width = 8.5, height = 6)
ggsave(
    PLOT_NONDISCRETE,
    get_facet_histogram(df_nondiscrete),
    width = 8.5,
    height = 6
)
