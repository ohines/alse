require(tidyr)
require(dplyr)
require(ggplot2)
require(glue)
require(cowplot)
require(latex2exp)

RESULTS_DIR <- "Output"
SIM_NAME <- "experiment1"
PLOT_OUTPUT <- file.path(RESULTS_DIR, "sim_plot.pdf")

results_path <- file.path(RESULTS_DIR, glue("{SIM_NAME}.parquet"))
res <- arrow::read_parquet(results_path)

true_vals <- tibble(
    estimand = c("Psi", "psi"),
    true_value = c(107 / 294, 1 / 2)
)

# process and tabulate results
df_r <- res |>
    pivot_longer(
        cols = c(
            starts_with("psi"),
            starts_with("max"),
            starts_with("min"),
            starts_with("std_err")
        ),
        names_to = c(".value", "method", "cv", "learner"),
        names_pattern = "(.{3,7})_(.{1})_(.{1})(.{1})$"
    ) |>
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
    ) |>
    left_join(true_vals) |>
    mutate(
        bias = psi - true_value,
        scaled_bias = bias * sqrt(n),
        coverage_yn = abs(bias) <= 1.959964 * std_err
    )

summary_stats <- df_r |>
    group_by(n, method, learner, estimand, cv) |>
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
    df_plot <- summary_stats |>
        rename(Learner = learner, Algorithm = cv, Coverage = coverage) |>
        filter(method == est)

    p1 <- df_plot |>
        ggplot(aes(
            x = n,
            y = root_n_bias,
            color = Learner,
            shape = Algorithm,
            ymin = bias_min,
            ymax = bias_max
        )) +
        ylab(TeX("$\\sqrt{n}$(Bias)")) +
        geom_hline(aes(yintercept = 0)) +
        geom_point() +
        geom_errorbar(width = 100) +
        theme_bw()

    p2 <- df_plot |>
        ggplot(aes(x = n, y = Coverage, color = Learner, shape = Algorithm)) +
        geom_hline(aes(yintercept = 0.95)) +
        geom_point() +
        ylim(0.62, 1) +
        theme_bw()

    p3 <- df_plot |>
        ggplot(aes(x = n, y = n_var, color = Learner, shape = Algorithm)) +
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
ggsave(PLOT_OUTPUT, plt, width = 8.5, height = 9)
