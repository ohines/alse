require(tidyverse)
require(SuperLearner)
require(arrow)
source("R/algorithms.R")
source("R/modelling_utilities.R")


clean_data_path <- "Data/iwpc_clean.parquet"

df <- read_parquet(clean_data_path, as_tibble = TRUE) %>%
    rename(A = Dose, Y = INR)

x <- select(df, -Y, -A)

# Basic parametric analysis
mod <- glm(Y ~ ., data = df, family = gaussian())
summary(mod)

# Randomly allocate folds
set.seed(123456)
n <- nrow(x)
n_folds <- 20
folds <- sample(rep_len(seq_len(n_folds), length.out = n))

# Basic Alse analysis using glm
ests_glm <- alse_analysis(df$Y, df$A, x, folds, fitfunc_glm)

# Computationally intensive analysis using super learner
ests_sl <- alse_analysis(
    df$Y,
    df$A,
    x,
    folds,
    fitfunc_sl1,
    parallel = TRUE
)
saveRDS(ests_sl, "Data/nondiscreteFit.rds")

# Computationally intensive analysis using discrete super learner
ests_sl2 <- ALSE_analysis(
    df$Y,
    df$A,
    x,
    folds,
    fitfunc_sl2,
    parallel = TRUE
)
saveRDS(ests_sl2, "Data/discreteFit.rds")

# collate results
results <- list(
    c(ests_glm, list(learner = "glm")),
    c(ests_sl, list(learner = "sl_non_discrete")),
    c(ests_sl2, list(learner = "sl_discrete"))
) %>%
    map(function(x) {
        df <- as.data.frame(x)
        colnames(df) <- names(x)
        df
    }) %>%
    bind_rows() %>%
    as_tibble()

df_r <- results %>%
    pivot_longer(
        cols = c(starts_with("psi"), starts_with("max"), starts_with("min"), starts_with("std_err")),
        names_to = c(".value", "method", "cv"),
        names_pattern = "(.{3,7})_(.{1})_(.{1})$"
    ) %>%
    rename(estimate = psi) %>%
    mutate(
        estimand = ifelse(method == "0", "Psi", "psi"),
        method = factor(
            method,
            levels = c("0", "r", "d"),
            labels = c("Psi", "r-learner", "direct")
        ),
        lower = estimate - 1.96 * std_err,
        upper = estimate + 1.96 * std_err,
        t_stat = abs(estimate / std_err),
        pval = pchisq(t_stat^2, df = 1, lower.tail = FALSE),
    )

# save results
write_parquet(df_r, "Output/data_illustration.parquet")
