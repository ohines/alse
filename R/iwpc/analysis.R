require(tidyverse)
require(SuperLearner)
source("R/alse.R")
source("R/fitfuncs.R")


clean_data_path <- "Data/iwpc_clean.parquet"

df <- arrow::read_parquet(clean_data_path, as_tibble = TRUE) %>%
    rename(A = Dose, Y = INR)

x <- select(df, -Y, -A)

## Basic parametric analysis
mod <- glm(Y ~ ., data = df, family = gaussian())
summary(mod)

set.seed(123456)
fold_list <- get_fold_list(nrow(x), 20)

## Basic Alse analysis using glm
ests_glm <- ALSE_analysis(df$Y, df$A, x, fitfunc_glm, fold_list = fold_list)

## Computationally intensive analysis using super learner
ests_sl <- ALSE_analysis(
    df$Y, df$A, x, fitfunc_SL1,
    parallel = TRUE, fold_list = fold_list
)
saveRDS(ests_sl, "Data/nondiscreteFit.rds")

## Computationally intensive analysis using discrete super learner
ests_sl2 <- ALSE_analysis(
    df$Y, df$A, x, fitfunc_SL2,
    parallel = TRUE, fold_list = fold_list
)
saveRDS(ests_sl2, "Data/discreteFit.rds")


a_sl <- as.data.frame(t(ests_sl2))
a_sl$lower <- a_sl$Estimate - 1.96 * a_sl$Std.err
a_sl$upper <- a_sl$Estimate + 1.96 * a_sl$Std.err
a_sl$Tsq <- (a_sl$Estimate / a_sl$Std.err)^2
a_sl$pval <- pchisq(a_sl$Tsq, df = 1, lower.tail = FALSE)
a_sl
signif(a_sl[, 1:4] * 1000, 3)
