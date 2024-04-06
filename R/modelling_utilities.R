fitfunc_gam <- function(y, x, x_new, weights = NULL) {
  fit <- mgcv::gam(
    Y ~ ti(X1) + ti(X2) + ti(X3) + ti(X1, X2) + ti(X1, X3) + ti(X2, X3),
    data = cbind(Y = y, x),
    weights = weights
  )
  list(pred = predict(fit, x_new), fit = fit)
}

fitfunc_ranger <- function(y, x, x_new, weights = NULL) {
  fit <- ranger::ranger(y = y, x = x, case.weights = weights, num.trees = 1000)
  list(pred = predict(fit, data = x_new)$predictions, fit = fit)
}

fitfunc_glm <- function(y, x, x_new, weights = NULL) {
  fit <- glm(y ~ ., data = x, weights, family = gaussian())
  list(pred = predict(fit, x_new), fit = fit)
}


RANGERlearners <- create.Learner("SL.ranger",
  params = list(num.trees = 2000),
  tune = list(mtry = c(3, 4)),
  name_prefix = "RANGER"
)

GAMlearners <- create.Learner("SL.gam",
  tune = list(deg.gam = c(2, 3, 4)),
  name_prefix = "GAM"
)

XGBlearners <- create.Learner("SL.xgboost",
  params = list(minobspernode = 10, ntrees = 2000, shrinkage = 0.01),
  tune = list(max_depth = c(2, 3)),
  name_prefix = "XGB"
)


fitfunc_sl1 <- function(y, x, x_new, weights = NULL) {
  # uses non-discrete SL
  sl_library <- c(
    "SL.glmnet", "SL.glm",
    RANGERlearners$names,
    GAMlearners$names,
    XGBlearners$names
  )
  sl_folds <- 20

  fit <- SuperLearner::SuperLearner(y, x,
    family = gaussian(),
    SL.library = sl_library,
    obsWeights = weights,
    cvControl = list(V = sl_folds)
  )
  preds <- predict(fit, x_new, onlySL = TRUE)$pred[, 1]
  list(pred = preds, fit = fit)
}


fitfunc_sl2 <- function(y, x, x_new, weights = NULL) {
  # uses discrete SL
  sl_library <- c(
    "SL.glmnet", "SL.glm",
    RANGERlearners$names,
    GAMlearners$names,
    XGBlearners$names
  )
  sl_folds <- 20

  fit <- SuperLearner::SuperLearner(y, x,
    family = gaussian(),
    SL.library = sl_library,
    obsWeights = weights,
    cvControl = list(V = sl_folds)
  )
  preds <- predict(fit, x_new)$library.predict[, which.min(fit$cvRisk)]
  list(pred = preds, fit = fit)
}
