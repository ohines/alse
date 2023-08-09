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
