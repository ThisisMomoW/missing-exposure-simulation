SL.study_interaction <- function(Y, X, newX, family, obsWeights, ...) {
  X <- as.data.frame(X); newX <- as.data.frame(newX)
  formula <- if (all(c("A", "X1", "X2", "X3") %in% names(X)))
    Y ~ A + X1 + X2 + X3 + A:I(X1 + X2 + X3) else Y ~ .
  fit <- stats::glm(formula, data = data.frame(Y = Y, X), family = family, weights = obsWeights)
  object <- structure(list(object = fit), class = "SL.study_interaction")
  list(pred = as.numeric(stats::predict(fit, newdata = newX, type = "response")), fit = object)
}

predict.SL.study_interaction <- function(object, newdata, ...) {
  as.numeric(stats::predict(object$object, newdata = as.data.frame(newdata), type = "response"))
}

SL.study_rf <- function(Y, X, newX, family, obsWeights, ...) {
  SuperLearner::SL.ranger(Y, X, newX, family, obsWeights,
    num.trees = 500L, mtry = floor(sqrt(ncol(X))),
    min.node.size = if (family$family == "gaussian") 5L else 1L,
    probability = family$family == "binomial", replace = TRUE, sample.fraction = 1,
    num.threads = 1L, write.forest = TRUE, verbose = FALSE)
}

# drcmd names the fully observed outcome predictor y; the outcome response is Y.
SL.study_without_y <- function(Y, X, newX, family, obsWeights, ...) {
  X <- as.data.frame(X); newX <- as.data.frame(newX)
  X <- X[, names(X) != "y", drop = FALSE]
  newX <- newX[, names(newX) != "y", drop = FALSE]
  SuperLearner::SL.glm(Y, X, newX, family, obsWeights, ...)
}

outcome_learner <- function(specification) {
  switch(specification, main = "SL.glm", glm = "SL.glm",
    interaction = "SL.study_interaction", rf = "SL.study_rf",
    stop("Unknown model specification"))
}
