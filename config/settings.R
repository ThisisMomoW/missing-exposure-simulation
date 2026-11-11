study_settings <- function() list(
  mice_m = 10L, mice_maxit = 10L, mice_rf_ntree = 500L,
  min_valid_imputations = 2L,
  rf_backend = "ranger", nuisance_rf_trees = 500L,
  outer_folds = 10L, aipw_internal_folds = 9L, mr_internal_folds = 5L,
  probability_bound = .025, imputation_probability_bound = 1e-6,
  mice_ci = "Rubin t", other_ci = "normal, critical value 1.96",
  rng = c("Mersenne-Twister", "Inversion", "Rejection"))

validate_imputation_minimum <- function(minimum, requested = study_settings()$mice_m) {
  if (!is.numeric(minimum) || length(minimum)!=1L || !is.finite(minimum) ||
      minimum<2L || minimum>requested || minimum!=as.integer(minimum))
    stop("min_valid_imputations must be an integer between 2 and the requested imputation count")
  as.integer(minimum)
}

replication_seeds <- function(design, analysis, replication) {
  method_id <- match(analysis$method, study_methods)
  if (design$dgp == "nonlinear" && analysis$specification == "rf") method_id <- method_id + 5L
  if (analysis$specification == "interaction") method_id <- method_id + 6L
  data <- 123L + design$seed_cell*100000L + replication
  fit <- data + method_id*10000000L + if (design$dgp == "nonlinear") 90000000L else 0L
  list(data = as.integer(data), fit = as.integer(fit),
       imputation = as.integer(if (design$dgp == "nonlinear") fit else 123L+replication))
}

design_seeds <- function(design) list(
  truth = if (design$dgp == "gaussian") 9001L else 761239L,
  calibration = if (design$dgp == "gaussian") 991L else 123L+design$seed_cell*1000000L+777L)
