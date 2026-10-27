# Model and inference settings; see report/ANALYSIS_NOTES.md.
execution_settings <- list(
  version = "missing_exposure_v1",
  aipw_outer_folds = 10L, aipw_inner_cv = 9L,
  drcmd_outer_folds = 10L, drcmd_inner_cv = 5L,
  sl_ranger_trees = 500L, sl_imputation_inner_cv = 10L,
  mice_rf_trees = 500L, mice_m = 10L, mice_maxit = 10L,
  trim = 0.025, imputation_probability_epsilon = 1e-6,
  tr_bootstrap = 100L, tr_min_valid = 80L, tr_shrink = 0.99,
  rng = c("Mersenne-Twister", "Inversion", "Rejection")
)
