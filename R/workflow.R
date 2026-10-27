load_implementation <- function(family) {
  # AIPW discovers custom SuperLearner wrappers in the global environment.
  for (f in c("dgp_functions", "trme_rd", "estimator_functions"))
    sys.source(file.path(root, "R", paste0(f, ".R")), .GlobalEnv)
  cfg <- load_config(family)
  list2env(as.list(cfg), .GlobalEnv)
  sys.source(file.path(root, "config/execution.R"), .GlobalEnv)
  if (family == "outcome_interaction") sys.source(file.path(root, "R/rf_estimators.R"), .GlobalEnv)
  if (family == "nonlinear_binary") {
    for (f in c("nonlinear_dgp", "nonlinear_estimators", "nonlinear_tr", "nonlinear_tr_wrapper"))
      sys.source(file.path(root, "R", paste0(f, ".R")), .GlobalEnv)
    assign("tr_settings", list(B = 100L, shrink_rate = .99), .GlobalEnv)
  }
  # Fail on changed software defaults instead of quietly changing the experiment.
  stopifnot(formals(AIPW::AIPW$public_methods$initialize)$k_split == 10L,
            formals(SuperLearner::SL.ranger)$num.trees == 500L,
            execution_settings$mice_rf_trees == 500L,
            formals(drcmd::drcmd)$cv_folds == 5L)
  invisible(cfg)
}
cell_constants <- function(family, sid, sc) {
  if (family == "nonlinear_binary") {
    truth <- calculate_truth(sc$gamma_int)
    alpha <- calibrate_r_intercept(sc$target_missing_share, sc$gamma_int,
      sc$gamma_y, seed = base_seed + sid * 1000000L + 777L)
    list(true_rd = truth$true_rd, R_intercept = alpha, truth = truth)
  } else {
    alpha <- do.call(calibrate_R_intercept, as.list(sc[c("target_missing_share", "gamma_int", "missing_mechanism", "gamma_y")]))
    list(true_rd = true_rd(sc$gamma_int), R_intercept = alpha)
  }
}
generate_replication <- function(family, sid, sim, sc, constants) {
  set.seed(base_seed + sid * 100000L + sim)
  if (family == "nonlinear_binary") {
    generate_one_dataset(sc$n, sc$gamma_int, sc$gamma_y, constants$R_intercept)
  } else {
    dat <- simulate_dataset(sc$n, sc$gamma_int, sc$missing_mechanism, sc$gamma_y, constants$R_intercept)
    attr(dat, "missing_mechanism") <- as.character(sc$missing_mechanism)
    dat
  }
}
seed_for <- function(family, sid, sim, method, spec) {
  if (family == "nonlinear_binary" && method == "tr_aipw") {
    fit <- as.integer(500000000 + sid * 1000000 + sim * 1000)
    return(c(fit = fit, internal = fit))
  }
  mid <- match(method, methods)
  offset <- 0L
  if (family == "nonlinear_binary") {
    mid <- match(method, methods[1:5]) + if (spec == "rf") 5L else 0L
    offset <- 90000000L
  } else if (spec == "rf") offset <- 90000000L
  else if (spec == "interaction") mid <- mid + 6L
  fit <- base_seed + sid * 100000L + sim + mid * 10000000L + offset
  c(fit = fit, internal = if (offset == 0L) base_seed + sim else fit)
}
fit_replication <- function(family, sid, sim, method, spec, dat, sc) {
  seeds <- seed_for(family, sid, sim, method, spec)
  set.seed(seeds[["fit"]])
  warnings <- character()
  rr <- withCallingHandlers(tryCatch({
    if (family == "nonlinear_binary") {
      if (method == "tr_aipw") fit_tr_binary(dat, sc, seeds[["internal"]])
      else fit_one(method, spec, dat, seeds[["internal"]], mice_rf_trees=execution_settings$mice_rf_trees)
    } else if (spec == "rf") {
      dd <- make_analysis_data(dat)
      switch(method,
        full_aipw = run_aipw_rf_result(dd$full, trim_cutoff),
        cc_aipw = run_aipw_rf_result(dd$complete_cases, trim_cutoff),
        single_imp_aipw = run_single_imp_rf_result(dd$missing, trim_cutoff),
        mice_aipw = run_mice_full_rf_result(dd$missing, mice_m, mice_maxit, seeds[["internal"]], trim_cutoff,
                                         ntree=execution_settings$mice_rf_trees),
        drcmd = run_drcmd(dd$missing, trim_cutoff))
    } else {
      fit_one_method(method, dat, seed = seeds[["internal"]], mice_m = mice_m,
        mice_maxit = mice_maxit, nuisance_spec = spec, trim = trim_cutoff,
        trme_boot_B = trme_boot_B, trme_shrink_rate = 0.99)
    }
  }, error = function(e) data.frame(estimate=NA_real_, se=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_, status="failed", error=conditionMessage(e))),
  warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  valid <- all(is.finite(unlist(rr[c("estimate", "se", "ci_lower", "ci_upper")]))) && rr$se >= 0 && rr$ci_lower <= rr$ci_upper
  if (!"status" %in% names(rr)) rr$status <- if (valid) "ok" else "failed"
  if (!"error" %in% names(rr)) rr$error <- if (valid) NA_character_ else paste(unique(warnings), collapse=" | ")
  rr$warnings <- paste(unique(warnings), collapse=" | ")
  rr$data_seed <- base_seed + sid * 100000L + sim
  rr$fit_seed <- seeds[["fit"]]; rr$internal_seed <- seeds[["internal"]]
  rr$drcmd_k <- if (method == "drcmd") 10L else NA_integer_
  rr$mice_rf_trees <- if (uses_rf_mice(family,method,spec)) execution_settings$mice_rf_trees else NA_integer_
  rr
}
