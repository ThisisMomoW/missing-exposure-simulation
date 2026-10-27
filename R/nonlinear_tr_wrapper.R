# Fit TR-AIPW to observed Scenario 6 variables and collect bootstrap diagnostics.
fit_tr_binary <- function(dat, scenario, seed) {
  rr <- empty_tr()
  rr <- tryCatch({
    d <- dat
    if (anyNA(d$Y) || !all(d$Y %in% c(0, 1))) stop("Y is not binary.")
    # Do not pass A_obs, latent Z, oracle probabilities, or extra NA columns:
    # the estimating-equation solver uses na.omit(data) internally.
    d <- data.frame(A = d$A_obs, Y = d$Y, d[c("X1", "X2", "X3")])
    fit <- trme_rd(covs = c("X1", "X2", "X3"), Y = "Y", A = "A",
      data = d, imp_model = TRUE, shrink_rate = tr_settings$shrink_rate,
      ci_alpha = 0.95, method = "TR-AIPW", bootstrap = TRUE,
      B = tr_settings$B, seed = seed, nuisance_spec = "main",
      missing_mechanism = as.character(scenario$missing_mechanism))
    boot <- as.numeric(fit$boot_est)
    boot <- boot[is.finite(boot)]
    rr$tau1 <- as.numeric(fit$tau1)
    rr$tau0 <- as.numeric(fit$tau0)
    rr$estimate <- rr$tau1 - rr$tau0
    if (!isTRUE(abs(rr$estimate - fit$Estimate) <= 1e-8)) stop("RD component mismatch.")
    rr$n_boot_valid <- length(boot)
    rr$n_boot_invalid <- tr_settings$B - length(boot)
    if (length(boot) < max(20, ceiling(0.8 * tr_settings$B))) stop("Too few valid bootstraps.")
    rr$se <- sd(boot)
    if (!is.finite(rr$se) || rr$se <= 0) stop("Invalid bootstrap SE.")
    ci <- as.numeric(quantile(boot, c(.025, .975), type = 1))
    rr$ci_lower <- rr$boot_q025 <- ci[1]
    rr$ci_upper <- rr$boot_q975 <- ci[2]
    rr$boot_median <- median(boot)
    rr$boot_min <- min(boot)
    rr$boot_max <- max(boot)
    ps <- pmin(pmax(as.numeric(fit$fit_ps), 1e-6), 1 - 1e-6)
    mw <- as.numeric(fit$miss_weights)
    rr$q01_treated_ps <- unname(quantile(ps, .01, na.rm = TRUE))
    rr$q99_treated_inv_weight <- unname(quantile(1 / ps, .99, na.rm = TRUE))
    rr$q99_control_inv_weight <- unname(quantile(1 / (1 - ps), .99, na.rm = TRUE))
    rr$q99_missingness_weight <- unname(quantile(mw, .99, na.rm = TRUE))
    rr$max_missingness_weight <- max(mw, na.rm = TRUE)
    rr$status <- "ok"
    rr
  }, error = function(e) {
    failed <- empty_tr()
    failed$error <- conditionMessage(e)
    failed
  })

  rr
}
