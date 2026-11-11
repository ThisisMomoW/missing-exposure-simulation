analysis_covariates <- c("X1", "X2", "X3")

empty_estimate <- function() c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)

capture_analysis <- function(expr) {
  warnings <- character(); error <- NA_character_
  result <- withCallingHandlers(tryCatch(expr, error = function(e) {
    error <<- conditionMessage(e); NULL
  }), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  list(result = result, warnings = warnings, error = error)
}

fit_aipw <- function(data, specification) {
  if (anyNA(data[c("A", "Y", analysis_covariates)])) stop("AIPW requires completed analysis data")
  rng <- .Random.seed
  fit <- AIPW::AIPW$new(Y = data$Y, A = data$A, W = data[analysis_covariates],
    Q.SL.library = outcome_learner(specification), g.SL.library = "SL.glm",
    k_split = 10L, verbose = FALSE)$fit()
  invisible(capture.output(fit$summary(g.bound = .025)))
  columns <- c("Estimate", "SE", "95% LCL", "95% UCL")
  row <- which(trimws(rownames(fit$result)) == "Risk Difference")
  if (length(row) != 1L || !all(columns %in% colnames(fit$result))) stop("Unexpected AIPW result schema")
  value <- setNames(as.numeric(fit$result[row, columns]), names(empty_estimate()))
  if (any(!is.finite(value)) || value['se'] <= 0 || value['ci_lower'] > value['ci_upper'])
    stop("AIPW returned invalid uncertainty estimates")
  list(value = value, diagnostics = list(n = nrow(data), fit_rng = rng,
    fold_indices = fit$libs$validation_index, outer_folds = 10L, internal_folds = 9L,
    propensity_range = range(fit$obs_est$raw_p_score)))
}

imputation_predictors <- function(data, specification, include_y) {
  x <- data[analysis_covariates]
  if (include_y) x$Y <- data$Y
  if (specification == "interaction") {
    if (!include_y) stop("Interaction imputation requires Y")
    for (j in 1:3) x[[paste0("YX", j)]] <- data$Y * data[[paste0("X", j)]]
  }
  x
}

fit_single_imputation <- function(data, specification, include_y) {
  x <- imputation_predictors(data, specification, include_y)
  observed <- !is.na(data$A)
  if (specification == "rf") {
    fit <- SuperLearner::SuperLearner(Y = data$A[observed], X = x[observed, ],
      family = stats::binomial(), SL.library = "SL.study_rf",
      cvControl = list(V = 10L, stratifyCV = FALSE, shuffle = TRUE))
    probability <- as.numeric(stats::predict(fit, newdata = x, onlySL = TRUE)$pred)
  } else {
    fit <- stats::glm(A ~ ., data = data.frame(A = data$A[observed], x[observed, ]),
                      family = stats::binomial())
    probability <- as.numeric(stats::predict(fit, newdata = x, type = "response"))
  }
  probability <- pmax(1e-6, pmin(1-1e-6, probability))
  data$A[!observed] <- stats::rbinom(sum(!observed), 1L, probability[!observed])
  result <- fit_aipw(data, specification)
  result$diagnostics$imputation_predictors <- names(x)
  result$diagnostics$imputation_probability_range <- range(probability)
  result
}

fit_mr <- function(data, specification, include_y) {
  q <- outcome_learner(specification)
  # Only the outcome model contains A and can receive the common AX interaction.
  observation <- if (include_y) "SL.glm" else "SL.study_without_y"
  pseudo <- if (!include_y) "SL.study_without_y" else if (specification == "rf") "SL.study_rf" else "SL.glm"
  fit <- drcmd::drcmd(Y = data$Y, A = data$A, X = data[analysis_covariates],
    m_learners = q, g_learners = "SL.glm", r_learners = observation, po_learners = pseudo,
    k = 10L, cv_folds = 5L, cutoff = .025, eem_ind = FALSE, tml = FALSE,
    Rprobs = NA, quiet = TRUE, parallel = FALSE, att = FALSE, atc = FALSE)
  estimate <- fit$results$estimates$psi_hat_ate
  se <- fit$results$ses$psi_hat_ate
  if (length(estimate) != 1L || length(se) != 1L || !is.finite(estimate) || !is.finite(se) || se <= 0)
    stop("MR returned an invalid estimate or SE")
  list(value = c(estimate = estimate, se = se, ci_lower = estimate - 1.96*se, ci_upper = estimate + 1.96*se),
    diagnostics = list(outer_folds = 10L, internal_folds = 5L,
      observation_predictors = c(analysis_covariates, if (include_y) "Y"),
      pseudo_outcome_family = "gaussian", outcome_learner = q,
      observation_probability_range = range(fit$results$nuis$kappa_hat),
      propensity_range = range(fit$results$nuis$g_hat)))
}

# Rubin's large-sample degrees of freedom use the number of usable Q/U pairs.
rubin_pool <- function(q, u) {
  if (!is.numeric(q) || !is.numeric(u) || length(q)!=length(u)) stop("Q and U must be numeric vectors of equal length")
  valid <- is.finite(q) & is.finite(u) & u >= 0
  q <- q[valid]; u <- u[valid]; m <- length(q)
  if (m < 2L) stop("Rubin pooling needs at least two valid estimate/variance pairs")
  qbar <- mean(q); ubar <- mean(u); between <- stats::var(q)
  total <- ubar + (1+1/m)*between
  df <- if (between <= .Machine$double.eps) Inf else (m-1)*(1+ubar/((1+1/m)*between))^2
  se <- sqrt(total); critical <- if (is.finite(df)) stats::qt(.975, df) else stats::qnorm(.975)
  value <- c(estimate = qbar, se = se, ci_lower = qbar-critical*se, ci_upper = qbar+critical*se)
  if (any(!is.finite(value)) || se <= 0) stop("Invalid pooled variance or confidence interval")
  list(value = value, components = list(m = m, qbar = qbar, ubar = ubar,
    between = between, total = total, df = df))
}

analyze_imputations <- function(completed, specification,
    minimum=study_settings()$min_valid_imputations, shared_rng=NULL, analyzer=fit_aipw) {
  m <- length(completed)
  minimum <- validate_imputation_minimum(minimum,m)
  if(!is.null(shared_rng) && (length(shared_rng)!=m || any(vapply(shared_rng,is.null,logical(1)))))
    stop("Paired analysis needs one RNG state per imputation")
  q <- u <- rep(NA_real_,m); fits <- vector("list",m); rng <- vector("list",m)
  for(j in seq_len(m)) {
    if(!is.null(shared_rng)) assign(".Random.seed",shared_rng[[j]],envir=.GlobalEnv)
    rng[[j]] <- .Random.seed
    fits[[j]] <- capture_analysis({
      data <- if(is.function(completed[[j]])) completed[[j]]() else completed[[j]]
      analyzer(data,specification)
    })
    value <- fits[[j]]$result$value
    if(!is.null(value)) {
      q[j] <- unname(value["estimate"])
      se <- unname(value["se"])
      if(length(se)==1L && is.finite(se) && se>=0) u[j] <- se^2
    }
  }
  valid <- is.finite(q)&is.finite(u)&u>=0
  diagnostics <- data.frame(imputation=seq_len(m),attempted=TRUE,estimate=q,variance=u,valid=valid,
    error=vapply(fits,function(z)z$error,character(1)))
  diagnostics$warnings <- I(lapply(fits,`[[`,"warnings"))
  pool <- if(sum(valid)>=minimum) capture_analysis(rubin_pool(q,u)) else
    list(result=NULL,warnings=character(),error=paste(sum(valid),"valid imputations;",minimum,"required"))
  success <- !is.null(pool$result)
  list(value=if(success)pool$result$value else empty_estimate(),
    diagnostics=list(mice_requested=m,mice_valid=sum(valid),mice_failed=sum(!valid),
      mice_pooling_status=if(!success)"failed" else if(all(valid))"complete" else "partial",
      mice_q=q,mice_u=u,component_fits=diagnostics,component_diagnostics=lapply(fits,function(z)z$result$diagnostics),
      component_rng=rng,pooling=pool$result$components,pooling_error=pool$error,
      min_valid_imputations=minimum))
}

fit_multiple_imputation <- function(data,specification,include_y,seed,
    minimum=study_settings()$min_valid_imputations,shared_rng=NULL) {
  minimum <- validate_imputation_minimum(minimum)
  # Y is retained as an analysis response even when it is excluded as a predictor.
  d <- data[c("A","Y",analysis_covariates)]
  if(specification=="interaction")for(j in 1:3)d[[paste0("YX",j)]]<-d$Y*d[[paste0("X",j)]]
  d$A <- factor(d$A,levels=0:1)
  method <- mice::make.method(d);method[]<-"";method["A"]<-if(specification=="rf")"rf" else "logreg"
  pm <- mice::make.predictorMatrix(d);pm[,]<-0
  predictors <- c(if(include_y)"Y",analysis_covariates,if(specification=="interaction")paste0("YX",1:3))
  pm["A",predictors]<-1
  arguments <- list(data=d,m=10L,maxit=10L,method=method,predictorMatrix=pm,seed=seed,printFlag=FALSE)
  if(specification=="rf")arguments<-c(arguments,list(ntree=500L,rfPackage="ranger"))
  imputation <- capture_analysis(do.call(mice::mice,arguments))
  if(is.null(imputation$result)) {
    result <- list(value=empty_estimate(),diagnostics=list(mice_requested=10L,mice_valid=0L,mice_failed=10L,
      mice_pooling_status="imputation_failed",mice_q=rep(NA_real_,10L),mice_u=rep(NA_real_,10L),
      component_fits=data.frame(imputation=1:10,attempted=FALSE,estimate=NA_real_,variance=NA_real_,valid=FALSE,
        error="Initial imputation did not complete"),imputation_error=imputation$error,
      min_valid_imputations=minimum))
    result$diagnostics$component_fits$warnings <- I(rep(list(character()),10L))
  } else {
    imp <- imputation$result
    completed <- lapply(1:10,function(j) local({
      index <- j
      function() { z <- mice::complete(imp,index); z$A<-as.integer(as.character(z$A)); z }
    }))
    result <- analyze_imputations(completed,specification,minimum,shared_rng)
    result$diagnostics$logged_events <- imp$loggedEvents
    result$diagnostics$chain_mean <- imp$chainMean
    result$diagnostics$chain_variance <- imp$chainVar
  }
  result$diagnostics$imputation_warnings <- imputation$warnings
  result$diagnostics$predictor_matrix <- pm
  result$diagnostics$mice_method <- method
  result$diagnostics$mice_rf_ntree <- if(specification=="rf")500L else NA_integer_
  result
}
