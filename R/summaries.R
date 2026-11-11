# Validity is recorded by the analysis; failed attempts remain in the denominator.
operating_characteristics <- function(records) {
  required <- c("analysis", "replication", "estimate", "se", "ci_lower", "ci_upper", "truth", "valid")
  if (!all(required %in% names(records))) stop("Missing replication fields")
  if (anyDuplicated(records[c("analysis", "replication")])) stop("Duplicate replication key")
  if (!is.logical(records$valid) || anyNA(records$valid)) stop("Each attempt needs a validity flag")
  finite <- is.finite(records$estimate) & is.finite(records$se) & records$se >= 0 &
    is.finite(records$ci_lower) & is.finite(records$ci_upper) & records$ci_lower <= records$ci_upper
  if (any(records$valid & !finite)) stop("A valid record has invalid uncertainty fields")
  rows <- lapply(split(records, records$analysis), function(x) {
    truth <- unique(x$truth)
    if (length(truth) != 1L || !is.finite(truth)) stop("Incompatible population ATEs")
    z <- x[x$valid, ]
    attempts <- nrow(x); successes <- nrow(z)
    err <- z$estimate - truth
    coverage <- if (successes) mean(z$ci_lower <= truth & truth <= z$ci_upper) else NA_real_
    pooling <- if ("mice_pooling_status" %in% names(x)) x$mice_pooling_status else rep(NA_character_,attempts)
    data.frame(analysis = x$analysis[1], attempted = attempts, valid = successes,
      failed = attempts - successes, failure_rate = (attempts - successes) / attempts,
      mice_diagnostics_available = sum(!is.na(pooling)),
      mice_complete = sum(pooling=="complete",na.rm=TRUE),
      mice_partial = sum(pooling=="partial",na.rm=TRUE),
      mice_pooling_failed = sum(pooling=="failed",na.rm=TRUE),
      mice_imputation_failed = sum(pooling=="imputation_failed",na.rm=TRUE),
      truth = truth, bias = if (successes) mean(err) else NA_real_,
      relative_bias = if (successes && truth != 0) 100 * mean(err) / truth else NA_real_,
      rmse = if (successes) sqrt(mean(err^2)) else NA_real_,
      coverage = coverage, mc_sd = if (successes > 1L) stats::sd(z$estimate) else NA_real_,
      mean_se = if (successes) mean(z$se) else NA_real_,
      coverage_mcse = if (successes) sqrt(coverage * (1 - coverage) / successes) else NA_real_,
      coverage_all = sum(z$ci_lower <= truth & truth <= z$ci_upper) / attempts)
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

figure_summaries <- function(summary, number, analyses = study_analyses()) {
  selection <- figure_design(number, analyses)
  if (anyDuplicated(summary$analysis)) stop("Duplicate analysis summaries")
  index <- match(selection$analysis, summary$analysis)
  if (anyNA(index)) stop("Figure ", number, " lacks results for ", sum(is.na(index)), " analyses")
  cbind(selection, summary[index, setdiff(names(summary), "analysis"), drop = FALSE])
}
