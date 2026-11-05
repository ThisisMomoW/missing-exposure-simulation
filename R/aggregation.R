validate_rows <- function(x, family, cells, sims, expected_arms) {
  keys <- c("scenario_id", "sim", "method", "nuisance_spec")
  if (anyDuplicated(x[keys])) stop("Duplicate replication keys: ", family)
  cfg <- load_config(family); grid <- cfg$scenario_grid()
  expected <- merge(merge(data.frame(scenario_id=cells),data.frame(sim=sims),by=NULL),expected_arms,by=NULL)
  key <- function(d) do.call(paste,c(d[keys],sep="/"))
  if(!setequal(key(x),key(expected))) stop("Incomplete or unexpected replication grid: ",family)
  for (sid in cells) {
    rows <- x[x$scenario_id == sid,,drop=FALSE]
    for(nm in names(grid)) if(any(is.na(rows[[nm]])) || any(rows[[nm]] != grid[[nm]][sid])) stop("Design mismatch: ",family," / ",nm)
    if(length(unique(rows$true_rd))!=1L || !all(is.finite(rows$true_rd))) stop("Inconsistent truth: ",family)
    if(length(unique(rows$R_intercept))!=1L) stop("Inconsistent missingness calibration: ",family)
  }
  invisible(TRUE)
}
validate_rf_mice_settings <- function(x, family, ntree) {
  selected <- uses_rf_mice(family,x$method,x$nuisance_spec)
  if(any(selected) && (!"mice_rf_trees" %in% names(x) ||
      anyNA(x$mice_rf_trees[selected]) || any(x$mice_rf_trees[selected]!=ntree))) {
    stop("Missing or incompatible RF-MICE tree counts: ",family)
  }
  if (any(!is.na(x$mice_rf_trees[!selected])))
    stop("RF-MICE tree count assigned to a different method: ",family)
  expected <- imputation_model(family,x$method,x$nuisance_spec)
  if (!"imputation" %in% names(x) || !identical(as.character(x$imputation),expected))
    stop("Incompatible imputation specification: ",family)
  mi <- x$method=="mice_aipw"
  expected_method <- ifelse(selected,"rf","logreg")
  if (any(mi) && (!all(c("mice_method","mice_m","mice_maxit") %in% names(x)) ||
      anyNA(x$mice_method[mi]) || any(x$mice_method[mi]!=expected_method[mi]) ||
      anyNA(x$mice_m[mi]) || any(x$mice_m[mi]!=10L) ||
      anyNA(x$mice_maxit[mi]) || any(x$mice_maxit[mi]!=10L)))
    stop("Incompatible MICE settings: ",family)
  invisible(TRUE)
}
validate_mice_intervals <- function(x, family, s6_method) {
  if(length(s6_method)!=1L || is.na(s6_method) || s6_method != "rubin_t")
    stop("Missing or unsupported Scenario 6 MICE interval method")
  selected <- x$method=="mice_aipw"
  expected <- mice_interval_method(family,x$method,s6_method)
  if(any(selected) && (!"mice_ci_method" %in% names(x) ||
      anyNA(x$mice_ci_method[selected]) || any(x$mice_ci_method[selected]!=expected[selected])))
    stop("Missing or incompatible MICE interval methods: ",family)
  invisible(TRUE)
}
valid_rows <- function(x) {
  ok <- is.finite(x$estimate) & is.finite(x$se) & x$se >= 0 &
    is.finite(x$ci_lower) & is.finite(x$ci_upper) & x$ci_lower <= x$ci_upper
  if("status" %in% names(x)) ok <- ok & !is.na(x$status) & x$status=="ok"
  ok
}
# Report joint validity and coverage denominators alongside operating characteristics.
failure_audit <- function(x) {
  x$valid <- valid_rows(x)
  x$covered <- x$ci_lower <= x$true_rd & x$true_rd <= x$ci_upper
  x |>
    dplyr::group_by(scenario_id,method,nuisance_spec) |>
    dplyr::summarise(n_attempted=dplyr::n(),n_estimate=sum(is.finite(estimate)),
      n_se=sum(is.finite(se)),n_ci=sum(is.finite(ci_lower)&is.finite(ci_upper)),
      n_valid=sum(valid),n_failed=sum(!valid),failure_rate=mean(!valid),
      coverage_valid=if(any(valid)) mean(covered[valid]) else NA_real_,
      coverage_all=sum(covered & valid,na.rm=TRUE)/dplyr::n(),.groups="drop")
}
report_summary <- function(x, family) {
  if(family!="nonlinear_binary") return(summarise_oc(x,group_vars=load_config(family)$group_vars))
  # The report requires status=ok, finite estimates/SE/CI, SE >= 0, and ordered CI.
  counts <- x |> dplyr::group_by(scenario_id,n,gamma_int,method,nuisance_spec) |>
    dplyr::summarise(n_present=dplyr::n(),n_valid=sum(valid_rows(dplyr::pick(dplyr::everything()))),.groups="drop")
  z <- x[valid_rows(x),]
  z |> dplyr::group_by(scenario_id,method,nuisance_spec) |>
    dplyr::summarise(true_rd=dplyr::first(true_rd),mean_estimate=mean(estimate),
      bias=mean(estimate-true_rd),mc_sd=sd(estimate),mean_se=mean(se),
      rmse=sqrt(mean((estimate-true_rd)^2)),
      coverage=mean(ci_lower<=true_rd & true_rd<=ci_upper),.groups="drop") |>
    dplyr::right_join(counts,by=c("scenario_id","method","nuisance_spec")) |>
    dplyr::mutate(se_ratio=ifelse(mc_sd>0,mean_se/mc_sd,NA_real_),n_invalid=n_present-n_valid)
}
