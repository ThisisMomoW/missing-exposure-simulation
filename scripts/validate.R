# Check bundled result grids, failure accounting and summary definitions.
arg0 <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
if (!file.exists(file.path(dirname(arg0), "../renv/setup-complete"))) stop("Run Rscript scripts/setup.R first.")
source(file.path(dirname(arg0), "../R/project.R"))
source(file.path(root,"R/summary_functions.R"))
source(file.path(root,"R/aggregation.R"))
if (length(commandArgs(TRUE))) stop("Usage: Rscript scripts/validate.R")
input_settings <- bundled_settings()
checks <- list()
for (family in families) {
  x <- readRDS(file.path(bundled_input_dir,paste0(family,".rds")))
  cfg <- load_config(family)
  cells <- seq_len(nrow(cfg$scenario_grid()))
  sims <- seq_len(cfg$nsim)
  aa <- arms(family)
  validate_rows(x,family,cells,sims,aa)
  sm <- report_summary(x,family)
  audit <- failure_audit(x)
  stopifnot(nrow(sm)==length(cells)*nrow(aa),sum(audit$n_attempted)==nrow(x),
            all(audit$n_valid+audit$n_failed==audit$n_attempted))
  checks[[family]] <- data.frame(family=family,records=nrow(x),
    method_cells=nrow(sm),valid=sum(audit$n_valid),failed=sum(audit$n_failed),mice_rf_trees=input_settings$mice_rf_trees)

  # A partial cell is not a complete run; duplicate and altered keys are errors.
  sample <- x[x$scenario_id==1L & x$sim==1L,,drop=FALSE]
  stopifnot(inherits(try(validate_rows(rbind(sample,sample[1,]),family,1L,1L,aa),silent=TRUE),"try-error"))
  stopifnot(inherits(try(validate_rows(sample[-1,],family,1L,1L,aa),silent=TRUE),"try-error"))
  bad <- sample; bad$n[1] <- bad$n[1]+1L
  stopifnot(inherits(try(validate_rows(bad,family,1L,1L,aa),silent=TRUE),"try-error"))
}
# Known values distinguish estimate, SE and interval denominators from attempts.
fixture <- data.frame(scenario_id=1,method="full_aipw",nuisance_spec="main",estimator="full_aipw__main",n=500,true_rd=.2,
  estimate=c(.1,.3,NA),se=c(.1,.2,NA),ci_lower=c(0,.25,NA),ci_upper=c(.2,.35,NA))
s <- summarise_oc(fixture,group_vars="n"); au <- failure_audit(fixture)
stopifnot(abs(s$bias)<1e-14,abs(s$rmse-.1)<1e-14,abs(s$mc_sd-sqrt(.02))<1e-14,
  abs(s$mean_se-.15)<1e-14,s$coverage==.5,s$n_total==3,s$n_finite==2,
  abs(au$failure_rate-1/3)<1e-14,abs(au$coverage_all-1/3)<1e-14)
out <- file.path(root,"outputs/validation")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
checks <- dplyr::bind_rows(checks)
write.csv(checks,file.path(out,"result_checks.csv"),row.names=FALSE)
writeLines("PASS: result grids; duplicate, incomplete and design mismatch rejection; summary and failure denominators.",file.path(out,"checks.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"session.txt"))
message("Validated ",sum(checks$records)," records in ",nrow(checks)," families; ",sum(checks$failed)," failed or invalid fits; bundled RF-MICE trees = ",input_settings$mice_rf_trees,". No models fitted.")
