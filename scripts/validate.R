# Check bundled result grids, failure accounting and summary definitions.
arg0 <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
if (!file.exists(file.path(dirname(arg0), "../renv/setup-complete"))) stop("Run Rscript scripts/setup.R first.")
source(file.path(dirname(arg0), "../R/project.R"))
source(file.path(root,"R/summary_functions.R"))
source(file.path(root,"R/aggregation.R"))
if (length(commandArgs(TRUE))) stop("Usage: Rscript scripts/validate.R")
input_settings <- bundled_settings()
checks <- list()
checksums <- read.csv(file.path(bundled_input_dir,"checksums.csv"))
stopifnot(all(unname(tools::md5sum(file.path(bundled_input_dir,checksums$file)))==checksums$md5))
inputs <- list()
for (family in families) {
  x <- readRDS(file.path(bundled_input_dir,paste0(family,".rds")))
  inputs[[family]] <- x
  cfg <- load_config(family)
  cells <- bundled_cells(family)
  sims <- seq_len(cfg$nsim)
  aa <- arms(family)
  validate_rows(x,family,cells,sims,aa)
  validate_rf_mice_settings(x,family,500L)
  validate_mice_intervals(x,family,"rubin_t")
  sm <- report_summary(x,family)
  audit <- failure_audit(x)
  stopifnot(nrow(sm)==length(cells)*nrow(aa),sum(audit$n_attempted)==nrow(x),
            all(audit$n_valid+audit$n_failed==audit$n_attempted))
  checks[[family]] <- data.frame(family=family,records=nrow(x),
    method_cells=nrow(sm),valid=sum(audit$n_valid),failed=sum(audit$n_failed),mice_rf_trees=input_settings$mice_rf_trees,
    mice_ci_method=if(family=="nonlinear_binary") input_settings$s6_mice_ci_method else "rubin_t")

  # A partial cell is not a complete run; duplicate and altered keys are errors.
  sample <- x[x$scenario_id==cells[1] & x$sim==1L,,drop=FALSE]
  stopifnot(inherits(try(validate_rows(rbind(sample,sample[1,]),family,cells[1],1L,aa),silent=TRUE),"try-error"))
  stopifnot(inherits(try(validate_rows(sample[-1,],family,cells[1],1L,aa),silent=TRUE),"try-error"))
  bad <- sample; bad$n[1] <- bad$n[1]+1L
  stopifnot(inherits(try(validate_rows(bad,family,cells[1],1L,aa),silent=TRUE),"try-error"))
}
# Known values distinguish estimate, SE and interval denominators from attempts.
fixture <- data.frame(scenario_id=1,method="full_aipw",nuisance_spec="main",estimator="full_aipw__main",n=500,true_rd=.2,
  estimate=c(.1,.3,NA),se=c(.1,.2,NA),ci_lower=c(0,.25,NA),ci_upper=c(.2,.35,NA))
s <- summarise_oc(fixture,group_vars="n"); au <- failure_audit(fixture)
stopifnot(abs(s$bias)<1e-14,abs(s$rmse-.1)<1e-14,abs(s$mc_sd-sqrt(.02))<1e-14,
  abs(s$mean_se-.15)<1e-14,s$coverage==.5,s$n_total==3,s$n_finite==2,
  abs(au$failure_rate-1/3)<1e-14,abs(au$coverage_all-1/3)<1e-14)
# Recompute every displayed manuscript cell directly from replication records.
reference <- read.csv(file.path(root,"validation/manuscript_reference.csv"))
metrics <- c("bias","rmse","coverage")
comparison <- lapply(seq_len(nrow(reference)),function(i) {
  ref <- reference[i,]; x <- inputs[[ref$family]]
  x <- x[x$scenario_id==ref$scenario_id & x$method==ref$method &
    x$nuisance_spec==ref$nuisance_spec,]
  stopifnot(nrow(x)==500L,setequal(x$sim,1:500))
  point <- (is.na(x$status)|x$status=="ok") & is.finite(x$estimate)
  interval <- point & is.finite(x$ci_lower) & is.finite(x$ci_upper) & x$ci_lower<=x$ci_upper
  values <- c(bias=mean(x$estimate[point]-x$true_rd[point]),
    rmse=sqrt(mean((x$estimate[point]-x$true_rd[point])^2)),
    coverage=mean(x$ci_lower[interval]<=x$true_rd[interval] & x$true_rd[interval]<=x$ci_upper[interval]))
  stopifnot(sum(point)==ref$valid_point,sum(interval)==ref$valid_interval)
  delta <- values-unlist(ref[metrics],use.names=FALSE)
  stopifnot(max(abs(delta))<1e-13)
  data.frame(ref[c("figure","family","scenario_id","method","nuisance_spec")],
    max_absolute_difference=max(abs(delta)))
})
# Every saved Scenario 6 MICE interval can be recalculated from its components.
x <- inputs$nonlinear_binary
x <- x[x$method=="mice_aipw",]
stopifnot(nrow(x)==3000L,all(x$mice_m_valid==10L))
tvar <- x$mice_ubar+(1+1/x$mice_m_valid)*x$mice_b
r <- (1+1/x$mice_m_valid)*x$mice_b/x$mice_ubar
df <- ifelse(x$mice_b<=.Machine$double.eps,Inf,(x$mice_m_valid-1)*(1+1/r)^2)
stopifnot(max(abs(x$mice_t-tvar))<1e-14,max(abs(x$se-sqrt(tvar)))<1e-14,
  max(abs(x$mice_df-df))<1e-10,
  max(abs(x$ci_lower-(x$estimate-qt(.975,df)*sqrt(tvar))))<1e-14,
  max(abs(x$ci_upper-(x$estimate+qt(.975,df)*sqrt(tvar))))<1e-14)
# Missing, incompatible and misassigned RF/MICE settings cannot pass validation.
z <- inputs$nonlinear_binary[1:33,]
rfmi <- which(z$method=="mice_aipw" & z$nuisance_spec=="rf")[1]
bad <- z; bad$mice_rf_trees[rfmi] <- 10L
stopifnot(inherits(try(validate_rf_mice_settings(bad,"nonlinear_binary",500L),silent=TRUE),"try-error"))
bad <- z; bad$mice_ci_method[bad$method=="mice_aipw"] <- "normal_1.96"
stopifnot(inherits(try(validate_mice_intervals(bad,"nonlinear_binary","rubin_t"),silent=TRUE),"try-error"))
bad <- z; bad$mice_rf_trees[bad$method=="mice_aipw" & bad$nuisance_spec=="glm"] <- 500L
stopifnot(inherits(try(validate_rf_mice_settings(bad,"nonlinear_binary",500L),silent=TRUE),"try-error"))
out <- file.path(root,"outputs/validation")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(dplyr::bind_rows(comparison),file.path(out,"manuscript_comparison.csv"),row.names=FALSE)
checks <- dplyr::bind_rows(checks)
write.csv(checks,file.path(out,"result_checks.csv"),row.names=FALSE)
writeLines("PASS: checksums and complete grids; configuration, duplicate and incomplete rejection; summary denominators; 105 manuscript cells; 3000 Rubin intervals.",file.path(out,"checks.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"session.txt"))
message("Validated ",sum(checks$records)," records in ",nrow(checks)," families; ",sum(checks$failed)," failed or invalid fits; bundled RF-MICE trees = ",input_settings$mice_rf_trees,". No models fitted.")
