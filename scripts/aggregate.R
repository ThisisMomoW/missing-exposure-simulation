arg0 <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
if (!file.exists(file.path(dirname(arg0), "../renv/setup-complete"))) stop("Run Rscript scripts/setup.R first.")
source(file.path(dirname(arg0), "../R/project.R"))
source(file.path(root,"R/summary_functions.R"))
source(file.path(root,"R/aggregation.R"))
a <- commandArgs(TRUE); run <- if(length(a)) a[1] else "simulation"
if(!grepl("^[a-zA-Z0-9_-]+$",run)) stop("Invalid run name")
if(run=="report") stop("Select a simulation run, or bundled_rf10 for the included 10-tree RF-MICE dataset.")
out <- file.path(root,"outputs",run); results_dir <- file.path(out,"Results")
dir.create(results_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"tables"),showWarnings=FALSE)
audits <- list()
run_settings <- if(run==bundled_run) bundled_settings() else NULL
for(family in families) {
  if(run==bundled_run) {
    x<-readRDS(file.path(bundled_input_dir,paste0(family,".rds")))
    cells<-seq_len(nrow(load_config(family)$scenario_grid())); sims<-1:500; aa<-arms(family)
    x$mice_rf_trees<-ifelse(uses_rf_mice(family,x$method,x$nuisance_spec),run_settings$mice_rf_trees,NA_integer_)
  } else {
    mp<-file.path(out,family,"manifest.rds")
    if(!file.exists(mp)) next
    m<-readRDS(mp); cells<-m$cells; sims<-m$sims; aa<-m$arms
    ntree<-m$settings$mice_rf_trees
    if(length(ntree)!=1L || is.na(ntree) || ntree<1L) stop("Missing RF-MICE tree count: ",family)
    settings<-list(run=run,mice_rf_trees=as.integer(ntree),s6_drcmd_folds=as.character(m$settings$drcmd_outer_folds),
                   code_fingerprint=m$code_fingerprint)
    if(is.null(run_settings)) run_settings<-settings
    else if(!identical(run_settings,settings)) stop("Incompatible model settings between families")
    fs<-list.files(dirname(mp),pattern="^cell_.*[.]rds$",full.names=TRUE)
    parts<-lapply(fs,readRDS)
    if(!length(parts) || any(!vapply(parts,function(z)
      identical(z$code_fingerprint,m$code_fingerprint) && identical(z$settings,m$settings),logical(1)))) stop("Mixed configurations: ",family)
    x<-dplyr::bind_rows(lapply(parts,`[[`,"result"))
  }
  validate_rf_mice_settings(x,family,run_settings$mice_rf_trees)
  validate_rows(x,family,cells,sims,aa)
  sm<-report_summary(x,family)
  sm$mice_rf_trees<-ifelse(uses_rf_mice(family,sm$method,sm$nuisance_spec),run_settings$mice_rf_trees,NA_integer_)
  au<-failure_audit(x);au$family<-family
  au$mice_rf_trees<-ifelse(uses_rf_mice(family,au$method,au$nuisance_spec),run_settings$mice_rf_trees,NA_integer_)
  audits[[family]]<-au
  write.csv(sm,file.path(out,"tables",paste0(family,".csv")),row.names=FALSE)
  if(family=="nonlinear_binary") {
    atomic_save(x,file.path(results_dir,family,"combined_with_tr_aipw/all_results.rds"))
  } else {
    atomic_save(list(results=x,raw=sm),file.path(results_dir,paste0(family,"_combined.rds")))
    write.csv(sm,file.path(results_dir,paste0(family,"_operating_characteristics.csv")),row.names=FALSE)
  }
}
if(!length(audits)) stop("No run manifests found: ",run)
write.csv(dplyr::bind_rows(audits),file.path(out,"tables/failure_audit.csv"),row.names=FALSE)
run_settings$families<-names(audits)
atomic_save(run_settings,file.path(results_dir,"settings.rds"))
message("Validated and aggregated ",length(audits)," families; RF-MICE trees = ",run_settings$mice_rf_trees,"; output: ",out)
