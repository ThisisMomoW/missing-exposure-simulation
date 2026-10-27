arg0 <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
if (!file.exists(file.path(dirname(arg0), "../renv/setup-complete"))) stop("Run Rscript scripts/setup.R first.")
source(file.path(dirname(arg0), "../R/project.R"))
source(file.path(root, "R/workflow.R"))
a <- commandArgs(TRUE)
if (!length(a) || !a[1] %in% c("scenario", "all")) stop("Usage: run.R scenario FAMILY | all [--cell=N --rep=N[,N...] --method=M --spec=S --run=NAME]")
opt <- function(key, default=NULL) { z <- a[startsWith(a,paste0("--",key,"="))]; if(length(z)) sub(paste0("^--",key,"="),"",z[1]) else default }
mode <- a[1]
selected <- if (mode == "scenario") a[2] else if (!is.null(opt("family"))) opt("family") else families
stopifnot(all(selected %in% families))
run_name <- opt("run", "simulation")
if (!grepl("^[a-zA-Z0-9_-]+$",run_name)) stop("Use a simple run name.")
if (run_name %in% c(bundled_run,"report")) stop("Choose a simulation run name; bundled_rf10 is reserved for the included dataset.")
run_dir <- file.path(root,"outputs",run_name)
dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
fingerprint <- code_fingerprint()
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
# Select replication IDs without changing the configured replication count.
rep_selection <- opt("rep")
if (!is.null(rep_selection)) {
  if (!grepl("^[0-9]+(,[0-9]+)*$", rep_selection)) stop("--rep must be an integer or comma-separated integers, e.g. --rep=1,2,3")
  rep_selection <- as.integer(strsplit(rep_selection, ",", fixed=TRUE)[[1]])
  if (anyNA(rep_selection) || anyDuplicated(rep_selection)) stop("--rep must contain distinct integer replication IDs")
}
for (family in selected) {
  cfg <- load_implementation(family)
  grid <- cfg$scenario_grid(); aa <- arms(family)
  cells <- if (!is.null(opt("cell"))) as.integer(opt("cell")) else seq_len(nrow(grid))
  sims <- if (!is.null(rep_selection)) rep_selection else seq_len(cfg$nsim)
  stopifnot(all(cells %in% seq_len(nrow(grid))),all(sims %in% seq_len(cfg$nsim)))
  if(!is.null(opt("method"))) aa<-aa[aa$method==opt("method"),,drop=FALSE]
  if(!is.null(opt("spec"))) aa<-aa[aa$nuisance_spec==opt("spec"),,drop=FALSE]
  stopifnot(nrow(aa)>0)
  family_dir <- file.path(run_dir,family)
  manifest <- list(family=family,mode=mode,cells=cells,sims=sims,arms=aa,
    config=as.list(cfg)[!vapply(as.list(cfg),is.function,logical(1))],
    settings=execution_settings,code_fingerprint=fingerprint,
    versions=sapply(c("AIPW","SuperLearner","drcmd","mice","ranger"),function(p)as.character(packageVersion(p))))
  mp<-file.path(family_dir,"manifest.rds")
  if(file.exists(mp) && !identical(readRDS(mp),manifest)) stop("Existing run has different configuration: ",family_dir,"; choose --run=NAME")
  atomic_save(manifest,mp)
  for(sid in cells) {
    sc<-grid[sid,,drop=FALSE]; constants<-cell_constants(family,sid,sc)
    for(sim in sims) {
      dat<-generate_replication(family,sid,sim,sc,constants)
      for(j in seq_len(nrow(aa))) {
        method<-as.character(aa$method[j]); spec<-as.character(aa$nuisance_spec[j])
        path<-file.path(family_dir,sprintf("cell_%03d_rep_%03d_%s_%s.rds",sid,sim,method,spec))
        if(file.exists(path)) {
          old<-readRDS(path)
          stopifnot(identical(old$code_fingerprint,fingerprint),identical(old$settings,execution_settings))
          if (uses_rf_mice(family,method,spec) &&
              !identical(old$result$mice_rf_trees,execution_settings$mice_rf_trees)) stop("RF-MICE tree count mismatch: ",path)
          next
        }
        message(family," cell ",sid," rep ",sim," ",method," / ",spec)
        rr<-fit_replication(family,sid,sim,method,spec,dat,sc)
        row<-cbind(data.frame(family=family,scenario_id=sid,sim=sim,method=method,nuisance_spec=spec,
          estimator=paste(method,spec,sep="__"),true_rd=constants$true_rd,R_intercept=constants$R_intercept),sc,rr)
        atomic_save(list(result=row,code_fingerprint=fingerprint,settings=execution_settings),path)
      }
    }
  }
}
writeLines(capture.output(sessionInfo()),file.path(run_dir,"session.txt"))
