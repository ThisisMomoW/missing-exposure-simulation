arg <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
source(file.path(dirname(arg),"../R/project.R"))
source(file.path(root,"R/workflow.R"))
opt <- parse_options(list(run="study",cell="all",method="all",specification="all",predictors="all",reps="1:500",workers="1",
                          `min-valid-imputations`=study_settings()$min_valid_imputations))
a <- select_analyses(opt$cell,opt$method,opt$specification,opt$predictors)
minimum <- validate_imputation_minimum(suppressWarnings(as.numeric(opt[['min-valid-imputations']])))
if (any(a$method=="multiple_imputation")) {
  # The predictor comparison uses identical AIPW fold RNG states in each pair.
  excluded <- a[a$method=="multiple_imputation" & !a$include_y,]
  for(id in excluded$analysis) {
    paired_id <- sub("without_y$","with_y",id)
    if(!paired_id %in% a$analysis) stop("Select both predictor variants for paired MICE folds: --predictors all")
  }
}
ids <- integer_range(opt$reps)
workers <- suppressWarnings(as.integer(opt$workers))
if (is.na(workers)||workers<1L) stop("workers must be a positive integer")
if (.Platform$OS.type=="windows"&&workers>1L) stop("Use --workers 1 on Windows")
folder <- run_directory(opt$run)
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
manifest <- list(analyses=a,design=study_design()[study_design()$cell %in% a$cell,],replications=ids,
                 settings=study_settings(),software=software_versions(),signature=implementation_signature())
manifest$settings$min_valid_imputations <- minimum
manifest$signature <- digest::digest(list(manifest$signature,minimum),algo="sha256")
path <- file.path(folder,"manifest.rds")
if(file.exists(path)) {
  if(!identical(readRDS(path),manifest)) stop("This run name has different settings or selection; use a new name")
} else atomic_rds(manifest,path)
load_estimators()
for(cell_id in unique(a$cell)) {
  cell <- study_design();cell <- cell[cell$cell==cell_id,]
  prepared <- prepare_cell(cell)
  atomic_rds(prepared,file.path(folder,"design",paste0(cell_id,".rds")))
  work <- data.frame(replication=ids)
  execute <- function(j) {
    rep <- work$replication[j]
    aa <- a[a$cell==cell_id,];aa<-aa[order(!aa$include_y),]
    paired_rng <- NULL
    for(k in seq_len(nrow(aa))) {
    analysis <- aa[k,]
    path <- record_path(folder,analysis$analysis,rep)
    if(file.exists(path)) {
      object <- readRDS(path)
      if(!identical(object$signature,manifest$signature)) stop("Checkpoint configuration differs")
      if(analysis$method=="multiple_imputation" && analysis$include_y)
        paired_rng <- object$record$diagnostics[[1]]$component_rng
      next
    }
    # A directory claim prevents two workers from fitting the same replication.
    claim <- paste0(path,".claim")
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    if(!dir.create(claim,showWarnings=FALSE)) stop("Replication is claimed by a worker: ",claim)
    record <- tryCatch({
      value <- fit_replication(prepared,analysis,rep,
        paired_rng=if(analysis$include_y)NULL else paired_rng,minimum=minimum)
      atomic_rds(list(signature=manifest$signature,record=value),path)
      value
    },finally=unlink(claim,recursive=TRUE))
    if(analysis$method=="multiple_imputation" && analysis$include_y)
      paired_rng <- record$diagnostics[[1]]$component_rng
    message(cell_id," / ",analysis$method," / ",analysis$specification," / ",rep,": ",record$status)
    }
    invisible(NULL)
  }
  if(workers==1L) lapply(seq_len(nrow(work)),execute) else {
    results <- parallel::mclapply(seq_len(nrow(work)),execute,mc.cores=workers,mc.set.seed=FALSE,mc.preschedule=FALSE)
    if(any(vapply(results,inherits,logical(1),"try-error"))) stop("One or more workers failed outside estimation")
  }
}
message("Run saved in ",folder)
