load_estimators <- function() {
  suppressPackageStartupMessages(library(SuperLearner))
  source(file.path(root,"R/learners.R"),local=.GlobalEnv)
  source(file.path(root,"R/estimators.R"),local=.GlobalEnv)
}

parse_options <- function(defaults) {
  args <- commandArgs(TRUE)
  if (length(args)%%2L) stop("Arguments must be --name value pairs")
  if (!length(args)) return(defaults)
  keys <- sub("^--","",args[seq(1,length(args),2)])
  if (anyDuplicated(keys) || any(!keys %in% names(defaults))) stop("Unknown or duplicate option")
  defaults[keys] <- as.list(args[seq(2,length(args),2)])
  defaults
}

integer_range <- function(x, upper=500L) {
  parts <- strsplit(x,",",fixed=TRUE)[[1]]
  ids <- unlist(lapply(parts,function(p) {
    r <- suppressWarnings(as.integer(strsplit(p,":",fixed=TRUE)[[1]]))
    if (anyNA(r) || !length(r) || length(r)>2L || any(r<1L | r>upper)) stop("Invalid replication range")
    if (length(r)==2L) seq.int(r[1],r[2]) else r
  }))
  if (anyDuplicated(ids)) stop("Duplicate replication IDs")
  sort(ids)
}

select_analyses <- function(cell="all", method="all", specification="all", predictors="all") {
  a <- study_analyses()
  if (cell!="all") a <- a[a$cell %in% strsplit(cell,",",fixed=TRUE)[[1]],]
  if (method!="all") a <- a[a$method %in% strsplit(method,",",fixed=TRUE)[[1]],]
  if (specification!="all") a <- a[a$specification %in% strsplit(specification,",",fixed=TRUE)[[1]],]
  if (predictors!="all") {
    if (!predictors %in% c("with_y","without_y")) stop("Choose with_y, without_y or all")
    a <- a[a$include_y==(predictors=="with_y"),]
  }
  if (!nrow(a)) stop("No analyses match the requested selection")
  a
}

software_versions <- function() {
  packages <- c("AIPW","drcmd","mice","SuperLearner","ranger","ggplot2","digest")
  sha <- utils::packageDescription("drcmd")$RemoteSha
  c(R=as.character(getRversion()),setNames(vapply(packages,function(p)as.character(utils::packageVersion(p)),""),packages),
    drcmd_commit=if(is.null(sha))NA_character_ else sha)
}

implementation_signature <- function() {
  paths <- file.path(root,c("R/dgp.R","R/estimators.R","R/learners.R","R/project.R","R/workflow.R",
                            "config/design.R","config/settings.R"))
  digest::digest(list(files=setNames(unname(tools::md5sum(paths)),basename(paths)),
    settings=study_settings(),software=software_versions()),algo="sha256")
}

prepare_cell <- function(design) {
  seeds <- design_seeds(design)
  set.seed(seeds$truth); truth <- population_ate(design)
  set.seed(seeds$calibration); calibration <- calibrate_observation(design)
  list(design=design, truth=truth, calibration=calibration, seeds=seeds)
}

run_directory <- function(name) {
  if (length(name)!=1L || !grepl("^[A-Za-z][A-Za-z0-9_-]*$",name) || name=="bundled")
    stop("Choose a run name starting with a letter; bundled is reserved")
  file.path(root,"outputs/runs",name)
}

record_path <- function(folder, analysis, replication) file.path(folder,"records",analysis,sprintf("%03d.rds",replication))

valid_estimate <- function(value) length(value)==4L && all(is.finite(value)) &&
  value["se"]>0 && value["ci_lower"]<=value["ci_upper"]

fit_replication <- function(prepared, analysis, replication, paired_rng=NULL,
    minimum=study_settings()$min_valid_imputations) {
  if(analysis$method=="multiple_imputation") minimum <- validate_imputation_minimum(minimum)
  seeds <- replication_seeds(prepared$design,analysis,replication)
  set.seed(seeds$data)
  data <- generate_data(prepared$design,prepared$calibration$intercept)
  set.seed(seeds$fit)
  start <- proc.time()[["elapsed"]]
  captured <- capture_analysis(switch(analysis$method,
    complete_data=fit_aipw(transform(data,A=A_true),analysis$specification),
    complete_case=fit_aipw(data[data$R==1L,,drop=FALSE],analysis$specification),
    single_imputation=fit_single_imputation(data,analysis$specification,analysis$include_y),
    multiple_imputation={
      if(!analysis$include_y && is.null(paired_rng)) stop("The paired Y-inclusive MICE analysis did not supply fold RNG states")
      fit_multiple_imputation(data,analysis$specification,analysis$include_y,seeds$imputation,minimum,paired_rng)
    },
    mr=fit_mr(data,analysis$specification,analysis$include_y),stop("Unknown method")))
  fit <- captured$result
  value <- if (is.null(fit)) empty_estimate() else fit$value
  valid <- valid_estimate(value)
  row <- data.frame(analysis=analysis$analysis,replication=replication,
    estimate=unname(value["estimate"]),se=unname(value["se"]),
    ci_lower=unname(value["ci_lower"]),ci_upper=unname(value["ci_upper"]),
    truth=prepared$truth$ate,valid=valid,status=if(valid)"success" else "failed",
    data_seed=seeds$data,fit_seed=seeds$fit,
    imputation_seed=if(analysis$method=="multiple_imputation")seeds$imputation else NA_integer_,
    observation_intercept=prepared$calibration$intercept,
    missing_fraction=mean(data$R==0L),elapsed_seconds=proc.time()[["elapsed"]]-start,
    error=captured$error,mice_requested=if(analysis$method=="multiple_imputation")10L else NA_integer_,
    mice_valid=NA_integer_,mice_failed=NA_integer_,mice_pooling_status=NA_character_)
  row$imputation <- if(analysis$method %in% c("single_imputation","multiple_imputation"))
    if(analysis$specification=="rf")"rf" else "logistic" else NA_character_
  row$mice_rf_ntree <- if(analysis$method=="multiple_imputation"&&analysis$specification=="rf")500L else NA_integer_
  row$mice_ubar <- row$mice_between <- row$mice_total <- row$mice_df <- NA_real_
  row$warnings <- I(list(captured$warnings))
  row$diagnostics <- I(list(if(is.null(fit))list() else fit$diagnostics))
  row$mice_q <- I(list(NULL));row$mice_u <- I(list(NULL));row$mice_components <- I(list(NULL))
  if(analysis$method=="multiple_imputation" && !is.null(fit)) {
    d <- fit$diagnostics
    row$mice_requested <- d$mice_requested;row$mice_valid <- d$mice_valid;row$mice_failed <- d$mice_failed
    row$mice_pooling_status <- d$mice_pooling_status
    row$mice_q <- I(list(d$mice_q));row$mice_u <- I(list(d$mice_u))
    row$mice_components <- I(list(d$component_fits))
    if(!is.null(d$pooling)) {
      row$mice_ubar <- d$pooling$ubar;row$mice_between <- d$pooling$between
      row$mice_total <- d$pooling$total;row$mice_df <- d$pooling$df
    }
    if(!valid)row$error <- if(!is.null(d$imputation_error))d$imputation_error else d$pooling_error
  }
  row
}

read_run <- function(name) {
  if (name=="bundled") {
    folder <- file.path(root,"data/bundled")
    manifest <- readRDS(file.path(folder,"manifest.rds"))
    checks <- utils::read.csv(file.path(folder,"checksums.csv"))
    if(!identical(unname(tools::md5sum(file.path(folder,checks$file))),checks$md5)) stop("Bundled input checksum differs")
    records <- readRDS(file.path(folder,"replications.rds"))
  } else {
    folder <- run_directory(name)
    manifest <- readRDS(file.path(folder,"manifest.rds"))
    files <- list.files(file.path(folder,"records"),pattern="[.]rds$",recursive=TRUE,full.names=TRUE)
    if (!length(files)) stop("No replication records")
    objects <- lapply(files,readRDS)
    if (!all(vapply(objects,function(z)identical(z$signature,manifest$signature),logical(1))))
      stop("Incompatible configuration or software in replication records")
    records <- do.call(rbind,lapply(objects,`[[`,"record"))
  }
  keys <- paste(records$analysis,records$replication)
  expected <- as.vector(outer(manifest$analyses$analysis,manifest$replications,paste))
  if (anyDuplicated(keys) || any(!keys %in% expected)) stop("Duplicate or unexpected replication key")
  list(records=records,manifest=manifest,missing=setdiff(expected,keys))
}
