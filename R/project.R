# Scripts locate the project from their own path, regardless of the shell's cwd.
project_root <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg) != 1L) stop("Run this entry point with Rscript.")
  normalizePath(file.path(dirname(sub("^--file=", "", arg)), ".."), mustWork = TRUE)
}
root <- project_root()
Sys.setenv(RENV_PATHS_LIBRARY=file.path(root,"renv/library"),
           RENV_CONFIG_SANDBOX_ENABLED="FALSE", RENV_PROJECT=root, RENV_PATHS_CACHE=file.path(root,"renv/cache"),
           RENV_PATHS_ROOT=file.path(root,"renv/local"), RENV_CONFIG_AUTO_SNAPSHOT="FALSE")
if (file.exists(file.path(root, "renv/activate.R"))) {
  source(file.path(root, "renv/activate.R"))
  renv::load(root, quiet = TRUE)
} else stop("Run Rscript scripts/setup.R first.")
options(stringsAsFactors = FALSE)
families <- c("missing_share", "gamma_y", "outcome_interaction",
              "missing_mechanism", "outcome_interaction_mcar", "nonlinear_binary")
methods <- c("full_aipw", "cc_aipw", "single_imp_aipw", "mice_aipw", "drcmd", "tr_aipw")
bundled_run <- "bundled_rf10"
bundled_input_dir <- file.path(root,"report/inputs/mice_rf_10")
uses_rf_mice <- function(family, method, spec) {
  method == "mice_aipw" & (family == "nonlinear_binary" |
    (family == "outcome_interaction" & spec == "rf"))
}
bundled_settings <- function() {
  meta <- read.csv(file.path(bundled_input_dir,"metadata.csv"))
  stopifnot(nrow(meta)==1L,meta$dataset==bundled_run,meta$mice_rf_trees==10L)
  list(run=bundled_run,mice_rf_trees=as.integer(meta$mice_rf_trees),
       s6_drcmd_folds=as.character(meta$s6_drcmd_folds))
}
load_config <- function(family) {
  stopifnot(family %in% families)
  e <- new.env(parent = baseenv())
  sys.source(file.path(root, "config", paste0(family, ".R")), e)
  e
}
arms <- function(family) {
  if (family == "nonlinear_binary") {
    return(rbind(expand.grid(method = methods[1:5], nuisance_spec = c("glm", "rf")),
                 data.frame(method = "tr_aipw", nuisance_spec = "glm")))
  }
  x <- expand.grid(method = methods, nuisance_spec = c("main", "interaction"))
  if (family == "outcome_interaction") x <- rbind(x, data.frame(method = methods[1:5], nuisance_spec = "rf"))
  x
}
code_fingerprint <- function() {
  paths <- sort(c(list.files(file.path(root, "R"), full.names = TRUE, pattern = "[.]R$"),
                  list.files(file.path(root, "config"), full.names = TRUE, pattern = "[.]R$"),
                  file.path(root, "renv.lock")))
  tmp <- tempfile(); on.exit(unlink(tmp))
  writeLines(paste(basename(paths), tools::md5sum(paths)), tmp)
  unname(tools::md5sum(tmp))
}
# Atomic publication plus a run manifest prevents partial writes from looking complete.
atomic_save <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".pending_", tmpdir = dirname(path)); on.exit(unlink(tmp))
  saveRDS(x, tmp, compress = "xz")
  if (!file.rename(tmp, path)) stop("Cannot publish ", path)
}
