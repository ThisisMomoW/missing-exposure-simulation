project_root <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg) != 1L) stop("Run an entry point with Rscript.")
  normalizePath(file.path(dirname(sub("^--file=", "", arg)), ".."), mustWork = TRUE)
}
root <- project_root()
Sys.setenv(RENV_PROJECT = root, RENV_PATHS_CACHE = file.path(root, "renv/cache"),
           RENV_PATHS_ROOT = file.path(root, "renv/local"),
           RENV_CONFIG_AUTO_SNAPSHOT = "FALSE", RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
if (!file.exists(file.path(root, "renv/activate.R"))) stop("Run scripts/setup.R first.")
source(file.path(root, "renv/activate.R"))
if (!requireNamespace("renv",quietly=TRUE) || !identical(normalizePath(renv::project()),root))
  stop("The project library was not activated; run scripts/setup.R")
source(file.path(root, "config/design.R"))
source(file.path(root, "config/settings.R"))
source(file.path(root, "R/dgp.R"))
source(file.path(root, "R/summaries.R"))
options(stringsAsFactors = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

atomic_rds <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temp <- tempfile(".pending_", tmpdir = dirname(file))
  on.exit(unlink(temp))
  saveRDS(object, temp, compress = "xz")
  if (!file.rename(temp, file)) stop("Could not save ", file)
}
