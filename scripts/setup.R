args <- commandArgs(trailingOnly = FALSE)
file <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(file), ".."), mustWork = TRUE)
Sys.setenv(RENV_PROJECT = root, RENV_PATHS_CACHE = file.path(root, "renv/cache"),
           RENV_PATHS_ROOT = file.path(root, "renv/local"),
           RENV_CONFIG_AUTO_SNAPSHOT = "FALSE", RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
source(file.path(root, "renv/activate.R"))
renv::restore(project = root, prompt = FALSE)
message("Dependencies restored in the project library.")
