# Installation is explicit and writes only inside this project.
a<-grep("^--file=",commandArgs(FALSE),value=TRUE)
root<-normalizePath(file.path(dirname(sub("^--file=","",a)),".."))
Sys.setenv(RENV_PATHS_LIBRARY=file.path(root,"renv/library"),
           RENV_CONFIG_SANDBOX_ENABLED="FALSE", RENV_PATHS_CACHE=file.path(root,"renv/cache"),RENV_PATHS_ROOT=file.path(root,"renv/local"),
           RENV_CONFIG_AUTO_SNAPSHOT="FALSE",RENV_CONFIG_PPM_ENABLED="FALSE",RENV_PROJECT=root)
options(repos=c(CRAN="https://cloud.r-project.org"))
source(file.path(root,"renv/activate.R"))
renv::restore(project=root,prompt=FALSE)
renv::status(project=root)

file.create(file.path(root,"renv/setup-complete"))
