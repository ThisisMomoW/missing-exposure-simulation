arg0 <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
if (!file.exists(file.path(dirname(arg0), "../renv/setup-complete"))) stop("Run Rscript scripts/setup.R first.")
source(file.path(dirname(arg0), "../R/project.R"))
a<-commandArgs(TRUE);format<-if(length(a))a[1] else "html"
if(!format %in% c("html","pdf"))stop("Choose html or pdf")
run<-if(length(a)>=2L)a[2] else "simulation"
if(!grepl("^[a-zA-Z0-9_-]+$",run)) stop("Invalid run name")
out<-file.path(root,"outputs",run)
if(!file.exists(file.path(out,"Results/settings.rds")))stop("Aggregate the complete run first")
settings<-readRDS(file.path(out,"Results/settings.rds"))
if(!identical(settings$run,run) || length(settings$mice_rf_trees)!=1L) stop("Missing or incompatible RF-MICE result settings")
if(!setequal(settings$families,families)) stop("Report rendering requires all six scenario families")
rmarkdown::render(file.path(root,"report/report.Rmd"),output_format=paste0(format,"_document"),
 output_file=paste0("report.",format),output_dir=out,
 intermediates_dir=out,knit_root_dir=root,
 params=list(project_root=root,results_dir=file.path("outputs",run,"Results")),envir=new.env(parent=globalenv()),quiet=TRUE)
