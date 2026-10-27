arg0 <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
if (!file.exists(file.path(dirname(arg0), "../renv/setup-complete"))) stop("Run Rscript scripts/setup.R first.")
source(file.path(dirname(arg0), "../R/project.R"))
a<-commandArgs(TRUE);run<-if(length(a))a[1] else "simulation"
if(!grepl("^[a-zA-Z0-9_-]+$",run))stop("Invalid run")
out<-file.path(root,"outputs",run);dir.create(file.path(out,"figures"),recursive=TRUE,showWarnings=FALSE)
settings_path<-file.path(out,"Results/settings.rds")
if(!file.exists(settings_path)) stop("Aggregate the named run first")
settings<-readRDS(settings_path)
if(length(settings$mice_rf_trees)!=1L || !identical(settings$run,run)) stop("Missing or incompatible RF-MICE result settings")
if(run==bundled_run) {
  # Knit the actual report chunks to export exactly its figures without fitting.
  params<-list(project_root=root,results_dir=file.path("outputs",run,"Results"),table_n=1000)
  knitr::opts_knit$set(root.dir=root)
  knitr::opts_chunk$set(fig.path=paste0(file.path(out,"figures"),"/"))
  knitr::knit(file.path(root,"report/report.Rmd"),output=file.path(out,"report.md"),envir=environment(),quiet=TRUE)
} else {
  for(family in families) {
    path<-file.path(out,"tables",paste0(family,".csv"));if(!file.exists(path))next
    x<-read.csv(path)
    for(metric in c("bias","rmse","coverage")) {
      p<-ggplot2::ggplot(x,ggplot2::aes(x=method,y=.data[[metric]],colour=nuisance_spec))+
        ggplot2::geom_point(position=ggplot2::position_dodge(width=.4))+
        ggplot2::facet_wrap(~n)+ggplot2::coord_flip()+ggplot2::theme_bw()+
        ggplot2::labs(subtitle=if(any(uses_rf_mice(family,x$method,x$nuisance_spec)))
          paste("RF-MICE imputation:",settings$mice_rf_trees,"trees") else NULL)
      ggplot2::ggsave(file.path(out,"figures",paste0(family,"_",metric,".png")),p,width=9,height=5,dpi=140)
    }
  }
}
