arg <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
source(file.path(dirname(arg),"../R/project.R"))
source(file.path(root,"R/workflow.R"))
source(file.path(root,"R/plots.R"))
opt <- parse_options(list(run="bundled",figures="1:5"))
ids <- integer_range(opt$figures,5L)
folder <- file.path(root,"outputs/figures",opt$run)
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
if(1L %in% ids) mechanism_diagram(file.path(folder,"figure_1.pdf"))
if(any(ids>1L)) {
  x <- read_run(opt$run)
  if(length(x$missing)) stop("Complete the selected run before generating performance figures")
  s <- operating_characteristics(x$records)
  for(i in ids[ids>1L]) {
    d <- figure_summaries(s,i)
    p <- performance_figure(d,i)
    dim <- figure_dimensions(i)
    ggplot2::ggsave(file.path(folder,paste0("figure_",i,".pdf")),p$plot,width=dim[1],height=dim[2],device=grDevices::cairo_pdf)
    utils::write.csv(p$data,file.path(folder,paste0("figure_",i,"_summary.csv")),row.names=FALSE)
    built <- ggplot2::ggplot_build(p$plot)
    layer <- which(vapply(p$plot$layers,function(l)inherits(l$geom,"GeomPoint"),logical(1)))
    points <- built$data[[layer]]
    stopifnot(nrow(points)==3L*nrow(d))
    utils::write.csv(points[c("PANEL","group","x","y","colour","shape")],
      file.path(folder,paste0("figure_",i,"_coordinates.csv")),row.names=FALSE)
  }
}
message("Figures saved in ",folder)
