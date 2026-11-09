# Project-only dependency installation, explicitly authorized 2026-10-22.
# Run using the extracted tools/R-4.5.2/bin/Rscript.exe from the project root.
stopifnot(getRversion()>="4.4",file.exists("R/sgb_association.R"))
root<-normalizePath(getwd(),winslash="/")
lib<-file.path(root,".r-library-sgb")
dir.create(lib,showWarnings=FALSE,recursive=TRUE)
.libPaths(c(lib,.Library));Sys.setenv(R_LIBS_USER=lib)
options(repos=c(CRAN="https://cloud.r-project.org"),timeout=900,Ncpus=4)
cran<-c("BiocManager","dplyr","plyr","pbapply","lmerTest","lme4","optparse","logging",
  "multcomp","ggplot2","RColorBrewer","patchwork","scales","rlang","tibble","ggnewscale",
  "jsonlite","svglite")
missing<-cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,lib=lib,type="binary")
bio<-c("Maaslin2","TreeSummarizedExperiment","SummarizedExperiment","BiocGenerics")
missing<-bio[!vapply(bio,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))BiocManager::install(missing,lib=lib,version="3.22",ask=FALSE,update=FALSE,type="binary")
dest<-file.path(root,"tools/downloads/maaslin3-paper.zip")
src<-file.path(root,"tools/maaslin3-paper-source")
if(!file.exists(dest))download.file("https://api.github.com/repos/biobakery/maaslin3/zipball/0.99.16-MaAsLin-3-paper",dest,mode="wb",method="libcurl")
if(!dir.exists(src)){dir.create(src);unzip(dest,exdir=src)}
pkg<-list.dirs(src,recursive=FALSE,full.names=TRUE)
pkg<-pkg[file.exists(file.path(pkg,"DESCRIPTION"))]
stopifnot(length(pkg)==1)
if(!requireNamespace("maaslin3",quietly=TRUE)||as.character(packageVersion("maaslin3"))!="0.99.16")
  install.packages(pkg,repos=NULL,type="source",lib=lib)
required<-c(cran,bio,"maaslin3")
stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)),
  as.character(packageVersion("maaslin3"))=="0.99.16")
dir.create(file.path(root,"inventory/sgb_runtime_20261022"),showWarnings=FALSE)
ip<-as.data.frame(installed.packages(lib.loc=lib),stringsAsFactors=FALSE)
write.csv(ip[,c("Package","Version","Built","LibPath")],file.path(root,"inventory/sgb_runtime_20261022/packages.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,"inventory/sgb_runtime_20261022/sessionInfo.txt"))
write.csv(data.frame(path=c(dest,file.path(pkg,"DESCRIPTION")),md5=unname(tools::md5sum(c(dest,file.path(pkg,"DESCRIPTION"))))),
  file.path(root,"inventory/sgb_runtime_20261022/paper_source_manifest.csv"),row.names=FALSE)
cat("READY: project-local R and package library; original R and data unchanged.\n")
