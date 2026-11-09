# Current standalone pipeline. Source in RStudio or invoke with project Rscript.
.v3_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.v3_root,'R/reanalysis_20261024.R'),encoding='UTF-8')
if(sys.nframe()==0L)v3_main(commandArgs(TRUE),root=.v3_root)else if(isTRUE(getOption('fatty_liver.autorun',TRUE)))v3_main(root=.v3_root)
