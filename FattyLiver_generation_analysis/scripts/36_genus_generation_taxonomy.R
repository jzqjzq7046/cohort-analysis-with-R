# Independently runnable: genus aggregation, family-cluster GEE and taxonomy map.
.gg_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.gg_root,'R/genus_generation_taxonomy.R'),encoding='UTF-8')
if(sys.nframe()==0L)gg_main(commandArgs(TRUE),.gg_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))gg_main(root=.gg_root)
