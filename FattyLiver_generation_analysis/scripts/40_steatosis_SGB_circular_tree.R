# Saved SGB association results and existing reference phylogeny only.
.sst_root <- Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.sst_root,'R/significant_sgb_circular_tree.R'),encoding='UTF-8')
if(sys.nframe()==0L)sst_main(commandArgs(TRUE),.sst_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))sst_main(root=.sst_root)
