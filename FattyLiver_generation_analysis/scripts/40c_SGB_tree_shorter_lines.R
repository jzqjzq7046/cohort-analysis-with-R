# Reuse the approved compact layout; halve its added coloured terminal guides.
.scs_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.scs_root,'R/sgb_tree_compact_variant.R'),encoding='UTF-8')
if(sys.nframe()==0L)scv_short_main(commandArgs(TRUE),.scs_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))scv_short_main(root=.scs_root)
