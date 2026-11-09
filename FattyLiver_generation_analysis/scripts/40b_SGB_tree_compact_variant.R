# Add a display-only variant beside Figure40. Never refit or replace originals.
.scv_root <- Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.scv_root,'R/sgb_tree_compact_variant.R'),encoding='UTF-8')
if(sys.nframe()==0L)scv_main(commandArgs(TRUE),.scv_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))scv_main(root=.scv_root)
