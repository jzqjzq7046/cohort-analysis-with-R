# A: within-family LDM results, presented as Chen 2024 Figure 2A-C.
# Reuses completed module 29; independently runnable. No refitting is hidden.
.lr_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.lr_root,'R/literature_ABC_20261028.R'),encoding='UTF-8')
if(sys.nframe()==0L)lr_main('A',commandArgs(TRUE),.lr_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))lr_main('A',root=.lr_root)
