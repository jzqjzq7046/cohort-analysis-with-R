# B: official MaAsLin2 ordinal and categorical ultrasound-grade models.
# Independently runnable; no A/C execution. Source: Mei et al. Nat Med 2024.
.lr_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.lr_root,'R/literature_ABC_20261028.R'),encoding='UTF-8')
if(sys.nframe()==0L)lr_main('B',commandArgs(TRUE),.lr_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))lr_main('B',root=.lr_root)
