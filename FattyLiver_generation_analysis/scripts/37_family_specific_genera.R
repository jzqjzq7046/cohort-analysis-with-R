# Independently runnable application of Valles-Colomer Figure 1g / Table S5.
.fg_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.fg_root,'R/family_specific_genera.R'),encoding='UTF-8')
if(sys.nframe()==0L)fg_main(commandArgs(TRUE),.fg_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))fg_main(root=.fg_root)
