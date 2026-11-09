.abc_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.abc_root,'R/abc_analysis.R'),encoding='UTF-8')
if(sys.nframe()==0L)abc_main('B',commandArgs(TRUE),.abc_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))abc_main('B',root=.abc_root)
