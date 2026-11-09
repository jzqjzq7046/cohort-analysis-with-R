.f4_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.f4_root,'R/family_four_stage.R'),encoding='UTF-8')
if(sys.nframe()==0L)f4_run('discovery',commandArgs(TRUE),.f4_root)else if(isTRUE(getOption('fatty_liver.autorun',TRUE)))f4_run('discovery',root=.f4_root)
