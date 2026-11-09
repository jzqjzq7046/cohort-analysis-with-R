.f6_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.f6_root,'R/fig6_family_LMM.R'),encoding='UTF-8')
if(sys.nframe()==0L)f6_main(commandArgs(TRUE),.f6_root)else if(isTRUE(getOption('fatty_liver.autorun',TRUE)))f6_main(root=.f6_root)
