# Standalone rebuild of module 07 / analysis 30. Run in project R automatically.
.ehr_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.ehr_root,'R/external_health_ranking.R'),encoding='UTF-8')
if(sys.nframe()==0L)ehr_main(commandArgs(TRUE),.ehr_root)else if(isTRUE(getOption('fatty_liver.autorun',TRUE)))ehr_main(root=.ehr_root)
