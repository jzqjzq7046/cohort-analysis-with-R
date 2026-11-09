# Numbered pipeline for the user-authorized A/B/C follow-up only.
.abc_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
if(getRversion()<'4.5'){
  .abc_status<-system2(file.path(.abc_root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(.abc_root,'scripts/31_ABC_pipeline.R'))))
  if(.abc_status!=0L)stop('ABC pipeline failed')
}else{
  source(file.path(.abc_root,'R/abc_analysis.R'),encoding='UTF-8')
  for(.abc_which in c('A','B','C'))abc_main(.abc_which,if(sys.nframe()==0L)commandArgs(TRUE)else character(),.abc_root)
}
