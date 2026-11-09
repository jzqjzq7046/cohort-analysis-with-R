# Independent D entry: observed HUMAnN community pathways, original local inputs.
.ade_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
source(file.path(.ade_root,'R/generation_ADE_extensions.R'),encoding='UTF-8')
if(sys.nframe()==0L)ade_main('D',commandArgs(TRUE),.ade_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))ade_main('D',root=.ade_root)
