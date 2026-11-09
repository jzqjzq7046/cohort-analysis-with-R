# Saved-result visualization only. No MaAsLin2/GEE fits or BH recomputation.
.sgm_root <- Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.sgm_root,'R/genus_main_text_figures.R'),encoding='UTF-8')
if(sys.nframe()==0L)sgm_main(commandArgs(TRUE),.sgm_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))sgm_main(root=.sgm_root)
