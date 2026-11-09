# C: partial-Spearman phenotype map following Asnicar et al. Nat Med 2021.
# Adults, fixed one-per-primary-family selection, six predeclared phenotypes.
.lr_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.lr_root,'R/literature_ABC_20261028.R'),encoding='UTF-8')
if(sys.nframe()==0L)lr_main('C',commandArgs(TRUE),.lr_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))lr_main('C',root=.lr_root)
