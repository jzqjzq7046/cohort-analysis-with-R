# Independent entry: unaffected younger members, observed older steatosis exposure.
# Does not run modules 24-26 or the deferred inflammation-phenotype analysis.
.fp_root <- Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.fp_root,'R/family_steatosis_microbiome.R'),encoding='UTF-8')
if(sys.nframe()==0L)fp_main(commandArgs(TRUE),.fp_root) else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))fp_main(root=.fp_root)
