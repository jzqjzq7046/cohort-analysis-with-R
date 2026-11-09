# Independent entry. Never launches module 25 or a cross-generation disease model.
.tf_root <- Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.tf_root,'R/taxa_followup.R'),encoding='UTF-8')
if(sys.nframe()==0L) tf_main('taxon_family_sharing',commandArgs(TRUE),.tf_root) else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))tf_main('taxon_family_sharing',root=.tf_root)
