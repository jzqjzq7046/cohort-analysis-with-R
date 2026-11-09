# Independent user-authorized pathway redesign; original inputs remain read-only.
.pfs_root <- Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.pfs_root,'R/pathway_function_sources.R'),encoding='UTF-8')
if(sys.nframe()==0L)pfs_main(commandArgs(TRUE),.pfs_root)else
  if(isTRUE(getOption('fatty_liver.autorun',TRUE)))pfs_main(root=.pfs_root)
