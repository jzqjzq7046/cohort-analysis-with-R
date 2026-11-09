# Independent entry: stored reference phylogeny, generation profiles and saved P2.
# Default/source() runs this module only; no model fitting or sequence inference.
.gp_root <- Sys.getenv('FATTYLIVER_ROOT', unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if (.Platform$OS.type=='windows') invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.gp_root,'R/generation_reference_phylogeny.R'),encoding='UTF-8')
if (sys.nframe()==0L) gp_main(commandArgs(TRUE),.gp_root) else
  if (isTRUE(getOption('fatty_liver.autorun',TRUE))) gp_main(root=.gp_root)
