# Independent entry: select genera by personal steatosis association, not abundance.
.sg_root <- Sys.getenv('FATTYLIVER_ROOT', unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', '.UTF-8'))
source(file.path(.sg_root, 'R/steatosis_selected_genera.R'), encoding='UTF-8')
if (sys.nframe() == 0L) sg_main(commandArgs(TRUE), .sg_root) else
  if (isTRUE(getOption('fatty_liver.autorun', TRUE))) sg_main(root=.sg_root)
