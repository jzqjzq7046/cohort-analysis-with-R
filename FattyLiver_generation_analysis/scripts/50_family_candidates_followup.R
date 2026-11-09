.fc_root <- Sys.getenv('FATTYLIVER_ROOT', unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', '.UTF-8'))
source(file.path(.fc_root, 'R/family_candidate_followup.R'), encoding='UTF-8')
if (sys.nframe() == 0L) fc_main(commandArgs(TRUE), .fc_root) else if (isTRUE(getOption('fatty_liver.autorun', TRUE))) fc_main(root=.fc_root)
