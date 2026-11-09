# Current C1/C2/B1/B2/P1/P2 entry, 2026-10-24.
# Original source retained at code_history/before_two_generation_20261024/scripts/00_prepare_confounders.R
# Historical method-specific filenames now delegate to the adopted protocol.
.current_root <- Sys.getenv("FATTYLIVER_ROOT", unset="E:/YFS/workspace/FattyLiver_generation_analysis")
if(.Platform$OS.type=="windows")invisible(Sys.setlocale("LC_CTYPE",".UTF-8"))
source(file.path(.current_root,"R/reanalysis_20261024.R"),encoding="UTF-8")
fl_current_entry <- function(args=character()) {
  
  pointer<-file.path(.current_root,"outputs/24_Two_generation_reanalysis/LATEST_RUN.txt")
  if(file.exists(pointer)&&!any(grepl("^--resume=",args)))args<-c(args,paste0("--resume=",readLines(pointer,n=1,warn=FALSE)))
  v3_main(args,stage="prepare",root=.current_root)
}
if(sys.nframe()==0L)fl_current_entry(commandArgs(TRUE))else if(isTRUE(getOption("fatty_liver.autorun",TRUE)))fl_current_entry()
