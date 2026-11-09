# Optional independent staging of the functional community analysis. It uses
# only outcome-blind design/filter metadata and never writes feature-fit state.
invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
sys.source(file.path(root,'R/pathway_function_sources.R'),envir=.GlobalEnv)
cfg<-pfs_setup(root)
sys.source(file.path(root,'R/pathway_function_sources_details.R'),envir=.GlobalEnv)
args<-commandArgs(TRUE);stopifnot(length(args)==1L)
out<-normalizePath(args[1],winslash='/',mustWork=TRUE)
boundary<-paste0(normalizePath(file.path(root,'outputs/44_Pathway_function_sources/_running'),winslash='/',mustWork=TRUE),'/')
stopifnot(startsWith(out,boundary))
options(fatty_liver.ade_log=file.path(out,'audit/composition_stage.log'))
x<-readRDS(file.path(out,'private/functional_input.rds'))
fits<-lapply(v3_micro_samples(x),function(sc){
  d<-v3_design(x,sc$mask,sc$fields)
  vari<-vapply(seq_len(ncol(d$relative)),function(j)length(unique(d$relative[,j]))>1L,logical(1))
  d$filter$retained<-d$filter$retained&vari
  list(filter=d$filter,sample_ids=as.character(d$metadata$sample_id))
})
result<-pfs_composition(list(fits=fits),cfg,out)
saveRDS(result,file.path(out,'private/composition_standalone.rds'))
ade_log('Independent functional composition staging COMPLETED')
