# Development preview from saved fresh estimates only; not a published result.
invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
sys.source(file.path(root,'R/pathway_function_sources.R'),envir=.GlobalEnv)
cfg<-pfs_setup(root)
sys.source(file.path(root,'R/pathway_function_sources_details.R'),envir=.GlobalEnv)
sys.source(file.path(root,'R/pathway_function_sources_figures.R'),envir=.GlobalEnv)
args<-commandArgs(TRUE);stopifnot(length(args)%in%1:2)
out<-normalizePath(args[1],winslash='/',mustWork=TRUE)
boundary<-paste0(normalizePath(file.path(root,'outputs/44_Pathway_function_sources'),winslash='/',mustWork=TRUE),'/')
stopifnot(startsWith(out,boundary))
r<-readRDS(file.path(out,'private/analysis_result.rds'))
s<-readRDS(file.path(out,'private/source_standalone.rds'))
p<-pfs_primary(r);cols<-c('feature_id','beta_log2','P_raw','P_BH')
stopifnot(identical(s$primary_signature,digest::digest(p[order(p$feature_id),cols],algo='sha256')))
for(k in c('classification','representatives','sources'))r[[k]]<-s[[k]]
if(length(args)==2L){
  stopifnot(args[2]%in%c('--strict-main','--strict-all'))
  if(args[2]=='--strict-main')pfs_main_figure(r,cfg,out)else pfs_figures(r,cfg,out)
  cat('Saved-result figure strict QA/export PASS; not adopted\n')
  quit(status=0)
}
figure<-pfs_main_figure(r,cfg,out,preview_only=TRUE)
ggplot2::ggsave(file.path(out,'audit/main_development_preview.png'),figure,
  width=183,height=268,units='mm',dpi=200,bg='white',type='cairo')
cat('Development preview saved\n')
