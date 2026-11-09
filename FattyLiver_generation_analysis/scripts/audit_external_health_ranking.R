invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
root<-'E:/YFS/workspace/FattyLiver_generation_analysis'
source(file.path(root,'R/external_health_ranking.R'))
source(file.path(root,'R/literature_ABC_20261028.R'))
cfg<-lr_setup(root)
out<-lr_adopted(root,'30_External_health_sets')
r<-readRDS(file.path(out,'private/analysis_result.rds'))
stopifnot(identical(r$protocol,ehr_protocol))
lr_check_hashes(r$hashes)
audit<-list()
for(model in unique(r$results$model))for(score in ehr_scores){
 f<-readRDS(file.path(out,'private',paste0(model,'_',score,'_fit.rds')))
 if(inherits(f,'geeglm')){
  X<-model.matrix(f);res<-residuals(f,type='response');B<-solve(crossprod(X))
  U<-rowsum(X*as.numeric(res),f$id,reorder=FALSE)
  V<-B%*%crossprod(U)%*%B
  a<-sqrt(V['exposure','exposure']);b<-sqrt(vcov(f)['exposure','exposure'])
  stopifnot(abs(a-b)<1e-8)
  audit[[paste(model,score)]]<-data.frame(model=model,score=score,manual_cluster_SE=a,geepack_SE=b,absolute_error=abs(a-b),passed=TRUE)
 }
}
fl_write_csv(ov_bind(audit),file.path(out,'audit/independent_cluster_sandwich_check.csv'))
# Snapshot exactly the functions used by the renderer, so the static checker
# sees fonts and devices defined in shared helpers instead of false missing-device failures.
functions<-c('nf_theme','tf_plot_theme','lr_theme','tf_export','lr_export','ehr_figures')
lines<-unlist(lapply(functions,function(n)c(paste0('# Runtime rendering dependency: ',n),paste0(n,' <- '),deparse(get(n)),'')))
writeLines(lines,file.path(out,'audit/plotting_source_bundle.R'))
stopifnot(all(vapply(list.files(file.path(out,'figures'),full.names=TRUE),function(p)file.info(p)$size>1000,logical(1))))
cat('Independent family-sandwich checks:',length(audit),'PASS\n')
cat('Original/reference hashes: PASS\n')
