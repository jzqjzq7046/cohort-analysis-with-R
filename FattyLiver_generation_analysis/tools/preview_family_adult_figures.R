# Unpublished visual QA using completed adult models only. No refits, no
# fabricated all-age row, and no adoption while the final sensitivity runs.
root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(root,'R/family_four_stage.R'),encoding='UTF-8');cfg<-f4_setup(root)
paths<-list.files(file.path(root,'outputs/45_Family_de_novo_discovery/_running'),pattern='incremental_fits.rds',recursive=TRUE,full.names=TRUE)
fits<-readRDS(paths[which.max(file.info(paths)$mtime)]);x<-abc_source(root)$x
stopifnot(all(c('B2','Adult_B1_same_sample')%in%names(fits)))
e<-ov_bind(lapply(fits,function(a)a$results));s<-ov_bind(lapply(fits,function(a)a$summary))
specs<-f4_models(x)
for(k in names(fits)) { e$model_display[e$model==k]<-specs[[k]]$display;s$model_display[s$model==k]<-specs[[k]]$display }
ann<-ade_labels(x$features);e<-merge(e,ann,by='feature_id',all.x=TRUE,sort=FALSE)
r<-list(effects=e,summary=s,kind='SGB')
out<-file.path(root,'outputs/family_four_stage_20261104/qa_preview')
for(d in c('audit','figures','private'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
saveRDS(lapply(specs,function(a)a$m),file.path(out,'private/model_populations.rds'))
writeLines('UNPUBLISHED QA: completed adult models only; full final sensitivity is still running',file.path(out,'run_status.txt'))
f4_figures(r,x,cfg,out)
