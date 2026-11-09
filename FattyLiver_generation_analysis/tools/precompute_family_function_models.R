# Independent function-model precomputation. Only complete verified official
# model checkpoints are atomically offered to the sequential pipeline; no figures
# or adopted outputs are published here and no species model is duplicated.
root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(root,'R/family_four_stage.R'),encoding='UTF-8')
cfg<-f4_setup(root);x<-abc_source(root)$x
module<-file.path(root,'outputs/48_Family_functional_explanation')
out<-file.path(module,'_precomputed',f4_protocol)
for(d in c('private','tables','audit'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
imp<-ade_import_paths(cfg,x,out);x<-imp$x
x$features$official_id<-x$features$feature_id;x$features$taxonomy<-''
x$relative<-x$relative[x$candidates,,drop=FALSE]
specs<-f4_models(x);cache<-file.path(out,'cache');target<-file.path(module,'_checkpoints',f4_protocol,'private')
dir.create(cache,recursive=TRUE,showWarnings=FALSE);dir.create(target,recursive=TRUE,showWarnings=FALSE)
for(k in names(specs)) {
  a<-specs[[k]];cat(format(Sys.time()),'PRECOMPUTE function',k,'\n');flush.console()
  fit<-abc_B_fit(x,a$m,a$fields,k,out,cache)
  cp<-file.path(cache,paste0(k,'.rds'));tmp<-file.path(target,paste0(k,'.precomputed.rds'))
  destination<-file.path(target,paste0(k,'.rds'))
  stopifnot(file.exists(cp));v<-readRDS(cp);stopifnot(!is.null(v$signature),!is.null(v$value$result))
  if(!file.exists(destination)) {
    stopifnot(file.copy(cp,tmp,overwrite=TRUE))
    if(!file.exists(destination))stopifnot(file.rename(tmp,destination))
  }
  cat(format(Sys.time()),'COMPLETE cached function',k,'\n');flush.console()
}
writeLines('COMPLETED_OFFICIAL_MODEL_PRECOMPUTATION',file.path(out,'precompute_status.txt'))
