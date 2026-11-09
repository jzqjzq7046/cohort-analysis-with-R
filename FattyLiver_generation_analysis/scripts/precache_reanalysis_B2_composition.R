# Independent numerical fit only; the pipeline performs all inference/publication.
root<-'E:/YFS/workspace/FattyLiver_generation_analysis'
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),.libPaths()))
for(f in c('config.R','R/common.R','R/confounders.R','R/reanalysis_numerical_helpers.R','R/reanalysis_20261024.R'))source(file.path(root,f),encoding='UTF-8')
cfg<-fl_config(root);run<-readLines(file.path(root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),n=1)
out<-file.path(root,'outputs/24_Two_generation_reanalysis/_running',run)
run_fit<-function(){
 x<-readRDS(file.path(out,'private/stage_prepare.rds'));nm<-'B2';sc<-v3_micro_samples(x)[[nm]]
 d<-v3_design(x,sc$mask,sc$fields);m<-d$metadata;md<-d$design;md$fatty_liver<-m$fatty_liver
 A<-x$relative[,match(m$sample_id,colnames(x$relative)),drop=FALSE]
 Z<-model.matrix(reformulate(c('family_id',d$nuisance)),md);q<-qr(Z,tol=1e-9);keep<-q$pivot[seq_len(q$rank)];Z<-Z[,sort(keep),drop=FALSE]
 kept_nuis<-d$nuisance[d$nuisance%in%colnames(Z)];condition<-c('family_id',kept_nuis);Z<-model.matrix(reformulate(condition),md)
 stopifnot(qr(Z)$rank==ncol(Z));signature<-digest::digest(list(protocol=v3_protocol,relative=A,design=Z,outcome=m$fatty_liver,vegan='2.7.6'),algo='sha256')
 path<-file.path(out,'private/composition_fit_checkpoint_B2.rds')
 if(file.exists(path)){stopifnot(identical(readRDS(path)$signature,signature));return(invisible(NULL))}
 v3_log('Precomputing independent B2 fit: ',nrow(m),' people')
 D<-vegan::vegdist(t(A),method='bray');form<-as.formula(paste('D ~ fatty_liver + Condition(',paste(condition,collapse=' + '),')'))
 fit<-vegan::dbrda(form,data=md,add='lingoes');stopifnot(fit$CCA$qrank==1L);setup<-fl_dv_fast_setup(fit,m,Z)
 tmp<-paste0(path,'.pending');saveRDS(list(signature=signature,fit=fit,setup=setup),tmp,compress='gzip')
 if(file.exists(path)){stopifnot(identical(readRDS(path)$signature,signature));unlink(tmp)}else stopifnot(file.rename(tmp,path))
 v3_log('Independent B2 fit checkpoint ready')
}
run_fit()
