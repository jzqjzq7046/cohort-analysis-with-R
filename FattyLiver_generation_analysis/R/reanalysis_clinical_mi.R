# Updated continuation of the existing 12/13 MI analyses, restricted to C2 targets.
v3_clinical_mi <- function(x,cfg,out) {
 for(p in c('mice','miceadds','posterior'))if(!requireNamespace(p,quietly=TRUE))stop('Missing MI dependency: ',p)
 assign('mice.impute.v3_2l_pmm',function(y,ry,x,type,...){miceadds::mice.impute.2l.pmm(y,ry,x,type,donors=5,draw.fixed=TRUE,match_sampled_pars=TRUE,glmer.warnings=TRUE,optimizer='bobyqa')},envir=.GlobalEnv)
 cfg$cg_alpha<-.05;cfg$cg_max_iterations<-100L;cfg$cg_working_correlation<-'independence';pooled<-flows<-diagnostics<-list();packages<-list()
 for(g in v3_pairs){
  key<-paste(g,collapse='-');d<-fl_cg_make_comparison(v3_scope(x$clinical,g),g)$data;d<-d[v3_cc(d,c('age','sex'))&d$age>=18,,drop=FALSE];d<-d[order(d$family_id,d$participant_id),,drop=FALSE]
  dat<-data.frame(cluster_id=as.integer(factor(d$family_id)),age=d$age,female=as.integer(d$sex=='Female'),fatty_liver=d$fatty_liver,older_any_case=d$older_any_case,n_older_observed=d$n_older_observed,alcohol_log=log1p(d$alkodpd18),BMI18=d$BMI18)
  method<-setNames(rep('',ncol(dat)),names(dat));targets<-c('alcohol_log','BMI18');missing<-targets[colSums(is.na(dat[,targets,drop=FALSE]))>0];repeated<-anyDuplicated(dat$cluster_id)>0
  method[missing]<-if(repeated)'v3_2l_pmm'else'pmm'
  pred<-matrix(0,ncol(dat),ncol(dat),dimnames=list(names(dat),names(dat)))
  for(v in missing){pred[v,]<-1;pred[v,v]<-0;pred[v,'cluster_id']<-if(repeated)-2 else 0;constant<-names(dat)[vapply(dat,function(z)length(unique(na.omit(z)))<2,logical(1))];pred[v,constant]<-0}
  cache<-file.path(out,'private',paste0('clinical_MI_',key,'.rds'));signature<-digest::digest(list(data=dat,method=method,predictor=pred,M=50,iterations=10),algo='sha256')
  if(file.exists(cache)){old<-readRDS(cache);stopifnot(identical(old$signature,signature));imp<-old$imputation;v3_log('Resume adult MI: ',key)}else{
   v3_log('Adult clinical MI ',key,': N=',nrow(d),'; targets ',paste(missing,collapse=','),'; method ',if(repeated)'two-level PMM'else'PMM (one younger per family)')
   if(length(missing)){
    imp<-mice::mice(dat,m=50,maxit=10,method=method,predictorMatrix=pred,seed=20261024L+match(key,c('G0-G1','G1-G2','G0-G2')),printFlag=FALSE,donors=5,remove.constant=TRUE,remove.collinear=TRUE)
   }else imp<-NULL
   saveRDS(list(signature=signature,imputation=imp),cache,compress='gzip')
  }
  completed<-vector('list',50);fits<-list();trace<-list()
  if(!is.null(imp)){
   for(v in missing)for(parameter in c('mean','sd')){
    a<-if(parameter=='mean')imp$chainMean[v,,,drop=TRUE]else sqrt(imp$chainVar[v,,,drop=TRUE]);rh<-tryCatch(posterior::rhat(a),error=function(e)NA_real_)
    diagnostics[[paste(key,v,parameter)]]<-data.frame(comparison=key,variable=v,parameter=parameter,iterations=imp$iteration,imputations=imp$m,Rhat=rh,status=if(is.finite(rh)&&rh>1.1)'REVIEW_REQUIRED'else'DIAGNOSTIC_RECORDED')
   }
   if(!is.null(imp$loggedEvents))fl_write_csv(imp$loggedEvents,file.path(out,'audit',paste0('MI_logged_events_',key,'.csv')))
  }
  for(i in 1:50){
   filled<-if(is.null(imp))dat else mice::complete(imp,i);z<-d
   for(v in targets){j<-is.na(dat[[v]]);if(v=='alcohol_log')z$alkodpd18[j]<-expm1(filled[[v]][j])else z$BMI18[j]<-filled[[v]][j]}
   j<-is.na(d$alkodpd18);z$alcohol[j]<-z$alkodpd18[j]
   stopifnot(all(is.finite(z$alkodpd18)),all(z$alkodpd18>=0),all(is.finite(z$BMI18)),all(z$BMI18>0))
   for(v in setdiff(names(d),c('alkodpd18','alcohol','BMI18')))stopifnot(identical(z[[v]],d[[v]]))
   for(v in c('alkodpd18','BMI18')){j<-!is.na(d[[v]]);stopifnot(identical(z[[v]][j],d[[v]][j]))}
   completed[[i]]<-z
   for(model in c('Adult_C1_MI','C2_alcohol_MI','C2_BMI_MI')){
    cov<-c('age','sex','n_older_observed',if(model!='Adult_C1_MI')'alkodpd18',if(model=='C2_BMI_MI')'BMI18')
    spec<-data.frame(model=model,sample_set='adult_MI_same_sample',covariates=paste(cov,collapse='+'),primary=FALSE,references='Long2019;White2011',study_adaptation='Updated observed two-generation scope; age/sex/exposure/outcome fixed; alcohol/BMI MI only')
    fit<-fl_cg_fit(z,spec,key,cfg)$result;fit$imputation<-i;fits[[paste(model,i)]]<-fit
   }
  }
  f<-ov_bind(fits);fl_write_csv(f,file.path(out,'audit',paste0('MI_per_imputation_models_',key,'.csv')))
  for(model in unique(f$model)){
   a<-f[f$model==model,];ok<-is.finite(a$log_OR)&is.finite(a$robust_SE)&grepl('^ESTIMATED',a$status);row<-data.frame(comparison=key,model=model,N=nrow(d),cases=sum(d$fatty_liver),families=length(unique(d$family_id)),imputations=50,valid_imputations=sum(ok),OR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P_raw=NA_real_,FMI=NA_real_,status='NE_INCOMPLETE_VALID_IMPUTATIONS')
   if(all(ok)){
    pp<-mice::pool.scalar(a$log_OR,a$robust_SE^2,n=Inf,k=1);se<-sqrt(pp$t);critical<-qt(.975,pp$df);row$OR<-exp(pp$qbar);row$CI_low<-exp(pp$qbar-critical*se);row$CI_high<-exp(pp$qbar+critical*se);row$P_raw<-2*pt(-abs(pp$qbar/se),pp$df);row$FMI<-pp$fmi;row$status<-'ESTIMATED_MI_POOLED'
   }
   row$warning_imputations<-sum(nzchar(a$warnings));row$warnings<-paste(unique(a$warnings[nzchar(a$warnings)]),collapse='; ')
   if(row$warning_imputations>0&&row$status=='ESTIMATED_MI_POOLED')row$status<-'ESTIMATED_MI_POOLED_WITH_WARNINGS'
   pooled[[paste(key,model)]]<-row
  }
  flows[[key]]<-data.frame(comparison=key,adult_age_sex_complete_N=nrow(d),cases=sum(d$fatty_liver),complete_case_N=sum(complete.cases(dat[,targets])),alcohol_missing=sum(is.na(dat$alcohol_log)),BMI_missing=sum(is.na(dat$BMI18)),imputation_method=if(repeated)'Two-level PMM'else'PMM: singleton families')
  packages[[key]]<-list(original=d,completed=completed,imputation=imp,methods=method,predictor_matrix=pred)
 }
 result<-list(models=ov_bind(pooled),flow=ov_bind(flows),diagnostics=ov_bind(diagnostics));saveRDS(list(protocol=v3_protocol,M=50,comparisons=packages),file.path(out,'private','adult_clinical_imputed_two_generation.rds'),compress='gzip')
 v3_table(result$models,'Supplementary_Clinical_C2_multiple_imputation_EN',out);v3_table(result$flow,'Adult_clinical_MI_sample_flow_EN',out);v3_table(result$diagnostics,'Adult_clinical_MI_diagnostics_EN',out);result
}
v3_mi_visualize <- function(clinical,mi,cfg,out){
 cc<-clinical$models[clinical$models$model=='C2_BMI',c('comparison','n','cases','OR','CI_low','CI_high','p_value')];names(cc)<-c('comparison','N','cases','OR','CI_low','CI_high','P_raw');cc$analysis<-'Complete cases'
 pooled<-mi$models[mi$models$model=='C2_BMI_MI',c('comparison','N','cases','OR','CI_low','CI_high','P_raw')];pooled$analysis<-'Multiple imputation';z<-rbind(cc,pooled);z$comparison<-factor(z$comparison,levels=rev(c('G0-G1','G1-G2','G0-G2')));z$analysis<-factor(z$analysis,levels=c('Complete cases','Multiple imputation'))
 pd<-ggplot2::position_dodge(.45);p<-ggplot2::ggplot(z,ggplot2::aes(OR,comparison,colour=analysis,shape=analysis))+ggplot2::geom_vline(xintercept=1,linetype=2,colour='#888888',linewidth=.3)+ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low,xmax=CI_high),orientation='y',width=0,position=pd,colour='black',linewidth=.3)+ggplot2::geom_point(position=pd,size=1.8)+ggplot2::scale_x_log10()+ggplot2::scale_colour_manual(values=setNames(unname(cfg$binary_colors),levels(z$analysis)))+ggplot2::labs(x='Odds ratio (log scale)',y=NULL,colour=NULL,shape=NULL)+v3_style()
 ov_export(p,'Supplement3_Clinical_C2_multiple_imputation',out,cfg,183,105);v3_table(z,'Supplement3_MI_figure_source_EN',out)
 data.frame(figure='Supplement3_Clinical_C2_multiple_imputation',width_mm=183,height_mm=105,formats='PDF;PNG;SVG;TIFF',dpi=600)
}
