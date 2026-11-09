# Final metadata-mask annotation and independent integrity gate; never re-imputes.
options(fatty_liver.autorun=FALSE)
root<-Sys.getenv("FATTYLIVER_ROOT",unset="E:/YFS/workspace/FattyLiver_generation_analysis")
source(file.path(root,"scripts/12_impute_adult_covariates.R"),encoding="UTF-8")
version<-readLines(file.path(root,"data/imputed/CURRENT"),n=1)
folder<-file.path(root,"data/imputed",version);path<-file.path(folder,"adult_covariates_mi.rds")
x<-readRDS(path);stopifnot(isTRUE(x$diagnostics_passed),x$m==50L)
expected<-is.na(x$original_metadata[,fl_mi_targets]);expected[-x$adult_rows,]<-FALSE
if(!identical(x$imputation_flags,expected)){
 # Retain the initial generated bundle outside the public gallery. Only flags
 # are clarified; all 50 completed datasets and the MICE object are unchanged.
 module<-file.path(root,"outputs/12_Adult_covariate_MI/current",version)
 backup<-file.path(module,"private/bundle_before_mask_annotation.rds")
 if(file.exists(backup)){
  if(unname(tools::md5sum(backup))!=unname(tools::md5sum(path)))stop("Existing annotation backup differs from current initial bundle")
 }else if(!file.copy(path,backup))stop("Cannot preserve initial generated bundle")
 x$original_missingness<-x$imputation_flags
 x$imputation_flags<-expected
 x$flag_definition<-"TRUE only for originally missing adult target cells actually completed; minor missingness remains in original_missingness"
 tmp<-paste0(path,".validated");saveRDS(x,tmp,version=3)
 y<-readRDS(tmp);old<-readRDS(backup)
 # Formula environments gain distinct identities after separate unserialization.
 # Compare every numerical/state field directly and formulas without environment identity.
 state<-setdiff(names(y$mids),c("formulas","call"))
 stopifnot(identical(y$completed_metadata,old$completed_metadata),identical(y$mids[state],old$mids[state]),
  isTRUE(all.equal(y$mids$formulas,old$mids$formulas,check.environment=FALSE)),identical(y$original_metadata,old$original_metadata))
 if(!file.copy(tmp,path,overwrite=TRUE))stop("Cannot publish annotated bundle")
 unlink(tmp)
 manifest<-read.csv(file.path(folder,"dataset_manifest.csv"),fileEncoding="UTF-8-BOM")
 manifest$md5<-unname(tools::md5sum(path));fl_write_csv(manifest,file.path(folder,"dataset_manifest.csv"))
 fl_write_csv(data.frame(action="Clarified original-missing versus actually-imputed mask; completed values and mids unchanged",before_md5=unname(tools::md5sum(backup)),after_md5=unname(tools::md5sum(path))),file.path(module,"audit/mask_annotation.csv"))
}
checks<-fl_mi_checks(x$original_metadata,x$completed_metadata,x$adult_rows)
stopifnot(!any(x$imputation_flags[x$original_metadata$age<18,]),identical(unname(tools::md5sum(x$source_files$path)),x$source_files$md5))
dg<-fl_mi_diagnostics(x$mids);last<-dg[dg$.it==max(dg$.it),]
stopifnot(all(is.finite(last$psrf)),max(last$psrf)<=1.1)
fl_write_csv(checks,file.path(root,"outputs/12_Adult_covariate_MI/current",version,"audit/final_integrity_checks.csv"))
module<-file.path(root,"outputs/12_Adult_covariate_MI/current",version)
dist<-list();groups<-list()
for(g in c("G0","G1","G2")){
 j<-x$adult_rows[x$original_metadata$generation[x$adult_rows]==g]
 z<-x$original_metadata[j,];missing_any<-!complete.cases(z[,fl_mi_targets])
 for(k in c(FALSE,TRUE)){
  q<-z[missing_any==k,];groups[[length(groups)+1]]<-data.frame(generation=g,group=if(k)"Any target missing"else"All targets complete",n=nrow(q),cases=sum(q$fatty_liver),mean_age=mean(q$age),female_n=sum(q$sex=="Female"))
 }
 for(v in fl_mi_targets){
  observed<-as.numeric(as.character(z[[v]]));miss<-is.na(observed)
  imputed<-unlist(lapply(x$completed_metadata,function(q)as.numeric(as.character(q[[v]][j][miss]))))
  for(s in c("observed","imputed")){
   val<-if(s=="observed")observed[!miss]else imputed
   qs<-if(length(val))quantile(val,c(0,.25,.5,.75,1))else rep(NA_real_,5)
   dist[[length(dist)+1]]<-data.frame(generation=g,variable=v,source=s,n_people=if(s=="observed")sum(!miss)else sum(miss),
    n_values=length(val),min=qs[1],q1=qs[2],median=qs[3],q3=qs[4],max=qs[5])
  }
 }
}
fl_write_csv(do.call(rbind,dist),file.path(module,"tables/distribution_by_generation.csv"))
fl_write_csv(do.call(rbind,groups),file.path(module,"tables/complete_vs_missing_adult_characteristics.csv"))
code<-file.path(root,c("config.R","config/imputation_design_20261020.md","R/adult_covariate_mi.R","scripts/12_impute_adult_covariates.R","scripts/13_clinical_multiple_imputation.R","scripts/validate_adult_mi_delivery.R"))
fl_write_csv(data.frame(file=code,md5=unname(tools::md5sum(code))),file.path(module,"audit/final_code_fingerprints.csv"))
file.copy(file.path(root,"README_multiple_imputation_zh.md"),file.path(folder,"README_zh.md"),overwrite=TRUE)
file.copy(file.path(root,"README_multiple_imputation_en.md"),file.path(folder,"README_en.md"),overwrite=TRUE)
cat("PASS: final bundle integrity; ",nrow(checks)," checks; maximum Rhat=",max(last$psrf),"; original RDS unchanged\n",sep="")
