# Adult-only multilevel MICE and family-GEE pooling. No original data are overwritten.
fl_mi_targets<-c("urbanicity","smoking","alcohol","diet_score","physical_activity","energy_log","BMI18")
fl_mi_root_data<-function(root){
 v<-readLines(file.path(root,"data/covariates/CURRENT"),n=1,warn=FALSE)
 path<-file.path(root,"data/covariates",v,"covariates.rds")
 pkg<-readRDS(path); m<-pkg$metadata
 stopifnot(nrow(m)==1275L,length(unique(m$family_id))==319L,all(table(m$family_id,m$generation)>0),!anyNA(m$age),!anyNA(m$fatty_liver))
 pv<-readLines(file.path(root,"data/processed/CURRENT"),n=1,warn=FALSE)
 protected<-c(path,file.path(root,"data/processed",pv,"microbiome_processed.rds"),
   file.path(root,"data/covariates/CURRENT"),file.path(root,"data/processed/CURRENT"))
 list(metadata=m,source=path,protected=data.frame(path=protected,md5=unname(tools::md5sum(protected))))
}
fl_mi_verify_source<-function(d){stopifnot(identical(unname(tools::md5sum(d$protected$path)),d$protected$md5));invisible(TRUE)}
fl_mi_prepare<-function(m){
 adult<-which(m$age>=18);adult<-adult[order(m$family_id[adult],m$participant_id[adult])];a<-m[adult,]
 stopifnot(length(adult)==1011L,all(na.omit(as.numeric(as.character(a$urbanicity)))%in%c(1,2)))
 # All fully observed auxiliaries are defined before imputation, from protected data.
 aux<-data.frame(age=a$age,sex=as.numeric(a$sex=="Male"),generation=factor(a$generation),case=a$fatty_liver,
  survey_center=factor(a$survey_center),shannon=a$shannon)
 for(g in c("G0","G1","G2")){
  z<-m[m$generation==g,];b<-split(seq_len(nrow(z)),z$family_id)
  sums<-do.call(rbind,lapply(names(b),function(f){q<-z[b[[f]],];data.frame(family_id=f,n=nrow(q),anycase=as.integer(any(q$fatty_liver==1)),meanage=mean(q$age),female=mean(q$sex=="Female"))}))
  k<-match(a$family_id,sums$family_id)
  for(v in c("n","anycase","meanage","female"))aux[[paste0(g,"_",v)]]<-sums[[v]][k]
 }
 stopifnot(!anyNA(aux))
 X<-model.matrix(~generation*(age+sex+case)+survey_center+shannon+
   G0_n+G0_anycase+G0_meanage+G0_female+G1_anycase+G1_meanage+G1_female+
   G2_n+G2_anycase+G2_meanage+G2_female+generation:G0_anycase+generation:G1_anycase,data=aux)
 dec<-qr(X,tol=1e-8);keep<-sort(dec$pivot[seq_len(dec$rank)]);dropped<-setdiff(colnames(X),colnames(X)[keep]);X<-X[,keep,drop=FALSE]
 X<-X[,colnames(X)!="(Intercept)",drop=FALSE];colnames(X)<-make.names(paste0("aux_",colnames(X)),unique=TRUE)
 X<-as.data.frame(X)
 for(v in names(X)){s<-sd(X[[v]]);if(s>0)X[[v]]<-(X[[v]]-mean(X[[v]]))/s}
 dat<-data.frame(cluster_id=match(a$family_id,unique(a$family_id)),urbanicity=as.numeric(as.character(a$urbanicity))-1,
  smoking=a$smoking,alcohol_log=log1p(a$alcohol),diet_score=a$diet_score,physical_activity=a$physical_activity,
  energy_log=a$energy_log,BMI18=a$BMI18,X,check.names=FALSE)
 targets<-c("urbanicity","smoking","alcohol_log","diet_score","physical_activity","energy_log","BMI18")
 # The miceadds manual explicitly demonstrates 2l.pmm for binary outcomes too.
 # Empirical donors preserve legal binary categories without perfect prediction
 # in a separate sparse logistic imputation model.
 method<-setNames(rep("",ncol(dat)),names(dat));method[targets]<-"fl2l_pmm"
 method[colSums(is.na(dat))==0]<-""
 pred<-matrix(0,ncol(dat),ncol(dat),dimnames=list(names(dat),names(dat)))
 for(v in targets[method[targets]!=""]){pred[v,]<-1;pred[v,v]<-0;pred[v,"cluster_id"]<--2}
 list(data=dat,method=method,predictor=pred,adult=adult,dropped=dropped)
}
fl_mi_register_methods<-function(){
 # Published miceadds methods; all random draws use mice's recorded seed streams.
 assign("mice.impute.fl2l_pmm",function(y,ry,x,type,...){
  miceadds::mice.impute.2l.pmm(y,ry,x,type,donors=5,draw.fixed=TRUE,match_sampled_pars=TRUE,
    glmer.warnings=TRUE,optimizer="bobyqa")
 },envir=.GlobalEnv)
}
fl_mi_decode<-function(completed,m,adult){
 d<-m
 for(v in fl_mi_targets){miss<-adult[is.na(m[[v]][adult])];if(!length(miss))next
  val<-if(v=="alcohol")expm1(completed$alcohol_log) else if(v=="urbanicity")completed$urbanicity+1 else completed[[v]]
  j<-match(miss,adult)
  if(v=="urbanicity")d[[v]][miss]<-as.character(val[j]) else d[[v]][miss]<-val[j]
 }
 # Keep harmonized aliases coherent, changing only originally missing adult cells.
 for(pair in list(c("alcohol","alkodpd18"),c("smoking","smoke18"),c("diet_score","dietscore_g18"))){
  j<-adult[is.na(m[[pair[2]]][adult])];d[[pair[2]]][j]<-d[[pair[1]]][j]
 }
 d
}
fl_mi_checks<-function(original,completed,adult){
 rows<-list();add<-function(name,ok){if(!isTRUE(ok))stop("Imputation integrity failed: ",name);rows[[length(rows)+1]]<<-data.frame(check=name,passed=TRUE)}
 editable<-c(fl_mi_targets,"alkodpd18","smoke18","dietscore_g18");minor<-which(original$age<18)
 for(i in seq_along(completed)){z<-completed[[i]]
  add(paste(i,"all original non-target fields preserved"),identical(z[,setdiff(names(z),editable)],original[,setdiff(names(z),editable)]))
  add(paste(i,"all minors unchanged"),identical(z[minor,],original[minor,]))
  for(v in editable){j<-which(!is.na(original[[v]]));add(paste(i,v,"observed values unchanged"),identical(z[[v]][j],original[[v]][j]))}
  add(paste(i,"adult targets complete"),all(complete.cases(z[adult,fl_mi_targets])))
  add(paste(i,"binary categories valid"),all(z$smoking[adult]%in%0:1)&&all(as.character(z$urbanicity[adult])%in%c("1","2")))
  for(v in c("diet_score","physical_activity","energy_log","BMI18"))add(paste(i,v,"observed donor support"),all(z[[v]][adult]%in%na.omit(original[[v]][adult])))
  add(paste(i,"alcohol valid"),all(is.finite(z$alcohol[adult])&z$alcohol[adult]>=0))
 }
 do.call(rbind,rows)
}
fl_mi_diagnostics<-function(imp){
 # Same published rank-normalized split/folded Rhat, without requiring Stan.
 vars<-names(imp$method)[nzchar(imp$method)];rows<-list()
 for(v in vars)for(par in c("mean","sd")){
  a<-if(par=="mean")imp$chainMean[v,,,drop=TRUE] else sqrt(imp$chainVar[v,,,drop=TRUE])
  for(k in seq.int(min(4L,imp$iteration),imp$iteration)){
   x<-a[seq_len(k),,drop=FALSE]
   ac<-apply(x,2,function(q)if(length(q)>2&&sd(q)>0&&sd(head(q,-1))>0&&sd(tail(q,-1))>0)cor(head(q,-1),tail(q,-1))else NA_real_)
   rows[[length(rows)+1]]<-data.frame(.it=k,vrb=v,parameter=par,psrf=posterior::rhat(x),
    ac=if(any(is.finite(ac)))mean(ac,na.rm=TRUE)else NA_real_,constant_chains=sum(apply(x,2,sd)==0),
    diagnostic="posterior rank-normalized split/folded Rhat")
  }
 }
 do.call(rbind,rows)
}
fl_mi_impute_main<-function(args,root){
 need<-c("mice","miceadds","lme4","jsonlite","posterior");stopifnot(all(vapply(need,requireNamespace,logical(1),quietly=TRUE)))
 opt<-fl_entry_options(args,fl_config(root),"adult_covariate_mi");if(opt$mode=="--help"){cat("12_impute_adult_covariates.R [--run | --check] [--run-id=NAME]\n");return(invisible())}
 if(opt$mode=="--check"){fl_mi_root_data(root);cat("PASS: sources and runtime\n");return(invisible())}
 ctx<-fl_output_start_run(opt$cfg);out<-ctx$run_dir;for(s in c("tables","reports","audit","private"))dir.create(file.path(out,s))
 writeLines("RUNNING",file.path(out,"run_status.txt"))
 tryCatch({
  d<-fl_mi_root_data(root);m<-d$metadata;prep<-fl_mi_prepare(m);fl_mi_register_methods()
  fl_write_csv(d$protected,file.path(out,"audit/original_file_fingerprints.csv"))
  fl_write_csv(data.frame(column=prep$dropped),file.path(out,"audit/QR_aliases.csv"))
  fl_write_csv(data.frame(variable=names(prep$method),method=prep$method),file.path(out,"audit/methods.csv"))
  fl_write_csv(data.frame(variable=rownames(prep$predictor),prep$predictor,check.names=FALSE),file.path(out,"audit/predictor_matrix.csv"))
  missing<-do.call(rbind,lapply(c("G0","G1","G2"),function(g){j<-prep$adult[m$generation[prep$adult]==g]
   do.call(rbind,lapply(fl_mi_targets,function(v)data.frame(generation=g,variable=v,n=length(j),missing=sum(is.na(m[[v]][j])),cases_missing=sum(m$fatty_liver[j][is.na(m[[v]][j])]))))}))
  fl_write_csv(missing,file.path(out,"tables/adult_missingness.csv"))
  warnings<-character();log_condition<-function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}
  cat("Starting multilevel MICE: 50 chains, adult N=",length(prep$adult),"; original minors and outcomes protected\n",sep="");flush.console()
  imp<-withCallingHandlers(mice::mice(prep$data,m=50,maxit=10,method=prep$method,predictorMatrix=prep$predictor,
    seed=20261020L,printFlag=FALSE,remove.constant=TRUE,remove.collinear=FALSE),warning=log_condition,message=function(x)invokeRestart("muffleMessage"))
  saveRDS(imp,file.path(out,"private/imputation_checkpoint.rds"));cat("MICE iteration 10 completed\n");flush.console()
  repeat{
   if(imp$iteration>=40){dg<-fl_mi_diagnostics(imp);last<-dg[dg$.it==max(dg$.it),];finite<-last$psrf[is.finite(last$psrf)]
    if(length(finite)&&max(finite)<=1.1||imp$iteration>=100)break}
   imp<-withCallingHandlers(mice::mice.mids(imp,maxit=10,printFlag=FALSE),warning=log_condition,message=function(x)invokeRestart("muffleMessage"))
   saveRDS(imp,file.path(out,"private/imputation_checkpoint.rds"));cat("MICE iteration ",imp$iteration," completed\n",sep="");flush.console()
  }
  dg<-fl_mi_diagnostics(imp);fl_write_csv(dg,file.path(out,"audit/chain_convergence.csv"))
  last<-dg[dg$.it==max(dg$.it),];psrf<-last$psrf[is.finite(last$psrf)];diag_ok<-length(psrf)>0&&max(psrf)<=1.1
  fl_write_csv(data.frame(warning=unique(warnings),occurrences=as.integer(table(factor(warnings,levels=unique(warnings))))),file.path(out,"audit/imputation_warnings.csv"))
  if(!is.null(imp$loggedEvents))fl_write_csv(imp$loggedEvents,file.path(out,"audit/mice_logged_events.csv"))
  complete<-lapply(seq_len(imp$m),function(i)fl_mi_decode(mice::complete(imp,i),m,prep$adult))
  checks<-fl_mi_checks(m,complete,prep$adult);fl_write_csv(checks,file.path(out,"audit/integrity_checks.csv"))
  dist<-do.call(rbind,lapply(fl_mi_targets,function(v){j<-prep$adult;obs<-as.numeric(as.character(m[[v]][j]));miss<-is.na(obs)
    vals<-unlist(lapply(complete,function(z)as.numeric(as.character(z[[v]][j][miss]))));obs<-obs[!is.na(obs)]
    do.call(rbind,lapply(c("observed","imputed"),function(k){x<-if(k=="observed")obs else vals;qs<-if(length(x))quantile(x,c(0,.25,.5,.75,1))else rep(NA_real_,5)
     data.frame(variable=v,source=k,n_values=length(x),min=qs[1],q1=qs[2],median=qs[3],q3=qs[4],max=qs[5])}))}))
  fl_write_csv(dist,file.path(out,"tables/observed_imputed_distribution.csv"))
  # Patient-level missingness masks and completed datasets remain outside the public gallery.
  version<-ctx$cfg$run_id;folder<-file.path(root,"data/imputed",version);if(dir.exists(folder))stop("Imputed version already exists")
  dir.create(folder,recursive=TRUE);writeLines("RUNNING",file.path(folder,"STATUS"))
  source_missing<-is.na(m[,fl_mi_targets]);actual_imputed<-source_missing;actual_imputed[-prep$adult,]<-FALSE
  pkg<-list(schema_version=1L,version=version,method="multilevel MICE; random-intercept predictive mean matching including binary donors",m=imp$m,iterations=imp$iteration,
    seed=20261020L,adult_rows=prep$adult,original_metadata=m,completed_metadata=complete,imputation_flags=actual_imputed,original_missingness=source_missing,
    mids=imp,source_files=d$protected,diagnostics_passed=diag_ok,maximum_final_psrf=if(length(psrf))max(psrf)else NA_real_,
    scope="Adult covariates only; fixed 319 three-generation families; minors, outcomes and microbiome unchanged",
    missingness_assumption="MAR conditional on predictors; MNAR not excluded",design=file.path(root,"config/imputation_design_20261020.md"))
  path<-file.path(folder,"adult_covariates_mi.rds");saveRDS(pkg,path,version=3);fl_mi_verify_source(d)
  fl_write_csv(data.frame(file="adult_covariates_mi.rds",md5=unname(tools::md5sum(path)),m=imp$m,adult_n=length(prep$adult),diagnostics_passed=diag_ok),file.path(folder,"dataset_manifest.csv"))
  writeLines(c("# 成人协变量多重插补 / Adult covariate multiple imputation","",
    "原始RDS未修改。新RDS保存original_metadata、50套completed_metadata及mids；不能只取一套作最终推断。",
    "Original RDS files are unchanged. Analyse all 50 completed_metadata datasets and pool estimates; never treat a single completion as observed truth.",
    "仅插补成人协变量，未成年人、结局及菌群不变。使用scripts/13_clinical_multiple_imputation.R分析。",
    "Adult covariates only. Minors, outcomes and microbiome measurements are preserved. See the linked project design for assumptions and limitations."),file.path(folder,"README_zh_en.md"),useBytes=TRUE)
  writeLines(if(diag_ok)"COMPLETED"else"DIAGNOSTICS_FAILED",file.path(folder,"STATUS"))
  fl_write_csv(data.frame(data_path=path,adult_n=length(prep$adult),adult_cases=sum(m$fatty_liver[prep$adult]),m=imp$m,iterations=imp$iteration,maximum_final_psrf=pkg$maximum_final_psrf,diagnostics_passed=diag_ok),file.path(out,"tables/imputation_summary.csv"))
  file.copy(file.path(root,"config/imputation_design_20261020.md"),file.path(out,"reports/methods_and_references_zh.md"))
  capture.output(sessionInfo(),file=file.path(out,"audit/sessionInfo.txt"))
  writeLines(c("# 成人协变量插补完成",paste0("新RDS：",path),paste0("成人1011人；插补",imp$m,"套；迭代",imp$iteration,"轮；最终最大PSRF=",signif(pkg$maximum_final_psrf,4)),
   "原数据MD5已核验，264名未成年人及所有已观察值保持不变。", "仅处理成人临床扩展所需协变量，不插补菌群、不增加病例。"),file.path(out,"reports/README_zh.md"),useBytes=TRUE)
  if(!diag_ok)stop("Convergence diagnostic gate not passed; data saved for inspection but not adopted")
  writeLines(version,file.path(root,"data/imputed/CURRENT"));writeLines("COMPLETED",file.path(out,"run_status.txt"));final<-fl_output_complete_run(ctx)
  cat("Completed imputation: ",path,"\nAdopted output: ",final,"\n",sep="");invisible(pkg)
 },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(e)})
}
fl_mi_pool<-function(z,m){
 good<-grepl("^ESTIMATED",z$status)&is.finite(z$log_OR)&is.finite(z$robust_SE)&z$robust_SE>0
 r<-data.frame(comparison=z$comparison[1],model=z$model[1],n=z$n[1],cases=z$cases[1],families=z$families[1],
  m_total=m,m_estimable=sum(good),status="NOT_POOLABLE",OR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p_raw=NA_real_,
  log_OR=NA_real_,within_variance=NA_real_,between_variance=NA_real_,total_variance=NA_real_,df=NA_real_,FMI=NA_real_,
  MCSE_log_OR=NA_real_,MCSE_fraction_total_SE=NA_real_,reason=paste(unique(z$status[!good]),collapse=";"))
 if(nrow(z)!=m||!all(good))return(r)
 p<-mice::pool.scalar(z$log_OR,z$robust_SE^2,n=Inf,k=1)
 crit<-qt(.975,p$df);r$status<-"POOLED";r$log_OR<-p$qbar;r$within_variance<-p$ubar;r$between_variance<-p$b;r$total_variance<-p$t;r$df<-p$df;r$FMI<-p$fmi
 r$OR<-exp(p$qbar);r$CI_low<-exp(p$qbar-crit*sqrt(p$t));r$CI_high<-exp(p$qbar+crit*sqrt(p$t));r$p_raw<-2*pt(-abs(p$qbar/sqrt(p$t)),df=p$df)
 r$MCSE_log_OR<-sqrt(p$b/m);r$MCSE_fraction_total_SE<-r$MCSE_log_OR/sqrt(p$t);r$reason<-"";r
}
fl_mi_clinical_main<-function(args,root){
 opt<-fl_entry_options(args,fl_cg_config(root),"clinical_mi");cfg<-opt$cfg
 pointer<-file.path(root,"data/imputed/CURRENT");if(!file.exists(pointer))stop("Run independent 12 imputation entry first")
 folder<-file.path(root,"data/imputed",readLines(pointer,n=1));path<-file.path(folder,"adult_covariates_mi.rds")
 if(readLines(file.path(folder,"STATUS"),n=1)!="COMPLETED")stop("Imputed data not validated")
 man<-read.csv(file.path(folder,"dataset_manifest.csv"),fileEncoding="UTF-8-BOM");stopifnot(unname(tools::md5sum(path))==man$md5)
 if(opt$mode=="--check"){cat("PASS: completed imputation available\n");return(invisible())}
 ctx<-fl_output_start_run(cfg);out<-ctx$run_dir;for(s in c("tables","reports","audit","private"))dir.create(file.path(out,s));writeLines("RUNNING",file.path(out,"run_status.txt"))
 tryCatch({
  pkg<-readRDS(path);stopifnot(isTRUE(pkg$diagnostics_passed));protected_before<-tools::md5sum(pkg$source_files$path)
  specs<-fl_cg_review_specs();specs<-specs[specs$sample_set=="adult_review",];fits<-list();cells<-list();count<-0L
  for(i in seq_len(pkg$m)){
   mm<-pkg$completed_metadata[[i]]
   for(g in cfg$cg_pairs){cmp<-paste(g,collapse="-");z<-fl_cg_review_older(mm,g);z<-z[z$age>=18,,drop=FALSE]
    for(s in seq_len(nrow(specs))){sp<-specs[s,,drop=FALSE];v<-strsplit(sp$complete_case_fields,"+",fixed=TRUE)[[1]];stopifnot(all(complete.cases(z[,v])))
     r<-fl_cg_fit(z,sp,cmp,cfg);count<-count+1L;r$result$imputation<-i;r$cells$imputation<-i;fits[[count]]<-r$result;cells[[count]]<-r$cells}
   }
   if(i%%5==0){cat("Clinical MI: ",i,"/",pkg$m," imputations evaluated\n",sep="");flush.console()}
  }
  all<-do.call(rbind,fits);pooled<-do.call(rbind,lapply(split(all,paste(all$comparison,all$model)),fl_mi_pool,m=pkg$m))
  fl_write_csv(all,file.path(out,"audit/per_imputation_models.csv"));fl_write_csv(do.call(rbind,cells),file.path(out,"tables/exposure_case_cells.csv"))
  fl_write_csv(pooled,file.path(out,"tables/pooled_adult_models.csv"))
  # Never discard separated imputations and pool only the successful subset.
  oldroot<-file.path(root,"outputs/06_Cross_generation_clinical/current");oldruns<-list.dirs(oldroot,recursive=FALSE,full.names=TRUE);stopifnot(length(oldruns)==1)
  old<-read.csv(file.path(oldruns,"clinical/review_models.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);old<-old[old$sample_set=="adult_review",]
  cc<-data.frame(comparison=old$comparison,model=old$model,n=old$n,cases=old$cases,OR=old$OR,CI_low=old$CI_low,CI_high=old$CI_high,p_raw=old$p_value,status=old$status,analysis="Complete case")
  mi<-pooled[,c("comparison","model","n","cases","OR","CI_low","CI_high","p_raw","status")];mi$analysis<-"Multiple imputation"
  comparison<-rbind(cc,mi);fl_write_csv(comparison,file.path(out,"tables/complete_case_vs_MI.csv"))
  validations<-list();add<-function(name,ok){if(!isTRUE(ok))stop("MI analysis validation failed: ",name);validations[[length(validations)+1L]]<<-data.frame(check=name,passed=TRUE)}
  add("All 450 model evaluations saved",nrow(all)==9*pkg$m)
  for(cmp in unique(all$comparison)){z<-all[all$comparison==cmp,];expected<-if(cmp=="G0-G1")c(319,88)else c(224,12)
   add(paste(cmp,"restored adult sample and observed cases"),all(z$n==expected[1])&&all(z$cases==expected[2])&&all(z$minors_n==0))}
  for(j in seq_len(nrow(pooled))){r<-pooled[j,];z<-all[all$comparison==r$comparison&all$model==r$model,]
   if(r$status=="POOLED"){
    q<-mean(z$log_OR);u<-mean(z$robust_SE^2);b<-var(z$log_OR);t<-u+(1+1/pkg$m)*b
    add(paste(r$comparison,r$model,"independent Rubin reconstruction"),abs(q-r$log_OR)<1e-12&&abs(t-r$total_variance)<1e-12)
   }else add(paste(r$comparison,r$model,"no selective pooling"),is.na(r$OR)&&r$m_estimable<pkg$m)
  }
  add("Original RDS and pointers unchanged",identical(unname(protected_before),pkg$source_files$md5)&&identical(unname(tools::md5sum(pkg$source_files$path)),pkg$source_files$md5))
  fl_write_csv(do.call(rbind,validations),file.path(out,"audit/validation_checks.csv"))
  fl_write_csv(data.frame(path=c(path,file.path(oldruns,"clinical/review_models.csv")),md5=unname(tools::md5sum(c(path,file.path(oldruns,"clinical/review_models.csv"))))),file.path(out,"audit/input_fingerprints.csv"))
  notes<-c("# 成人临床多重插补重分析", "",paste0("方法：多层MICE，",pkg$m,"套；Rubin合并log OR及家系稳健方差。"),
   "仅成人协变量插补；原始RDS、结局、未成年人及菌群未改变。原06全龄主结果仍保留。",
   "G0-G1：319人、88例；G1-G2及G0-G2：224人、12例。", "",
   "| 比较 | 模型 | N/病例 | 可估套数 | 合并OR（95%CI） | 原始P | 状态 |", "|---|---|---:|---:|---|---:|---|")
  notes<-c(notes,vapply(seq_len(nrow(pooled)),function(j){r<-pooled[j,];val<-if(r$status=="POOLED")sprintf("%.3f (%.3f–%.3f)",r$OR,r$CI_low,r$CI_high)else"NE"
   paste0("| ",r$comparison," | ",r$model," | ",r$n,"/",r$cases," | ",r$m_estimable,"/",r$m_total," | ",val," | ",if(is.na(r$p_raw))"NE"else format(r$p_raw,digits=4)," | ",r$status," |")},character(1)),"",
   "NOT_POOLABLE表示至少一套发生分离等问题，不能只挑可估套数合并；不代表无关联。插补不能创造病例，模型稀疏问题可能仍存在。",
   "MI使用MAR假设。供体匹配用于插补连续、评分及二分类协变量；后续仍是原GEE，不换成忽略家系的普通Firth。",
   "完整资料与MI的样本范围不同，变化不能全部解释为协变量调整作用。所有P原始、双侧，95%CI逐项计算。")
  writeLines(enc2utf8(notes),file.path(out,"reports/results_summary_zh.md"),useBytes=TRUE)
  fl_ap_page(file.path(out,"reports/results.html"),"成人协变量插补后临床重分析",notes[3:5],list("MI合并结果"=pooled,"完整资料与MI并列"=comparison))
  file.copy(file.path(root,"config/imputation_design_20261020.md"),file.path(out,"reports/methods_and_references_zh.md"))
  capture.output(sessionInfo(),file=file.path(out,"audit/sessionInfo.txt"))
  writeLines("COMPLETED",file.path(out,"run_status.txt"));final<-fl_output_complete_run(ctx)
  cat("Completed clinical MI analysis: ",final,"\n",sep="");print(pooled[,c("comparison","model","n","cases","m_estimable","status","OR","CI_low","CI_high","p_raw")]);invisible(pooled)
 },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(e)})
}
