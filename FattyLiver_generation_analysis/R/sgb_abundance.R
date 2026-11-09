# Literature-based SGB abundance analysis; original inputs are read-only.
fl_ab_prepare <- function(cfg) {
  x<-fl_sgb_prepare(cfg)
  x$reference_mean<-colMeans(x$log10);x$reference_sd<-apply(x$log10,2,sd)
  stopifnot(all(is.finite(x$reference_sd)),all(x$reference_sd>0))
  x$standardized<-sweep(sweep(x$log10,2,x$reference_mean,"-"),2,x$reference_sd,"/")
  stopifnot(max(abs(colMeans(x$standardized)))<1e-8)
  x
}

fl_ab_scenarios <- function(x) {
  m<-x$metadata;adult<-m$age>=18&complete.cases(m[,c(fl_cv_life(),"BMI18"),drop=FALSE])
  list(primary=list(mask=rep(TRUE,nrow(m)),extra=character()),
    no_antibiotics=list(mask=m$abx_recent==0,extra=character()),
    adult_core=list(mask=adult,extra=character()),
    adult_lifestyle=list(mask=adult,extra=fl_cv_life()),
    adult_BMI=list(mask=adult,extra=c(fl_cv_life(),"BMI18")))
}

fl_ab_sample_audit <- function(x,scenarios,out) {
  m<-x$metadata;adult<-m$age>=18;included<-scenarios$adult_core$mask
  flow<-data.frame(stage=c("Core all ages","Core adults","Adult common complete cases",
    "Adults excluded for missing lifestyle covariates or BMI"),
    n=c(nrow(m),sum(adult),sum(included),sum(adult&!included)),
    cases=c(sum(m$fatty_liver),sum(m$fatty_liver[adult]),
      sum(m$fatty_liver[included]),sum(m$fatty_liver[adult&!included])))
  fields<-c(fl_cv_life(),"BMI18")
  missing<-data.frame(variable=fields,adult_n=sum(adult),
    missing_n=vapply(fields,function(v)sum(is.na(m[[v]][adult])),integer(1)))
  missing$missing_percent<-100*missing$missing_n/missing$adult_n
  strata<-do.call(rbind,lapply(c("G0","G1","G2"),function(g){
    a<-adult&m$generation==g;i<-included&m$generation==g
    data.frame(generation=g,eligible_adults=sum(a),eligible_cases=sum(m$fatty_liver[a]),
      included_adults=sum(i),included_cases=sum(m$fatty_liver[i]))
  }))
  fl_write_csv(flow,file.path(out,"tables/SGB_adult_complete_case_flow_EN.csv"))
  fl_write_csv(missing,file.path(out,"tables/SGB_adult_covariate_missingness_EN.csv"))
  fl_write_csv(strata,file.path(out,"tables/SGB_adult_inclusion_by_generation_EN.csv"))
  invisible(flow)
}

fl_ab_design <- function(x,scenario) {
  ii<-which(scenario$mask);m<-droplevels(x$metadata[ii,,drop=FALSE])
  z<-x$design[ii,x$mapping$model_column,drop=FALSE]
  original<-setNames(x$mapping$original_column,x$mapping$model_column)
  for(v in scenario$extra){z[[v]]<-m[[v]];original[v]<-v}
  stopifnot(!anyNA(z))
  constant<-names(z)[vapply(z,function(a)length(unique(a))<2L,logical(1))]
  z<-z[,setdiff(names(z),constant),drop=FALSE]
  for(v in names(z))z[[v]]<-(z[[v]]-mean(z[[v]]))/sd(z[[v]])
  Z<-model.matrix(~.,z);q<-qr(Z);keep<-sort(q$pivot[seq_len(q$rank)])
  aliases<-setdiff(colnames(Z),colnames(Z)[keep]);Z<-Z[,keep,drop=FALSE]
  Z<-Z[,colnames(Z)!="(Intercept)",drop=FALSE]
  names_new<-sprintf("nuis_%03d",seq_len(ncol(Z)))
  mapping<-data.frame(model_column=names_new,original_column=unname(original[colnames(Z)]))
  colnames(Z)<-names_new
  md<-data.frame(family_id=droplevels(m$family_id),fatty_liver=m$fatty_liver,Z,row.names=m$sample_id,check.names=FALSE)
  for(g in c("G0","G1","G2"))md[[paste0("case_",g)]]<-as.integer(m$generation==g)*m$fatty_liver
  for(terms in list("fatty_liver",paste0("case_",c("G0","G1","G2")))){
    A<-model.matrix(reformulate(c(terms,names_new)),md)
    if(qr(A)$rank<ncol(A))stop("Protected disease design is rank deficient")
  }
  list(index=ii,metadata=m,design=md,covariates=names_new,mapping=mapping,
    dropped=c(constant,aliases),n=nrow(m),cases=sum(m$fatty_liver),families=nlevels(md$family_id))
}

fl_ab_extract <- function(ans,terms,scenario,stage) {
  effects<-diagnostics<-hetero<-list()
  stopifnot(all(c("results","fits")%in%names(ans)),!is.null(names(ans$fits)))
  for(feature in names(ans$fits)){
    fit<-ans$fits[[feature]];state<-fl_sgb_status(fit)
    cs<-if(inherits(fit,"merMod"))tryCatch(coef(summary(fit)),error=function(e)NULL)else NULL
    warn<-if(inherits(fit,"merMod"))paste(fit@optinfo$conv$lme4$messages,collapse="; ")else""
    for(term in terms){
      ok<-startsWith(state,"ESTIMATED")&&!is.null(cs)&&term%in%rownames(cs)
      b<-se<-df<-p<-lo<-hi<-NA_real_
      if(ok){
        b<-cs[term,"Estimate"];se<-cs[term,"Std. Error"];df<-cs[term,"df"];p<-cs[term,"Pr(>|t|)"]
        ok<-all(is.finite(c(b,se,df,p)))&&se>0&&df>0
        if(ok){critical<-qt(.975,df);lo<-b-critical*se;hi<-b+critical*se}
        rr<-ans$results
        rr<-rr[rr$feature==feature&rr$metadata==term,,drop=FALSE]
        if(nrow(rr)==1L&&ok&&(!isTRUE(all.equal(unname(b),rr$coef[1],tolerance=1e-7))||
          !isTRUE(all.equal(unname(p),rr$pval[1],tolerance=1e-7))))stop("Official package/fit extraction mismatch: ",feature," / ",term)
      }
      effects[[length(effects)+1L]]<-data.frame(feature_id=feature,scenario=scenario,stage=stage,
        term=term,estimate=if(ok)b else NA_real_,SE=if(ok)se else NA_real_,df=df,
        CI_low=if(ok)lo else NA_real_,CI_high=if(ok)hi else NA_real_,p_raw=if(ok)p else NA_real_,
        status=if(ok)state else if(startsWith(state,"ESTIMATED"))"INVALID_COEFFICIENT"else state,
        warning=warn,stringsAsFactors=FALSE)
    }
    icc<-skew<-NA_real_
    if(inherits(fit,"merMod")){
      v<-as.data.frame(lme4::VarCorr(fit));vf<-sum(v$vcov[v$grp!="Residual"]);vr<-sigma(fit)^2
      icc<-vf/(vf+vr);r<-residuals(fit);skew<-mean((r-mean(r))^3)/sd(r)^3
    }
    diagnostics[[length(diagnostics)+1L]]<-data.frame(feature_id=feature,scenario=scenario,stage=stage,
      status=state,family_ICC=icc,residual_skewness=skew,warning=warn)
    if(stage=="generation"){
      statistic<-p_int<-df2<-NA_real_;hstate<-state
      if(startsWith(state,"ESTIMATED")){
        h<-tryCatch({
          b<-lme4::fixef(fit);L<-matrix(0,2,length(b),dimnames=list(NULL,names(b)))
          L[,terms]<-rbind(c(1,-1,0),c(0,1,-1))
          lmerTest::contestMD(fit,L,rhs=c(0,0),ddf="Satterthwaite")
        },error=function(e)NULL)
        if(is.null(h))hstate<-"NOT_ESTIMABLE_CONTRAST"else{
          statistic<-h[["F value"]][1];p_int<-h[["Pr(>F)"]][1];df2<-h[["DenDF"]][1]
          if(!all(is.finite(c(statistic,p_int,df2))))hstate<-"NOT_ESTIMABLE_CONTRAST"
        }
      }
      hetero[[length(hetero)+1L]]<-data.frame(feature_id=feature,statistic=statistic,df1=2,df2=df2,
        p_interaction_raw=p_int,status=hstate,test="Satterthwaite 2-df F")
    }
  }
  list(effects=do.call(rbind,effects),diagnostics=do.call(rbind,diagnostics),heterogeneity=do.call(rbind,hetero))
}

fl_ab_stage <- function(x,d,cfg,out,scenario,stage="overall",feature_index=seq_len(ncol(x$standardized))) {
  terms<-if(stage=="generation")paste0("case_",c("G0","G1","G2"))else"fatty_liver"
  Y<-x$standardized[d$index,feature_index,drop=FALSE]
  stopifnot(identical(rownames(Y),rownames(d$design)),all(is.finite(Y)))
  batches<-split(seq_len(ncol(Y)),ceiling(seq_len(ncol(Y))/20L));pieces<-ledger<-list()
  destination<-file.path(out,"private",scenario,stage);dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  fl_write_csv(d$mapping,file.path(out,"audit",paste0(scenario,"_nuisance_mapping.csv")))
  fl_write_csv(data.frame(dropped=d$dropped),file.path(out,"audit",paste0(scenario,"_dropped_nuisance.csv")))
  for(i in seq_along(batches)){
    ids<-batches[[i]];folder<-file.path(destination,sprintf("batch_%03d",i));dir.create(folder,showWarnings=FALSE)
    checkpoint<-file.path(folder,"statistics.rds")
    signature<-digest::digest(list(Y[,ids,drop=FALSE],d$design,terms,d$covariates,
      packageVersion("Maaslin2"),packageVersion("lmerTest"),body(fl_ab_extract)),algo="sha256")
    start<-Sys.time();cat(format(start)," ",scenario,"/",stage," batch ",i,"/",length(batches),"\n",sep="");flush.console()
    if(file.exists(checkpoint)){
      cache<-readRDS(checkpoint);if(!identical(cache$signature,signature))stop("Checkpoint signature mismatch")
      pieces[[i]]<-cache$statistics;rm(cache)
    }else{
      fit<-Maaslin2::Maaslin2(input_data=as.data.frame(Y[,ids,drop=FALSE]),input_metadata=d$design,
        output=folder,fixed_effects=c(terms,d$covariates),random_effects="family_id",
        normalization="NONE",transform="NONE",standardize=FALSE,
        min_prevalence=0,min_abundance=-Inf,min_variance=0,analysis_method="LM",
        cores=cfg$ab_cores,plot_heatmap=FALSE,plot_scatter=FALSE,save_models=TRUE)
      pieces[[i]]<-fl_ab_extract(fit,terms,scenario,stage)
      tmp<-paste0(checkpoint,".tmp");saveRDS(list(signature=signature,statistics=pieces[[i]]),tmp)
      if(!file.rename(tmp,checkpoint))stop("Could not commit batch checkpoint")
      rm(fit);gc(FALSE)
    }
    ledger[[i]]<-data.frame(scenario=scenario,stage=stage,batch=i,completed_features=max(ids),
      total_features=ncol(Y),elapsed_seconds=as.numeric(difftime(Sys.time(),start,units="secs")),status="COMPLETE")
    fl_write_csv(do.call(rbind,ledger),file.path(out,"audit",paste0(scenario,"_",stage,"_progress.csv")))
  }
  setNames(lapply(names(pieces[[1]]),function(v){z<-lapply(pieces,`[[`,v);do.call(rbind,z)}),names(pieces[[1]]))
}

fl_ab_linda <- function(x,d,cfg,out) {
  meta<-d$design;meta$fatty_liver<-factor(meta$fatty_liver,levels=c(0,1),labels=c("No","Yes"))
  formula<-paste("~",paste(c("fatty_liver",d$covariates,"(1 | family_id)"),collapse=" + "))
  cat(format(Sys.time())," LinDA global sensitivity: all 695 features\n");flush.console()
  fit<-MicrobiomeStat::linda(feature.dat=t(x$relative[d$index,,drop=FALSE]),meta.dat=meta,formula=formula,
    feature.dat.type="proportion",prev.filter=0,mean.abund.filter=0,max.abund.filter=0,
    is.winsor=TRUE,outlier.pct=.03,adaptive=FALSE,p.adj.method="none",
    n.cores=if(.Platform$OS.type=="windows")1L else cfg$ab_cores,verbose=TRUE)
  saveRDS(fit,file.path(out,"private/LinDA_complete.rds"))
  term<-grep("^fatty_liver",names(fit$output),value=TRUE)
  if(length(term)!=1L)stop("Unexpected LinDA disease coefficient")
  rr<-fit$output[[term]];critical<-qt(.975,rr$df)
  data.frame(feature_id=rownames(rr),estimate_log2=rr$log2FoldChange,SE=rr$lfcSE,df=rr$df,
    CI_low=rr$log2FoldChange-critical*rr$lfcSE,CI_high=rr$log2FoldChange+critical*rr$lfcSE,
    p_raw=rr$pvalue,status=ifelse(is.finite(rr$pvalue),"ESTIMATED","NOT_ESTIMABLE"))
}

fl_ab_complete <- function(x,results,linda,out) {
  full<-x$features;features<-full$feature_id
  for(name in c("primary","no_antibiotics","adult_core","adult_lifestyle","adult_BMI")){
    r<-results[[name]]$effects;j<-match(features,r$feature_id)
    for(v in c("estimate","CI_low","CI_high","p_raw","status"))full[[paste(name,v,sep="_")]]<-r[[v]][j]
  }
  generation<-results$generation$effects
  for(g in c("G0","G1","G2")){
    r<-generation[generation$term==paste0("case_",g),];j<-match(features,r$feature_id)
    for(v in c("estimate","CI_low","CI_high","p_raw","status"))full[[paste(g,v,sep="_")]]<-r[[v]][j]
  }
  h<-results$generation$heterogeneity;j<-match(features,h$feature_id)
  full$p_generation_difference_raw<-h$p_interaction_raw[j];full$generation_test_status<-h$status[j]
  j<-match(features,linda$feature_id)
  for(v in c("estimate_log2","CI_low","CI_high","p_raw","status"))full[[paste0("LinDA_",v)]]<-linda[[v]][j]
  for(g in c("Overall","G0","G1","G2")){
    a<-x$support[x$support$generation==g,];j<-match(features,a$feature_id)
    for(v in c("n","cases","detected","detected_cases","detected_controls","prevalence"))full[[paste(g,v,sep="_")]]<-a[[v]][j]
  }
  full$reference_log10_mean<-x$reference_mean[features];full$reference_log10_SD<-x$reference_sd[features]
  full$LinDA_direction_agreement<-ifelse(is.finite(full$primary_estimate)&is.finite(full$LinDA_estimate_log2),
    sign(full$primary_estimate)==sign(full$LinDA_estimate_log2),NA)
  ok<-startsWith(full$primary_status,"ESTIMATED")&is.finite(full$primary_p_raw)&
    full$Overall_detected_cases>=10&full$Overall_detected_controls>=10
  selected<-which(ok);selected<-head(selected[order(full$primary_p_raw[selected],features[selected])],12L)
  full$shown_in_figure<-seq_len(nrow(full))%in%selected
  stopifnot(nrow(full)==695L,!anyDuplicated(full$feature_id),length(selected)>0)
  fl_write_csv(full,file.path(out,"tables/SGB_abundance_complete_results_EN.csv"))
  fl_write_csv(full[selected,,drop=FALSE],file.path(out,"tables/SGB_abundance_figure_selection_EN.csv"))
  fl_write_csv(generation,file.path(out,"tables/SGB_generation_effects_EN.csv"))
  fl_write_csv(h,file.path(out,"tables/SGB_generation_difference_EN.csv"))
  fl_write_csv(linda,file.path(out,"tables/SGB_LinDA_sensitivity_EN.csv"))
  fl_write_csv(x$support,file.path(out,"tables/SGB_detection_support_EN.csv"))
  full
}

fl_ab_sensitivity_report <- function(r,out) {
  f<-r$complete;shown<-f$shown_in_figure
  stopifnot(nrow(f)==695L,sum(shown)==12L,!anyDuplicated(f$feature_id))
  scenarios<-c("primary","no_antibiotics","adult_core","adult_lifestyle","adult_BMI","LinDA")
  tab<-do.call(rbind,lapply(scenarios,function(v){
    b<-f[[paste0(v,if(v=="LinDA")"_estimate_log2"else"_estimate")]]
    p<-f[[paste0(v,"_p_raw")]];valid<-is.finite(b)&is.finite(p)
    same<-sign(b)==sign(f$primary_estimate)
    data.frame(analysis=v,estimable_SGBs=sum(valid),nominal_P_lt_0_05=sum(p<.05,na.rm=TRUE),
      negative_effects=sum(valid&b<0),positive_effects=sum(valid&b>0),
      displayed_SGBs_nominal_P_lt_0_05=sum(shown&p<.05,na.rm=TRUE),
      same_direction_as_primary=sum(valid&same,na.rm=TRUE),
      unit=if(v=="LinDA")"Bias-corrected log2 effect"else"SD of core-reference log10 relative abundance")
  }))
  fl_write_csv(tab,file.path(out,"tables/SGB_sensitivity_summary_EN.csv"))
  h<-sum(f$p_generation_difference_raw<.05,na.rm=TRUE)
  hs<-sum(shown&f$p_generation_difference_raw<.05,na.rm=TRUE)
  j<-which.min(f$primary_p_raw)
  format_p<-function(p)format(p,digits=4,scientific=TRUE)
  lines<-c("# 脂肪肝相关SGB丰度分析结果", "",
    "本轮分析使用基线三代均有成员的主家系；全龄核心完整样本1267人、232例脂肪肝、319个家系。G0/G1/G2分别465/315/487人、128/87/17例，保留未成年人。原始数据和既有临床/群落分析结果保留。", "",
    paste0("全部695个SGB进入分析，总体有",tab$nominal_P_lt_0_05[1],"个SGB的原始P<0.05。效应是脂肪肝组减去非脂肪肝组的标准化对数相对丰度差异，单位为全龄核心参考样本SD，不是OR或绝对菌量。"),
    sprintf("总体最小原始P对应%s：差异%.3f SD（95%%CI %.3f至%.3f），P=%s。完整分类注释保留在英文全表中；多重或未命名注释不强行改写为确定的单一物种。",
      f$official_id[j],f$primary_estimate[j],f$primary_CI_low[j],f$primary_CI_high[j],format_p(f$primary_p_raw[j])),"",
    paste0("全695个SGB中，",h,"个的2自由度代际差异检验原始P<0.05。按预定总体P和检出支持规则展示的12个SGB中，",hs,"个达到该代际差异阈值。单代显著、另一代不显著不能代替正式的代际差异检验。G2仅17例，不能视为独立验证。"),"",
    "## 主要调整", "",
    "代际、代内年龄自然样条、性别、调查地区、建库批次、测序池、log10输入矩阵总量、采便季节、冻存间隔、采便与超声时间差、近期抗生素；所有丰度模型包含主家系随机截距。常量或完全别名列的处理见audit中的映射。", "",
    "## 稳健性", "",
    "抗生素排除分析：1165人、210例。成人核心、生活方式及生活方式+BMI模型均为同一736人、158例。成人完整样本来自1004名核心成人，其中268人因生活方式协变量或BMI缺失未进入此敏感性分析；各变量缺失及各代筛选人数另表列出。本轮成人SGB结果未复用先前临床MI，不能声称已经消除完整病例选择偏倚。",
    paste0("LinDA在相同1267人、全部695个SGB上共同运行，使用同一核心协变量空间和家系随机截距；",tab$same_direction_as_primary[tab$analysis=="LinDA"],"个SGB的估计方向与主要模型一致。其系数为偏倚校正log2单位，与主模型SD单位不同，不能直接相减或解释为独立验证。"),
    "各模型名义P<0.05数量、方向和主图特征支持情况见SGB_sensitivity_summary_EN.csv。全龄与成人属于不同人群，只有三种同样本成人模型之间可直接描述增加协变量后的变化。", "",
    "## 解释范围", "",
    "按既有要求，所有公开P均为原始P。695个SGB的名义显著只用于探索，不能称为已控制多重比较的发现。SGB丰度关联及家系随机截距不能证明遗传、同住或菌株传播。主图只展示预定的12个特征，全表保留全部695个SGB及模型状态。", "",
    "图形：Fig06_SGB_abundance、Supp05_SGB_abundance_extensions、Supp06_SGB_composition_sensitivity；由19号独立脚本生成，PDF/PNG/TIFF/SVG。")
  writeLines(lines,file.path(out,"reports/Results_summary_zh.md"))
  invisible(tab)
}

fl_ab_run <- function(cfg,mode="--run") {
  x<-fl_ab_prepare(cfg);scenarios<-fl_ab_scenarios(x);designs<-lapply(scenarios,function(z)fl_ab_design(x,z))
  counts<-do.call(rbind,lapply(names(designs),function(v){d<-designs[[v]]
    data.frame(scenario=v,n=d$n,cases=d$cases,families=d$families,minors=sum(d$metadata$age<18),nuisance_columns=length(d$covariates))}))
  if(mode=="--check"){print(counts);cat("PASS: checksums, fixed reference scaling and all protected designs\n");return(invisible(counts))}
  resume<-Sys.getenv("FATTYLIVER_ABUNDANCE_RESUME_RUN",unset="")
  if(nzchar(resume)){
    if(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$",resume))stop("Invalid resume run ID")
    cfg$run_id<-resume;module<-file.path(cfg$project_root,"outputs",cfg$output_module)
    out<-file.path(module,"_running",resume)
    if(!dir.exists(out)||!readLines(file.path(out,"run_status.txt"),n=1)%in%c("RUNNING","CANCELLED_BY_USER"))stop("Only interrupted runs can resume")
    ctx<-list(cfg=cfg,module=module,run_dir=out)
  }else{ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir}
  for(v in c("tables","private","audit","reports","figures"))dir.create(file.path(out,v),showWarnings=FALSE)
  writeLines("RUNNING",file.path(out,"run_status.txt"))
  tryCatch({
    if(nzchar(resume)){
      old<-read.csv(file.path(out,"audit/input_manifest.csv"),check.names=FALSE)
      stopifnot(identical(as.character(old$path),x$source_manifest$path),identical(as.character(old$md5),x$source_manifest$md5))
    }
    fl_write_csv(x$source_manifest,file.path(out,"audit/input_manifest.csv"))
    fl_write_csv(counts,file.path(out,"tables/SGB_scenario_sample_counts_EN.csv"))
    fl_write_csv(x$sample_flow,file.path(out,"tables/SGB_sample_flow_EN.csv"))
    fl_write_csv(x$sample_counts,file.path(out,"tables/SGB_primary_generation_counts_EN.csv"))
    fl_ab_sample_audit(x,scenarios,out)
    saveRDS(list(metadata=x$metadata,mean=x$reference_mean,sd=x$reference_sd),file.path(out,"private/reference_data.rds"))
    fl_write_csv(x$knots,file.path(out,"audit/age_spline_knots.csv"))
    results<-list()
    # Outcome-blind small benchmark; exact completed features are reused below.
    start<-Sys.time();pilot<-fl_ab_stage(x,designs$primary,cfg,out,"primary","overall",seq_len(20L))
    fl_write_csv(data.frame(features=20,seconds=as.numeric(difftime(Sys.time(),start,units="secs"))),file.path(out,"audit/runtime_benchmark.csv"))
    results$primary<-fl_ab_stage(x,designs$primary,cfg,out,"primary")
    results$generation<-fl_ab_stage(x,designs$primary,cfg,out,"primary","generation")
    for(v in setdiff(names(designs),"primary"))results[[v]]<-fl_ab_stage(x,designs[[v]],cfg,out,v)
    linda<-fl_ab_linda(x,designs$primary,cfg,out)
    full<-fl_ab_complete(x,results,linda,out)
    diagnostics<-do.call(rbind,lapply(results,`[[`,"diagnostics"))
    fl_write_csv(diagnostics,file.path(out,"tables/SGB_model_diagnostics_EN.csv"))
    derived<-list(results=results,linda=linda,complete=full,counts=counts)
    saveRDS(derived,file.path(out,"private/derived_results.rds"))
    fl_ab_sensitivity_report(derived,out)
    summary<-data.frame(item=c("SGBs","Primary raw P < 0.05","Generation difference raw P < 0.05","Valid primary effects","Valid LinDA effects","LinDA direction agreement"),
      value=c(nrow(full),sum(full$primary_p_raw<.05,na.rm=TRUE),sum(full$p_generation_difference_raw<.05,na.rm=TRUE),
        sum(is.finite(full$primary_p_raw)),sum(is.finite(full$LinDA_p_raw)),mean(full$LinDA_direction_agreement,na.rm=TRUE)))
    fl_write_csv(summary,file.path(out,"tables/SGB_analysis_summary_EN.csv"))
    text<-c("# Standardized SGB abundance associations", "",
      "Primary population: 1267 participants from eligible three-generation primary families, including minors; 232 steatosis cases. All 695 SGBs retained.",
      "Outcome: log10 relative abundance with the existing fixed pseudocount, standardized once in the full core sample. Binary steatosis is not standardized. Effects are SD differences, not ORs or absolute microbial counts.",
      "Official MaAsLin 2 LM with family random intercepts (REML/lmerTest). Nuisance adjustment: generation, within-generation age splines, sex, survey centre, library batch, sequencing pool, log10 input-matrix total, season sine/cosine, recorded freezing interval, stool-ultrasound interval, recent antibiotics. Constant/aliased nuisance columns are recorded.",
      "Generation effects use three disease-by-generation indicators and a joint two-df Satterthwaite F test. G2 has 17 cases and is not an independent validation set.",
      "Adult extensions use the same 736 complete-case adults (158 cases) for core, lifestyle and lifestyle+BMI models. This is a complete-case sensitivity; missing adult lifestyle data were not imputed. Previous clinical imputations were not assumed compatible with these new outcomes. Adult selection may bias these estimates.",
      "LinDA uses all 695 proportions jointly with the same core design, official proportion-zero handling and default 3% winsorization. Its bias-corrected log2 effect has a different estimand/unit; direction agreement is sensitivity evidence, not independent validation.",
      "All public P values are raw. Nominal significance is exploratory and does not establish FDR-controlled discoveries. Primary-family membership is not verified pedigree or cohabitation; no transmission/mediation inference.",
      "Display selection: 12 smallest primary P among estimable features with >=10 detected cases and >=10 detected controls; ID tie-break. All features remain in the complete table.",
      "Sources: Zhou et al. Nat Commun (2026), doi:10.1038/s41467-026-71981-0; Mallick et al. PLOS Comput Biol (2021), doi:10.1371/journal.pcbi.1009442; Zhou et al. Genome Biol (2022), doi:10.1186/s13059-022-02655-5. Exact covariate coding and generation contrasts are cohort adaptations.")
    writeLines(text,file.path(out,"reports/Methods_EN.md"))
    for(f in c("scripts/18_sgb_abundance.R","R/sgb_abundance.R","config/sgb_abundance_v3_20261022.md"))file.copy(file.path(cfg$project_root,f),file.path(out,"audit",basename(f)))
    capture.output(sessionInfo(),file=file.path(out,"audit/sessionInfo.txt"))
    stopifnot(identical(unname(tools::md5sum(x$source_manifest$path)),x$source_manifest$md5))
    writeLines("COMPLETED",file.path(out,"run_status.txt"));final<-fl_output_complete_run(ctx)
    key<-file.path(cfg$project_root,"Results_表格与图片/Results_重点表格/tables/English/SGB_abundance")
    dir.create(key,recursive=TRUE,showWarnings=FALSE)
    file.copy(list.files(file.path(final,"tables"),full.names=TRUE),key,overwrite=TRUE)
    cat("COMPLETED: ",final,"\n",sep="");invisible(final)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(e)})
}
