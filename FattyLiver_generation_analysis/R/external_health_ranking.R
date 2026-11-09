# Rebuild of module 30: fixed external health-ranking scores, not old top/bottom sets.
# All definitions are inert until ehr_main() is called.
ehr_protocol <- 'Asnicar_fixed661_mapped_health_scores_v2_20261029'
ehr_scores <- c('Unweighted','Abundance_weighted')
ehr_dir <- function(root)file.path(root,'literature_reviews/external_rank_redesign_20261029')

ehr_compute <- function(cfg,out){
  src<-tf_source(cfg$project_root);x<-tf_load(src);m<-x$metadata
  paths<-c(src$files[1],file.path(ehr_dir(cfg$project_root),c('frozen_health_rank_full_precision.csv','public_score_reproduction_audit.csv','ZOE_Score_Public_Data_Fundamental_Analyses.py','upstream_commit.txt')),
    file.path(cfg$project_root,'literature_reviews/ABC_primary_20261027/Asnicar2026_Supplementary_Tables.xlsx'))
  hashes<-lr_hashes(paths);fl_write_csv(hashes,file.path(out,'audit/input_fingerprints.csv'))
  validation<-read.csv(paths[3]);stopifnot(nrow(validation)==2L,all(validation$N==4332),all(validation$max_absolute_error<1e-10))
  fl_write_csv(validation,file.path(out,'audit/published_score_reproduction.csv'))
  ranks<-read.csv(paths[2],stringsAsFactors=FALSE);stopifnot(nrow(ranks)==661L,!anyDuplicated(ranks$SGB))
  fixed<-1-2*(ranks$HEALTH_ranks-min(ranks$HEALTH_ranks))/diff(range(ranks$HEALTH_ranks))
  stopifnot(max(abs(fixed-ranks$fixed_weight))<1e-12,min(fixed)==-1,max(fixed)==1)
  ft<-x$features[match(rownames(x$relative),x$features$feature_id),]
  stopifnot(!anyNA(ft$feature_id),!anyDuplicated(ft$official_id))
  j<-match(ranks$SGB,ft$official_id);map<-ranks
  map$feature_id<-ft$feature_id[j];map$mapped<-!is.na(j)
  map$status<-ifelse(map$mapped,'Exact official identifier match','Unavailable: not interpreted as biological absence')
  map$observed_prevalence<-NA_real_;map$observed_mean_relative_abundance<-NA_real_
  map$observed_prevalence[map$mapped]<-rowMeans(x$relative>0)[j[map$mapped]]
  map$observed_mean_relative_abundance[map$mapped]<-rowMeans(x$relative)[j[map$mapped]]
  v3_table(map,'Table30_Frozen_661_SGB_mapping_EN',out)
  A<-x$relative[map$feature_id[map$mapped],,drop=FALSE];w<-map$fixed_weight[map$mapped]
  stopifnot(identical(colnames(A),as.character(m$sample_id)),all(A>=0),all(A<=1))
  m$Unweighted<-as.numeric(crossprod(w,A>0))
  m$Abundance_weighted<-as.numeric(crossprod(w,asin(sqrt(A))))
  m$mapped_SGBs_detected<-colSums(A>0);m$mapped_relative_abundance<-colSums(A)
  # Independent per-sample summation, plus invariance to species ordering.
  for(i in unique(c(1L,ncol(A),round(ncol(A)/2)))){
    stopifnot(abs(m$Unweighted[i]-sum(w[A[,i]>0]))<1e-10,
      abs(m$Abundance_weighted[i]-sum(w*asin(sqrt(A[,i]))))<1e-10)
  }
  stopifnot(max(abs(as.numeric(crossprod(rev(w),A[nrow(A):1,,drop=FALSE]>0))-m$Unweighted))<1e-10)
  coverage<-data.frame(published_SGBs=661L,exactly_mapped=sum(map$mapped),unavailable=sum(!map$mapped),identifier_coverage_percent=100*mean(map$mapped),
    available_group_identifiers=sum(map$mapped&grepl('_group$',map$SGB)),missing_group_identifiers=sum(!map$mapped&grepl('_group$',map$SGB)),
    median_observed_total_abundance_percent=100*median(m$mapped_relative_abundance),
    qualification='Partial-coverage score; missing taxa are unknown; no within-mapped reranking or abundance renormalization')
  v3_table(coverage,'Table30_External_ranking_coverage_EN',out)
  b1<-v3_cc(m,c('fatty_liver','family_id',v3_B1))
  adult<-!is.na(m$age)&m$age>=18&v3_cc(m,c('fatty_liver','family_id',v3_B2))
  source<-!is.na(m$age)&m$age>=18&v3_cc(m,c('fatty_liver','family_id','age','sex','BMI18'))
  specs<-list(B1=list(mask=b1,fields=v3_B1,spline=TRUE,ols=FALSE),
    Adult_B1_same_sample=list(mask=adult,fields=v3_B1,spline=TRUE,ols=FALSE),
    B2=list(mask=adult,fields=v3_B2,spline=TRUE,ols=FALSE),
    Source_age_sex_BMI_GEE=list(mask=source,fields=c('age','sex','BMI18'),spline=FALSE,ols=FALSE),
    Source_age_sex_BMI_OLS=list(mask=source,fields=c('age','sex','BMI18'),spline=FALSE,ols=TRUE))
  for(g in c('G0','G1','G2'))specs[[paste0('Generation_',g)]]<-list(mask=b1&m$generation==g,fields=setdiff(v3_B1,'generation'),spline=TRUE,ols=FALSE)
  rows<-list();counts<-list();designs<-list();residuals<-list()
  for(key in names(specs)){
    sp<-specs[[key]];z<-droplevels(m[sp$mask,]);z<-z[order(z$family_id,z$sample_id),]
    stopifnot(nrow(z)>0,all(c(0,1)%in%z$fatty_liver))
    des<-abc_design(z,sp$fields,age_spline=sp$spline);Z<-des$Z
    stopifnot(qr(cbind(z$fatty_liver,Z))$rank==ncol(Z)+1L)
    fl_write_csv(des$mapping,file.path(out,'audit',paste0(key,'_design_columns.csv')))
    designs[[key]]<-list(sample_ids=z$sample_id,fields=sp$fields,dropped=des$dropped,rank=ncol(Z))
    counts[[key]]<-v3_counts(z,key)
    tf_log('External health ranking: ',key,'; N=',nrow(z),'; cases=',sum(z$fatty_liver))
    for(score in ehr_scores){
      if(sp$ols){dd<-data.frame(y=z[[score]],exposure=z$fatty_liver,Z);fit<-lm(y~.-1,dd)
        co<-coef(summary(fit))['exposure',];ci<-confint(fit,'exposure')
        rr<-data.frame(beta=co[1],SE=co[2],CI_low=ci[1],CI_high=ci[2],P_raw=co[4],row.names=NULL)
      }else{ans<-abc_gee(z[[score]],z$fatty_liver,Z,z);rr<-ans$row;fit<-ans$fit}
      # Independent OLS coefficient check: working-independence Gaussian GEE
      # must agree in point estimate, but uses family-robust standard errors.
      ols<-lm.fit(cbind(exposure=z$fatty_liver,Z),z[[score]])
      stopifnot(abs(ols$coefficients[1]-rr$beta)<1e-7)
      source_d<-source_se<-NA_real_
      if(sp$ols){n0<-sum(z$fatty_liver==0);n1<-sum(z$fatty_liver==1);N<-n0+n1
        source_d<-(rr$beta/rr$SE)*N/(sqrt(n0*n1)*sqrt(N-2))
        source_se<-sqrt(((N-1)/(N-3))*((4/N)*(1+source_d^2/8)))}
      rows[[paste(key,score)]]<-cbind(model=key,score=score,N=nrow(z),cases=sum(z$fatty_liver),families=length(unique(z$family_id)),
        minors=sum(z$age<18),nuisance_rank=ncol(Z),rr,source_code_d=source_d,source_code_d_SE=source_se,
        inference=if(sp$ols)'Source OLS comparator; ignores family dependence'else'Gaussian GEE; family-robust sandwich SE',
        effect_unit='Native partial-score points; steatosis minus no steatosis')
      saveRDS(fit,file.path(out,'private',paste0(key,'_',score,'_fit.rds')))
      if(key=='B1')residuals[[score]]<-data.frame(sample_id=z$sample_id,fitted=fitted(fit),residual=residuals(fit),score=score)
    }
  }
  stopifnot(identical(designs$Adult_B1_same_sample$sample_ids,designs$B2$sample_ids))
  t<-ov_bind(rows);t$FDR_family<-ifelse(grepl('^Generation_',t$model),'Six_generation_score_tests',paste0(t$model,'_two_score_tests'))
  t$P_BH<-ave(t$P_raw,t$FDR_family,FUN=function(p)p.adjust(p,'BH'));t$BH_lt_005<-t$P_BH<.05
  t$source_code_d_CI_low<-t$source_code_d-1.96*t$source_code_d_SE;t$source_code_d_CI_high<-t$source_code_d+1.96*t$source_code_d_SE
  v3_table(t,'Table30_All_health_score_associations_EN',out)
  v3_table(t[t$model%in%c('B1','Adult_B1_same_sample','B2'),],'Table30_Main_health_score_associations_EN',out)
  v3_table(ov_bind(counts),'Table30_Analysis_populations_EN',out)
  desc<-ov_bind(lapply(c('Overall','G0','G1','G2'),function(g)ov_bind(lapply(ehr_scores,function(score)ov_bind(lapply(0:1,function(e){
    z<-m[b1&m$fatty_liver==e&(g=='Overall'|m$generation==g),];y<-z[[score]]
    data.frame(group=g,score=score,steatosis=e,N=nrow(z),mean=mean(y),SD=sd(y),median=median(y),Q1=unname(quantile(y,.25)),Q3=unname(quantile(y,.75))) }))))))
  v3_table(desc,'Table30_Health_score_distributions_EN',out)
  m$included_B1<-b1;m$included_adult_common<-adult;m$included_source_adult<-source
  saveRDS(list(protocol=ehr_protocol,metadata=m,mapping=map,description='Partial fixed-rank scores; original source data unchanged'),file.path(out,'private/new_external_health_scores.rds'))
  saveRDS(designs,file.path(out,'private/model_designs.rds'));saveRDS(residuals,file.path(out,'private/B1_model_diagnostics.rds'))
  lr_check_hashes(hashes)
  list(protocol=ehr_protocol,results=t,summary=desc,coverage=coverage,mapping=map,metadata=m,hashes=hashes,designs=designs)
}

ehr_figures<-function(r,cfg,out){
  bin<-c('No steatosis'='#79CAFB','Steatosis'='#F7A24F')
  d<-r$metadata[r$metadata$included_B1,];d$status<-factor(ifelse(d$fatty_liver==1,'Steatosis','No steatosis'),levels=names(bin))
  titles<-c(Unweighted='Unweighted health score',Abundance_weighted='Abundance-weighted health score')
  # ECDF is an unsmoothed description of every observed value; no independent-
  # subject KS test is applied to clustered family observations.
  ecdf_plot<-function(score)ggplot2::ggplot(d,ggplot2::aes(x=.data[[score]],colour=status,linetype=status))+
    ggplot2::stat_ecdf(geom='step',linewidth=.55,pad=FALSE)+ggplot2::scale_colour_manual(values=bin)+
    ggplot2::scale_linetype_manual(values=c('solid','22'))+
    ggplot2::scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1),labels=c('0','25','50','75','100'))+
    ggplot2::labs(title=titles[[score]],x='Partial score (higher = more favourable)',y='Cumulative participants (%)',colour=NULL,linetype=NULL)+
    lr_theme()+ggplot2::theme(legend.position='bottom')
  forest_plot<-function(score){
    z<-r$results[r$results$score==score&r$results$model%in%c('B1','Adult_B1_same_sample','B2'),]
    z$label<-factor(z$model,levels=c('B2','Adult_B1_same_sample','B1'),labels=c('Adult B2','Adult B1','All-age B1'))
    ggplot2::ggplot(z,ggplot2::aes(beta,label))+
      ggplot2::geom_vline(xintercept=0,linetype=2,colour='#888888',linewidth=.25)+
      ggplot2::geom_segment(ggplot2::aes(x=CI_low,xend=CI_high,yend=label),linewidth=.45,colour='black')+
      ggplot2::geom_point(shape=21,fill='#777777',colour='black',stroke=.25,size=2)+
      ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=.6))+
      ggplot2::labs(title=titles[[score]],x='Adjusted score difference (95% CI)',y=NULL)+lr_theme()
  }
  one<-lr_export(list(ecdf_plot('Unweighted'),ecdf_plot('Abundance_weighted'),forest_plot('Unweighted'),forest_plot('Abundance_weighted')),
    'Figure30_External_health_scores',out,cfg,width=183,height=140,nrow=2,uniform_y_axes=TRUE)
  gen_plot<-function(score){
    z<-r$results[r$results$score==score&grepl('^Generation_',r$results$model),]
    z$generation<-factor(sub('Generation_','',z$model),levels=c('G2','G1','G0'))
    cols<-c(G0='#F1B656',G1='#397FC7',G2='#040676')
    ggplot2::ggplot(z,ggplot2::aes(beta,generation))+
      ggplot2::geom_vline(xintercept=0,linetype=2,colour='#888888',linewidth=.25)+
      ggplot2::geom_segment(ggplot2::aes(x=CI_low,xend=CI_high,yend=generation,colour=generation),linewidth=.5)+
      ggplot2::geom_point(ggplot2::aes(fill=generation),shape=21,colour='black',stroke=.25,size=2)+
      ggplot2::scale_colour_manual(values=cols,guide='none')+ggplot2::scale_fill_manual(values=cols,guide='none')+
      ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=.6))+
      ggplot2::labs(title=titles[[score]],x='Adjusted score difference (95% CI)',y=NULL)+lr_theme()
  }
  two<-lr_export(lapply(ehr_scores,gen_plot),'FigureS30_Generation_health_scores',out,cfg,width=183,height=85,nrow=1,uniform_y_axes=TRUE)
  v3_table(ov_bind(list(one,two)),'Figure30_Export_manifest_EN',out)
  # Every plotting value is also available in distribution and effect tables;
  # participant-level records remain private rather than entering the gallery.
  writeLines(c('# Figure legends','',
    '**Figure 30. External health-ranking scores and hepatic steatosis.** (a,b) Empirical cumulative distributions of unweighted and abundance-weighted scores in the all-age B1 sample. Each curve includes all eligible participants; higher scores indicate a more favourable external rank profile. (c,d) Adjusted steatosis-minus-no-steatosis score differences and family-robust 95% confidence intervals. All-age B1 includes generation, generation-specific age splines, sex, library batch, sequencing pool, sampling season and recent antibiotic use. Adult B1 and B2 use identical complete-case adults; B2 additionally includes smoking, alcohol, diet score, physical activity, energy intake and continuous BMI. Dashed vertical lines mark zero. Raw P and BH q values (two score tests per model) and sample counts are reported in Table30. Scores use exact matches to 578 of 661 externally ranked SGBs; unavailable SGBs are not treated as absent. These are partial-coverage scores, not complete ZOE scores. Fixed weights and AST abundance weighting follow Asnicar et al., Nature 650, 450-458 (2026), Fig. 3f. Native score effects, family GEE and cumulative distributions are adaptations; no cross-cohort meta-analysis was performed.',
    '**Supplementary Figure 30. Generation-specific associations of partial external health scores with hepatic steatosis.** Points and horizontal lines show within-generation adjusted score differences and family-robust 95% confidence intervals. Models use B1 covariates with within-generation age splines and without a generation main effect. BH correction covers all six generation-by-score tests. Differences between generation point estimates do not establish effect heterogeneity. Exact sample counts, P and q values are in Table30.'),file.path(out,'reports/Figure_legends_EN.md'))
}

ehr_reports<-function(r,cfg,out){
  method<-file.path(cfg$project_root,'analysis_plan/external_health_ranking_redesign_20261029.md')
  file.copy(method,file.path(out,'reports/Methods_and_references_EN.md'),overwrite=TRUE)
  t<-r$results[r$results$model%in%c('B1','Adult_B1_same_sample','B2'),]
  writeLines(c('# 外部健康排名分析：重做结果','',
    sprintf('固定外部661个SGB排名；精确匹配%d个，%d个无法匹配。仅计算可匹配菌的贡献，不能视作完整ZOE评分。',r$coverage$exactly_mapped,r$coverage$unavailable),
    '评分越高，外部排名所对应的健康有利特征越多。差值为脂肪肝者减去非脂肪肝者。',
    '',capture.output(print(t[,c('model','score','N','cases','beta','CI_low','CI_high','P_raw','P_BH')],row.names=FALSE)),
    '', '成人B1和B2使用同一批完整样本；B2增加生活方式及BMI。所有主要置信区间按主家系聚类处理。',
    '旧版最有利/最不利50菌集合的结果不再采纳，且未用于这次重新计算。原始临床和菌群数据保留。',
    '结果是横断面关联，不能表示疾病风险预测、致病菌鉴定或菌株传播。未成年人属于原文成人背景之外的队列扩展。'),file.path(out,'reports/Results_summary_zh.md'),useBytes=TRUE)
  capture.output(sessionInfo(),file=file.path(out,'audit/sessionInfo.txt'))
}

ehr_main<-function(args=character(),root=Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')){
  if(getRversion()<'4.5'){
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts/30_external_health_sets.R')),vapply(args,shQuote,character(1))))
    if(status!=0L)stop('Project R external-ranking entry failed');return(invisible(NULL))
  }
  sys.source(file.path(root,'R/literature_ABC_20261028.R'),envir=.GlobalEnv);cfg<-lr_setup(root)
  cfg$output_module<-'30_External_health_sets';cfg$run_id<-'external_health_ranking_v2_20261029';cfg$run_id_explicit<-FALSE
  if(any(!args%in%c('--run','--check','--figures-only')))stop('Allowed: --run --check --figures-only')
  if('--check'%in%args){tf_load(tf_source(root));stopifnot(file.exists(file.path(ehr_dir(root),'frozen_health_rank_full_precision.csv')));cat('Input/dependency check passed\n');return(invisible(NULL))}
  module<-file.path(root,'outputs',cfg$output_module)
  if('--figures-only'%in%args){
    # Prefer saved calculations from this protocol, including an export-failed run.
    paths<-unlist(lapply(c('_running','history','current'),function(s)list.dirs(file.path(module,s),recursive=FALSE,full.names=TRUE)))
    paths<-paths[grepl('external_health_ranking_v2_',basename(paths))&file.exists(file.path(paths,'private/analysis_result.rds'))]
    if(!length(paths))stop('No saved v2 calculation to render')
    out<-paths[which.max(file.info(file.path(paths,'private/analysis_result.rds'))$mtime)]
    r<-readRDS(file.path(out,'private/analysis_result.rds'));stopifnot(identical(r$protocol,ehr_protocol));lr_check_hashes(r$hashes)
    ehr_figures(r,cfg,out);ehr_reports(r,cfg,out)
    if(basename(dirname(out))!='current'){
      newpath<-file.path(module,'_running',basename(out));if(out!=newpath){stopifnot(!dir.exists(newpath),file.rename(out,newpath));out<-newpath}
      cfg$run_id<-basename(out);ctx<-list(cfg=cfg,module=module,run_dir=out)
      writeLines(c('COMPLETED',ehr_protocol),file.path(out,'run_status.txt'));out<-fl_output_complete_run(ctx)
      writeLines(basename(out),file.path(module,'LATEST_RUN.txt'))
    }
    fl_publish_organized_results(cfg);tf_log('Completed saved-score rendering: ',out);return(invisible(out))
  }
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(p in c('private','tables','figures','reports','audit'))dir.create(file.path(out,p),showWarnings=FALSE)
  writeLines(c('RUNNING',ehr_protocol),file.path(out,'run_status.txt'))
  tryCatch({r<-ehr_compute(cfg,out);saveRDS(r,file.path(out,'private/analysis_result.rds'));ehr_reports(r,cfg,out);ehr_figures(r,cfg,out);lr_check_hashes(r$hashes)
    writeLines(c('COMPLETED',ehr_protocol),file.path(out,'run_status.txt'));target<-fl_output_complete_run(ctx)
    writeLines(basename(target),file.path(module,'LATEST_RUN.txt'));fl_publish_organized_results(cfg);tf_log('COMPLETED: ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
