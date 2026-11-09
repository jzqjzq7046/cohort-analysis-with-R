# Family-aware, all-adult revision of the exploratory Fig. 6 host map.
f6_protocol <- 'fig6_all_adults_family_LMM_20261105_v1'
f6_transform <- function(v) {
  ok <- is.finite(v); z <- rep(NA_real_,length(v)); n <- sum(ok)
  if(n<2L || length(unique(v[ok]))<2L) return(z)
  z[ok] <- as.numeric(scale(qnorm((rank(v[ok],ties.method='average')-.5)/n)))
  z
}
f6_fit <- function(id,d,X,candidate_ids,conditional,keep_model=FALSE) {
  dd <- d; dd$microbe <- X[id,dd$sample_id]
  fixed <- c('microbe','age_z','sex',if(dd$trait_id[1]!='BMI18')'BMI_z',if(conditional)'fatty_liver')
  form <- as.formula(paste('y ~',paste(fixed,collapse=' + '),'+ (1 | family_id)'))
  warnings <- character(); retries <- FALSE
  result <- tryCatch({
    fit <- withCallingHandlers(lmerTest::lmer(form,data=dd,REML=TRUE,
      control=lme4::lmerControl(optimizer='nloptwrap',check.conv.singular='ignore')),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')},
      message=function(m)invokeRestart('muffleMessage'))
    conv <- fit@optinfo$conv$lme4$messages
    if(length(conv)>0L || !is.null(fit@optinfo$conv$opt) && any(fit@optinfo$conv$opt!=0)) {
      retries <- TRUE
      fit <- withCallingHandlers(lmerTest::lmer(form,data=dd,REML=TRUE,
        control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000),check.conv.singular='ignore')),
        warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')},
        message=function(m)invokeRestart('muffleMessage'))
      conv <- fit@optinfo$conv$lme4$messages
    }
    b <- coef(summary(fit,ddf='Satterthwaite'))['microbe',]
    vc <- as.data.frame(lme4::VarCorr(fit)); fv<-vc$vcov[vc$grp=='family_id'];rv<-sigma(fit)^2
    converged <- !length(conv) && (is.null(fit@optinfo$conv$opt)||all(fit@optinfo$conv$opt==0))
    out <- data.frame(feature_id=id,N=nrow(dd),families=nlevels(dd$family_id),
      beta=unname(b['Estimate']),SE=unname(b['Std. Error']),df=unname(b['df']),
      CI_low=unname(b['Estimate']-qt(.975,b['df'])*b['Std. Error']),
      CI_high=unname(b['Estimate']+qt(.975,b['df'])*b['Std. Error']),P_raw=unname(b['Pr(>|t|)']),
      family_variance=fv,residual_variance=rv,ICC=fv/(fv+rv),
      singular=lme4::isSingular(fit,tol=1e-4),converged=converged,retried=retries,
      warning=paste(unique(c(warnings,conv)),collapse=' | '),
      status=if(converged)'ESTIMATED'else'NOT_ESTIMABLE: convergence failure',
      formula=paste(deparse(form),collapse=' '),stringsAsFactors=FALSE)
    diag <- NULL
    if(id %in% candidate_ids) {
      rr<-residuals(fit)/sigma(fit);ff<-fitted(fit)
      diag<-data.frame(feature_id=id,sample_id=dd$sample_id,trait_id=dd$trait_id[1],
        residual_z=rr,fitted=ff,observed=dd$y,microbe=dd$microbe,age_z=dd$age_z,BMI_z=dd$BMI_z)
    }
    list(row=out,diagnostics=diag,model=if(keep_model)fit else NULL)
  },error=function(e)list(row=data.frame(feature_id=id,N=nrow(dd),families=nlevels(dd$family_id),
    beta=NA_real_,SE=NA_real_,df=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P_raw=NA_real_,
    family_variance=NA_real_,residual_variance=NA_real_,ICC=NA_real_,singular=NA,converged=FALSE,
    retried=retries,warning=conditionMessage(e),status='NOT_ESTIMABLE: fit error',formula=paste(deparse(form),collapse=' ')),
    diagnostics=NULL,model=NULL))
  result
}
f6_prepare <- function(root,out) {
  parent<-fc_adopted(root,'50_Family_candidate_characterization')
  old<-readRDS(file.path(parent,'private/analysis_result.rds'))
  ade_protect(old$protected)
  x<-readRDS(old$protected$path[1]);m<-x$metadata
  ok<-!is.na(m$age)&m$age>=18&v3_cc(m,c('family_id','age','sex','BMI18'))&m$BMI18>0
  m<-droplevels(m[ok,]);m<-m[order(m$family_id,m$sample_id),]
  stopifnot(!anyDuplicated(m$participant_id),!anyDuplicated(m$sample_id))
  A<-x$relative[,match(m$sample_id,colnames(x$relative)),drop=FALSE]
  stopifnot(identical(colnames(A),m$sample_id),all(is.finite(A)),all(A>=0))
  prevalence<-rowMeans(A>0);ids<-rownames(A)[prevalence>.20]
  X<-t(vapply(ids,function(id)f6_transform(as.numeric(A[id,])),numeric(ncol(A))))
  dimnames(X)<-list(ids,m$sample_id)
  m$age_z<-as.numeric(scale(m$age));m$BMI_z<-as.numeric(scale(m$BMI18))
  filter<-data.frame(feature_id=rownames(A),prevalence=prevalence,eligible=prevalence>.20,
    candidate=rownames(A)%in%old$candidates$feature_id,rule='Detection >20% of all eligible adults before trait missingness')
  originals<-c(old$protected$path,file.path(parent,'private/analysis_result.rds'),
    list.files(file.path(parent,'figures'),full.names=TRUE),
    file.path(root,'11.4 fatty liver .pptx'))
  protected<-data.frame(path=originals,md5=unname(tools::md5sum(originals)))
  stopifnot(!anyNA(protected$md5))
  private<-list(metadata=m,abundance_transformed=X,relative_candidate=A[old$candidates$feature_id,,drop=FALSE],
    fields=old$host$fields,filter=filter,candidates=old$candidates,old_host=old$host,
    protected=protected,old_parent=parent,baseline_N=nrow(x$metadata),protocol=f6_protocol)
  saveRDS(private,file.path(out,'private/input_population.rds'))
  fl_write_csv(m[,c('sample_id','participant_id','family_id','generation','age','sex','BMI18')],
    file.path(out,'private/eligible_adult_roster.csv'))
  private
}
f6_analyze <- function(inp,out,workers=4L) {
  m<-inp$metadata;X<-inp$abundance_transformed;ids<-rownames(X);fields<-inp$fields
  cl<-parallel::makePSOCKcluster(workers);on.exit(parallel::stopCluster(cl),add=TRUE)
  libs<-.libPaths();parallel::clusterCall(cl,function(p){.libPaths(p);Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1');NULL},libs)
  parallel::clusterExport(cl,c('f6_fit'),envir=.GlobalEnv)
  candidate_ids<-inp$candidates$feature_id;parallel::clusterExport(cl,c('X','candidate_ids'),envir=environment())
  rows<-list();diagnostics<-list();missing<-list();nblock<-0L
  for(mode in c('Primary','Steatosis-conditional')) for(field in names(fields)) {
    conditional<-mode=='Steatosis-conditional';nblock<-nblock+1L
    checkpoint<-file.path(out,'private',paste0('block_',mode,'_',field,'.rds'))
    y<-m[[field]];lower<-if(field%in%c('ALAT18','ASAT18'))0 else .Machine$double.eps
    invalid<-!is.na(y)&(!is.finite(y)|y<lower);y[invalid]<-NA_real_
    yy<-f6_transform(y);ok<-is.finite(yy)&if(conditional)!is.na(m$fatty_liver) else TRUE
    d<-droplevels(m[ok,]);d$y<-yy[ok];d$trait_id<-field
    roster<-data.frame(sample_id=m$sample_id,included=ok,trait_observed=is.finite(yy),
      missing_steatosis=is.na(m$fatty_liver),trait=field,mode=mode)
    fl_write_csv(roster,file.path(out,'private',paste0('roster_',mode,'_',field,'.csv')))
    missing[[nblock]]<-data.frame(mode=mode,phenotype=field,phenotype_label=unname(fields[field]),
      eligible_adults=nrow(m),observed_trait=sum(is.finite(yy)),missing_or_invalid_trait=sum(!is.finite(yy)),
      invalid_trait=sum(invalid),N=nrow(d),families=nlevels(d$family_id),
      singleton_families=sum(table(d$family_id)==1),min_family=min(table(d$family_id)),max_family=max(table(d$family_id)))
    signature<-digest::digest(list(protocol=f6_protocol,ids=ids,d=d,X=X),algo='sha256')
    if(file.exists(checkpoint)) {
      block<-readRDS(checkpoint);stopifnot(identical(block$signature,signature))
      cat('RESUMED',mode,field,'\n');flush.console()
    } else {
      cat(format(Sys.time()),'FITTING',mode,field,'N=',nrow(d),'SGBs=',length(ids),'\n');flush.console()
      parallel::clusterExport(cl,c('d','conditional'),envir=environment())
      fits<-parallel::parLapplyLB(cl,ids,function(id)f6_fit(id,d,X,candidate_ids,conditional))
      block<-list(signature=signature,rows=ov_bind(lapply(fits,`[[`,'row')),
        diagnostics=ov_bind(lapply(fits,`[[`,'diagnostics')))
      saveRDS(block,checkpoint)
      cat('SAVED',mode,field,'estimated=',sum(block$rows$status=='ESTIMATED'),
        'singular=',sum(block$rows$singular,na.rm=TRUE),'\n');flush.console()
    }
    z<-block$rows;z$mode<-mode;z$phenotype<-field;z$phenotype_label<-unname(fields[field]);rows[[nblock]]<-z
    dg<-block$diagnostics;if(nrow(dg)){dg$mode<-mode;diagnostics[[nblock]]<-dg}
  }
  tests<-ov_bind(rows);tests$P_BH<-NA_real_
  for(mode in unique(tests$mode)) {
    ix<-which(tests$mode==mode);good<-is.finite(tests$P_raw[ix])&tests$status[ix]=='ESTIMATED'
    pp<-rep(1,length(ix));pp[good]<-tests$P_raw[ix][good]
    tests$P_BH[ix[good]]<-p.adjust(pp,'BH')[good]
  }
  tests$candidate<-tests$feature_id%in%candidate_ids
  tests$FDR_family<-'All eligible background SGBs x nine traits, separately within each panel'
  tests$effect_unit<-'SD of inverse-normal transformed trait per SD of inverse-normal transformed SGB'
  display<-expand.grid(feature_id=candidate_ids,phenotype=names(fields),mode=unique(tests$mode),stringsAsFactors=FALSE)
  display<-merge(display,tests,by=c('feature_id','phenotype','mode'),all.x=TRUE,sort=FALSE)
  display$status[is.na(display$status)]<-'NOT_TESTED: prevalence <=20%'
  display$strict<-is.finite(display$P_BH)&display$P_BH<.05
  display$exploratory<-is.finite(display$P_BH)&display$P_BH>=.05&display$P_BH<.10
  list(tests=tests,display=display,flow=ov_bind(missing),diagnostics=ov_bind(diagnostics),protocol=f6_protocol)
}
f6_tables_reports <- function(r,inp,out) {
  tests<-r$tests;d<-r$display
  counts<-ov_bind(lapply(c('Primary','Steatosis-conditional'),function(mode){z<-d[d$mode==mode,]
    data.frame(mode=mode,displayed=nrow(z),estimated=sum(z$status=='ESTIMATED'),
      candidate_links_q005=sum(z$strict),candidate_links_q010=sum(z$strict|z$exploratory))}))
  for(k in c('tests','display','flow'))v3_table(r[[k]],paste0('Fig6_',k,'_EN'),out)
  v3_table(inp$filter,'Fig6_background_filter_EN',out);v3_table(inp$candidates,'Fig6_frozen_candidates_EN',out)
  v3_table(counts,'Fig6_candidate_association_counts_EN',out)
  v3_table(d[d$strict|d$exploratory,],'Fig6_candidate_links_q010_EN',out)
  qcheck<-ov_bind(lapply(unique(tests$mode),function(mode){z<-tests[tests$mode==mode,];good<-z$status=='ESTIMATED'&is.finite(z$P_raw)
    pp<-rep(1,nrow(z));pp[good]<-z$P_raw[good];delta<-max(abs(z$P_BH[good]-p.adjust(pp,'BH')[good]))
    stopifnot(delta<1e-12,nrow(z)==sum(inp$filter$eligible)*9)
    data.frame(mode=mode,planned=nrow(z),estimated=sum(good),failed=sum(!good),max_BH_difference=delta,
      singular=sum(z$singular,na.rm=TRUE),nonconverged=sum(!z$converged))}))
  v3_table(qcheck,'Fig6_model_validation_EN',out)
  stopifnot(all(tests$CI_low[tests$status=='ESTIMATED']<=tests$beta[tests$status=='ESTIMATED']),
    all(tests$CI_high[tests$status=='ESTIMATED']>=tests$beta[tests$status=='ESTIMATED']),
    nrow(d)==180L,all(is.finite(inp$abundance_transformed)),
    max(abs(rowMeans(inp$abundance_transformed)))<1e-10)
  old<-inp$old_host$display;old$mode<-ifelse(old$mode=='Source-style','Primary',old$mode)
  comparison<-merge(d,old[,c('feature_id','phenotype','mode','rho','P_raw','P_BH','N')],
    by=c('feature_id','phenotype','mode'),all=TRUE,suffixes=c('_LMM','_one_per_family'))
  v3_table(comparison,'Fig6_prior_one_per_family_comparison_EN',out)
  dg<-r$diagnostics
  diagnostics<-ov_bind(lapply(split(dg,list(dg$mode,dg$trait_id,dg$feature_id),drop=TRUE),function(z){
    qq<-qnorm((seq_len(nrow(z))-.5)/nrow(z));a<-z$residual_z;ss<-sd(a)
    data.frame(mode=z$mode[1],phenotype=z$trait_id[1],feature_id=z$feature_id[1],N=nrow(z),
      residual_skewness=mean((a-mean(a))^3)/ss^3,residual_excess_kurtosis=mean((a-mean(a))^4)/ss^4-3,
      residual_normal_QQ_correlation=cor(sort(a),qq),max_absolute_residual=max(abs(a)),
      abs_residual_fitted_spearman=cor(abs(a),z$fitted,method='spearman'),
      residual_age_squared_correlation=cor(a,z$age_z^2),
      residual_BMI_squared_correlation=if(z$trait_id[1]=='BMI18')NA_real_ else cor(a,z$BMI_z^2))
  }))
  v3_table(diagnostics,'Fig6_candidate_model_diagnostics_EN',out)
  saveRDS(dg,file.path(out,'private/candidate_residual_diagnostics.rds'))
  nrange<-range(r$flow$N);frange<-range(r$flow$families);bg<-sum(inp$filter$eligible)
  en<-c('Figure 6 | Candidate SGBs and metabolic traits in adults, accounting for family clustering.',
    sprintf('All %s eligible adults from %s primary families were retained; model-specific N=%s-%s and families=%s-%s, depending on observed measurements.',
      nrow(inp$metadata),nlevels(inp$metadata$family_id),nrange[1],nrange[2],frange[1],frange[2]),
    'Rows show the ten frozen adult Model 4 discovery candidates (original BH q<0.10); columns show nine traits. (a) Linear mixed models adjust for age, sex and BMI, omitting BMI adjustment when BMI is the outcome. (b) Models additionally adjust for ultrasound steatosis. Both panels include a primary-family random intercept.',
    'Trait and SGB abundances undergo average-rank inverse-normal transformation followed by standardization; zeros are retained as tied observations. Colour shows adjusted beta, in SD of transformed trait per SD of transformed SGB, on one shared symmetric scale. These coefficients are not Spearman correlations or abundance ratios.',
    sprintf('Two-sided Satterthwaite t tests and t-based 95%% CIs are used. BH correction includes all %s background SGBs detected in >20%% of eligible adults x nine traits (%s planned tests), separately in each panel. * q<0.05; + 0.05<=q<0.10; x not tested/estimable. Full effects, CIs, P/q, per-model participant/family counts and diagnostic flags are tabulated.',bg,bg*9),
    'Per-trait complete cases are analyzed without imputation. This post hoc exploratory analysis is not independent validation, a purely within-family effect, or causal evidence. Primary-family labels do not establish exact kinship. Generation and technical/lifestyle covariates are not included in this reduced host-characterization model.')
  writeLines(en,file.path(out,'reports/Fig6_legend_EN.txt'),useBytes=TRUE)
  zh<-c('Fig. 6｜候选SGB与代谢指标的家系混合模型关联',
    sprintf('保留全部%s名合格成人、%s个主家系，不再每家抽1人。各模型实际N为%s–%s，家系数为%s–%s；按指标及协变量完整情况纳入，不插补。',nrow(inp$metadata),nlevels(inp$metadata$family_id),nrange[1],nrange[2],frange[1],frange[2]),
    '展示原成人Model 4发现的10个候选SGB与9项指标。a图调整年龄、性别和BMI，BMI自身分析不调整BMI；b图额外调整超声脂肪变。两图均加入主家系随机截距。固定协变量沿用原图的简化宿主表型模型，未额外调整世代和技术/生活方式因素。',
    '菌丰度和指标均采用平均秩逆正态转换后标准化；零丰度保留。颜色表示每增加1个标准差的转换后菌丰度，对应转换后指标的标准差变化β，不再是偏Spearman ρ。',
    sprintf('每个面板分别对%s个检出率>20%%的背景SGB×9项指标（%s项计划检验）进行BH校正；* q<0.05，+ 0.05≤q<0.10，x表示未检验/不可估计。完整β、95%%CI、P/q及各模型人数与家系数见表。',bg,bg*9),
    sprintf('候选菌关联：a图%s条q<0.05，累计%s条q<0.10；b图%s条q<0.05，累计%s条q<0.10。',counts$candidate_links_q005[1],counts$candidate_links_q010[1],counts$candidate_links_q005[2],counts$candidate_links_q010[2]),
    '这是同一队列中的事后探索性关联，不是独立验证或因果/中介证据。家系随机截距处理簇内相关，不等同于纯家系内比较，也不代表已知精确遗传亲缘关系。未达到q阈值不等同于无关联。',
    '原每家1人结果保留为单独敏感性比较；本次样本、转换及推断方法共同变化，不能将新旧差异归结为人数增加。')
  writeLines(zh,file.path(out,'reports/Fig6_legend_and_results_ZH.txt'),useBytes=TRUE)
  writeLines(c('# Methods, diagnostics and limitations',en,'',capture.output(print(qcheck,row.names=FALSE)),
    '', 'Convergence and zero variance flags are retained in the complete model table. A singular family intercept is a boundary estimate, not a failed fit. Complete-case missingness and the linear covariate/random-intercept assumptions limit inference. Rank transforms address marginal skewness, not every residual or compositional-data assumption. Associations are with relative abundance and cannot establish absolute microbial load effects.',
    '', 'Applied precedents: Beaumont et al. Genome Biology 2016 doi:10.1186/s13059-016-1052-7; Zierer et al. Nature Genetics 2018 doi:10.1038/s41588-018-0135-7; Visconti et al. Nature Communications 2019 doi:10.1038/s41467-019-12476-z. All are TwinsUK studies; their sampling/kinship models are not literally copied.',
    'Reporting: Althouse et al. Circulation 2021 doi:10.1161/CIRCULATIONAHA.121.055393. Project-specific choices are documented in the adopted plan.'),
    file.path(out,'reports/Methods_results_limits_EN.md'),useBytes=TRUE)
  counts
}
f6_figures <- function(r,inp,out,cfg) {
  ids<-inp$candidates$feature_id
  labels<-setNames(vapply(inp$candidates$plot_label,function(v)paste(strwrap(v,34),collapse='\n'),character(1)),ids)
  d<-r$display;d$row<-factor(d$feature_id,levels=rev(ids),labels=rev(labels))
  d$beta[d$status!='ESTIMATED']<-NA_real_
  d$trait<-factor(d$phenotype,levels=names(inp$fields),labels=unname(inp$fields))
  d$mark<-ifelse(d$strict,'*',ifelse(d$exploratory,'+',ifelse(d$status!='ESTIMATED','x','')))
  lim<-max(.05,ceiling(max(abs(d$beta),na.rm=TRUE)/.05)*.05)
  heat<-lapply(seq_along(c('Primary','Steatosis-conditional')),function(i){
    z<-d[d$mode==c('Primary','Steatosis-conditional')[i],]
    g<-ggplot2::ggplot(z,ggplot2::aes(trait,row,fill=beta))+
      ggplot2::geom_tile(colour='white',linewidth=.35)+
      ggplot2::geom_text(ggplot2::aes(label=mark),size=6/ggplot2::.pt,colour='black')+
      ggplot2::scale_fill_gradient2(low='#2166AC',mid='#FFFFFF',high='#B2182B',limits=c(-lim,lim),
        na.value='#D9D9D9',name='Standardized beta',breaks=c(-lim,0,lim),
        guide=ggplot2::guide_colourbar(display='rectangles',barwidth=grid::unit(25,'mm'),barheight=grid::unit(2.7,'mm'),title.position='top'))+
      ggplot2::scale_x_discrete(expand=ggplot2::expansion(add=0))+
      ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=0))+
      ggplot2::labs(x=NULL,y=NULL,tag=letters[i])
    fc_themed(g)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55,hjust=1,vjust=1),
      axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),legend.position='bottom')
  })
  capsule<-file.path(out,'audit/figure_source_capsule.R')
  writeLines(unlist(lapply(c('nf_theme','tf_export','fc_themed','f6_figures'),function(name)
    paste0(name,' <- ',paste(deparse(get(name,envir=.GlobalEnv),width.cutoff=120),collapse='\n')))),capsule,useBytes=TRUE)
  invisible(parse(capsule,encoding='UTF-8'))
  st<-system2(cfg$python,c(shQuote(file.path(nf_skill_root(),'scripts/validate_figure.py')),shQuote(capsule)),
    stdout=file.path(out,'audit/figure_static_preflight.txt'),stderr=file.path(out,'audit/figure_static_preflight_stderr.txt'))
  if(st!=0L)stop('Figure source preflight failed')
  manifest<-tf_export(list(plot=patchwork::wrap_plots(heat,nrow=1),width=183,height=105,panels=2,
    oriented_text=TRUE,exemptions=list()),'Figure6_Candidate_metabolic_family_LMM_20261105',out,cfg)
  v3_table(manifest,'Figure_export_manifest_EN',out)
  v3_table(data.frame(lower=-lim,midpoint=0,upper=lim,unit='Standardized inverse-normal beta',
    rule='One shared symmetric scale, no row or column rescaling'),'Fig6_heatmap_colour_scale_EN',out)
  # Candidate residual diagnostics: private R-rendered QA pages, never a manuscript figure.
  grDevices::cairo_pdf(file.path(out,'audit/candidate_residual_diagnostics.pdf'),width=10,height=7,family='Arial')
  on.exit(grDevices::dev.off(),add=TRUE)
  dg<-r$diagnostics
  for(field in names(inp$fields)) {
    z<-dg[dg$trait_id==field & dg$feature_id==ids[1],]
    par(mfrow=c(2,2),mar=c(4,4,3,1))
    for(mode in c('Primary','Steatosis-conditional')) {
      v<-z[z$mode==mode,];qqnorm(v$residual_z,main=paste(field,mode,'QQ'));qqline(v$residual_z,col='red')
      plot(v$fitted,v$residual_z,pch=16,cex=.25,col=grDevices::adjustcolor('grey30',alpha.f=.3),
        xlab='Fitted',ylab='Standardized residual',main=paste(field,mode));abline(h=0,col='red')
      lines(lowess(v$fitted,v$residual_z),col='blue',lwd=1)
    }
  }
}
f6_main <- function(args=commandArgs(TRUE),root='E:/YFS/workspace/FattyLiver_generation_analysis') {
  stopifnot(all(args%in%c('--run','--resume','--figures-only')))
  options(fatty_liver.autorun=FALSE);invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
  source(file.path(root,'R/family_candidate_followup.R'),encoding='UTF-8')
  source(file.path(root,'R/family_four_stage.R'),encoding='UTF-8');cfg<-f4_setup(root)
  cfg$output_module<-'51_Fig6_family_mixed_models';cfg$run_id<-'fig6_family_LMM_20261105_v1';cfg$run_id_explicit<-FALSE
  existing<-list.dirs(file.path(root,'outputs',cfg$output_module,'_running'),recursive=FALSE,full.names=TRUE)
  do_resume<-any(c('--resume','--figures-only')%in%args)
  if(do_resume && length(existing)) {
    out<-existing[which.max(file.info(existing)$mtime)];cfg$run_id<-basename(out)
    ctx<-list(cfg=cfg,module=dirname(dirname(out)),run_dir=out)
  } else {ctx<-fl_output_start_run(cfg);out<-ctx$run_dir}
  for(k in c('tables','reports','private','audit','figures'))dir.create(file.path(out,k),showWarnings=FALSE)
  if(do_resume && !file.exists(file.path(out,'private/input_population.rds'))) {
    current_base<-file.path(root,'outputs',cfg$output_module)
    if(file.exists(file.path(current_base,'LATEST_RUN.txt'))) {
      adopted<-fc_adopted(root,cfg$output_module)
      for(k in c('input_population.rds','analysis_result.rds'))
        stopifnot(file.copy(file.path(adopted,'private',k),file.path(out,'private',k)))
    } else if('--figures-only'%in%args)stop('No saved analysis available for figures-only')
  }
  writeLines(c('RUNNING',f6_protocol),file.path(out,'run_status.txt'))
  tryCatch({
    input<-file.path(out,'private/input_population.rds')
    inp<-if(file.exists(input))readRDS(input)else f6_prepare(root,out)
    ade_protect(inp$protected)
    result<-file.path(out,'private/analysis_result.rds')
    if('--figures-only'%in%args && !file.exists(result))stop('No complete saved analysis; use --resume to finish inference')
    r<-if(file.exists(result))readRDS(result)else f6_analyze(inp,out)
    stopifnot(identical(r$protocol,f6_protocol));saveRDS(r,result)
    counts<-f6_tables_reports(r,inp,out);f6_figures(r,inp,out,cfg)
    ade_protect(inp$protected);fl_write_csv(inp$protected,file.path(out,'audit/protected_source_fingerprints.csv'))
    writeLines(capture.output(sessionInfo()),file.path(out,'audit/sessionInfo.txt'))
    writeLines(c('COMPLETED',f6_protocol),file.path(out,'run_status.txt'))
    target<-fl_output_complete_run(ctx);writeLines(basename(target),file.path(ctx$module,'LATEST_RUN.txt'))
    cat('COMPLETED',target,'\n');print(counts);invisible(target)
  },error=function(e){writeLines(c('FAILED_RETAINED_FOR_RESUME',conditionMessage(e)),file.path(out,'run_status.txt'));stop(e)})
}
