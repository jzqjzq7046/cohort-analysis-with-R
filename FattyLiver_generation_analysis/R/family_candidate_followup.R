# Inert post-discovery characterization. No discovery refits, no strain inference.
fc_protocol <- 'family_candidate_characterization_20261104_v1'
fc_adopted <- function(root, module) {
  base <- file.path(root, 'outputs', module)
  id <- readLines(file.path(base, 'LATEST_RUN.txt'), 1, warn=FALSE)
  stopifnot(length(id)==1L, grepl('^[A-Za-z0-9_.-]+$', id))
  path <- file.path(base, 'current', id)
  stopifnot(readLines(file.path(path, 'run_status.txt'), 1)=='COMPLETED')
  path
}
fc_refs <- function(include_historical=FALSE) {
  tab<-data.frame(
  reference=c("Chen W, Li Y, Wang W, et al. Enhanced microbiota profiling in patients with quiescent Crohn's disease through comparison with paired healthy first-degree relatives. Cell Reports Medicine. 2024;5(7):101624. doi:10.1016/j.xcrm.2024.101624",
    'Asnicar F, Manghi P, Fackelmann G, et al. Gut micro-organisms associated with health, nutrition and dietary interventions. Nature. 2026;650:450-458. doi:10.1038/s41586-025-09854-7',
    'Ho J, Tumkaya T, Aryal S, Choi H, Claridge-Chang A. Moving beyond P values: data analysis with estimation graphics. Nature Methods. 2019;16:565-566. doi:10.1038/s41592-019-0470-3',
    'Leonardi I, Paramsothy S, Doron I, et al. Fungal trans-kingdom dynamics linked to responsiveness to fecal microbiota transplantation (FMT) therapy in ulcerative colitis. Cell Host & Microbe. 2020;27(5):823-829.e3. doi:10.1016/j.chom.2020.03.006'),
  URL=c('https://pmc.ncbi.nlm.nih.gov/articles/PMC11293350/','https://www.nature.com/articles/s41586-025-09854-7',
    'https://doi.org/10.1038/s41592-019-0470-3','https://pmc.ncbi.nlm.nih.gov/articles/PMC8647676/'),
  adopted_content=c('Figures 2-4: differential-taxon profiles, contrast across controls and correlation with host markers',
    'Partial Spearman using relative abundance including zeros; age/sex/BMI adjustment; phenotype-association heatmaps',
    'Estimation graphics: raw data plus bootstrap effect distribution, mean contrast and BCa 95% confidence interval',
    'Extended Data 2/4: microbial relative-abundance estimation graphics with mean difference and bootstrap confidence interval'),
  adaptation=c('Steatosis/SGBs/unequal primary-family sets replace CD/genera/sibling pairs; descriptive state means are not independent paired individuals',
    'R rank-residual equivalent estimator; one adult per primary family; nine available traits; full-background BH and a separate steatosis-conditional sensitivity',
    'Transpose and compact all ten SGBs; equal-weight primary-family disease-state mean differences; 5000 family bootstraps implemented in R boot, not a DABEST software run',
    'Graphical approach only: family-paired SGB contrasts replace unpaired fungal-group comparisons; no fungal, FMT or intervention inference'),
  accessed='2026-11-04')
  if(include_historical)tab else tab[1:2,,drop=FALSE]
}
fc_family_data <- function(parent,ids,out) {
  discovery<-readRDS(file.path(parent,'private/analysis_result.rds'))
  stopifnot(identical(discovery$protocol,f4_protocol))
  input<-discovery$protected$path[1];fingerprint<-unname(tools::md5sum(input))
  stopifnot(identical(fingerprint,discovery$protected$md5[1]))
  source<-readRDS(input);m<-readRDS(file.path(parent,'private/model_populations.rds'))$B2
  ix<-match(m$sample_id,colnames(source$relative))
  stopifnot(!anyNA(ix),!anyDuplicated(ix),all(ids%in%rownames(source$relative)))
  A<-source$relative[ids,ix,drop=FALSE]
  stopifnot(all(is.finite(A)),all(A>=0),identical(colnames(A),as.character(m$sample_id)))
  fm<-ov_bind(lapply(ids,function(id){
    d<-data.frame(family_id=as.character(m$family_id),state=m$fatty_liver,abundance=100*as.numeric(A[id,]))
    a<-aggregate(abundance~family_id+state,d,mean)
    stopifnot(all(table(a$family_id)==2L),setequal(a$family_id,as.character(m$family_id)))
    a$feature_id<-id;a
  }))
  saveRDS(list(metadata=m,family_means=fm),file.path(out,'private/family_mean_plot_data.rds'))
  summary<-ov_bind(lapply(split(fm,list(fm$feature_id,fm$state),drop=TRUE),function(z)
    data.frame(feature_id=z$feature_id[1],state=z$state[1],families=nrow(z),mean=mean(z$abundance),
      median=median(z$abundance),Q1=unname(quantile(z$abundance,.25)),Q3=unname(quantile(z$abundance,.75)))))
  v3_table(summary,'Family_state_descriptive_summary_EN',out)
  fl_write_csv(data.frame(source_path=input,md5=fingerprint,participants=nrow(m),families=length(unique(m$family_id)),
    missing_sample_matches=sum(is.na(ix)),source_relative_abundances_preserved=TRUE),file.path(out,'audit/exact_family_display_input.csv'))
  stopifnot(identical(unname(tools::md5sum(input)),fingerprint))
  list(means=fm,summary=summary)
}
fc_read_sources <- function(root, out) {
  source <- tf_source(root); x <- tf_load(source)
  parent <- fc_adopted(root, '45_Family_de_novo_discovery')
  old <- fc_adopted(root, '25_Taxa_robustness')
  files <- c(source$files, file.path(parent, 'tables/All_tested_features_EN.csv'),
    file.path(parent, 'private/model_populations.rds'),
    file.path(old, 'tables', paste0('Disease_',c('B1','B2'),'_all_SGBs_EN.csv')))
  stopifnot(all(file.exists(files)))
  protected <- data.frame(path=files, md5=unname(tools::md5sum(files)))
  e <- tf_read(files[5]); primary <- e[e$model=='B2', ]
  stopifnot(!anyDuplicated(primary$feature_id), all(is.finite(primary$P_BH)),
    max(abs(primary$P_BH-p.adjust(primary$P_raw,'BH')))<1e-12)
  candidates <- primary[primary$P_BH < .10, ]
  candidates <- candidates[order(candidates$P_BH,candidates$P_raw,candidates$feature_id), ]
  stopifnot(nrow(candidates)>0, all(candidates$feature_id %in% rownames(x$relative)))
  candidates$discovery_tier <- ifelse(candidates$P_BH<.05,'Primary: BH q < 0.05','Exploratory: 0.05 <= BH q < 0.10')
  candidates$display_order <- seq_len(nrow(candidates))
  candidates$SGB_id <- sub('^t__','',candidates$feature_id)
  candidates$plot_label <- candidates$feature_label
  # Preserve explicit species annotations, retain SGB IDs for unresolved/conflicting names.
  named <- !is.na(candidates$scientific_name) & nzchar(candidates$scientific_name)
  candidates$plot_label[named] <- paste0(candidates$scientific_name[named],' (',candidates$SGB_id[named],')')
  m <- readRDS(file.path(parent,'private/model_populations.rds'))$B2
  stopifnot(!anyDuplicated(m$sample_id), !anyDuplicated(m$participant_id),
    all(tapply(m$fatty_liver,m$family_id,function(z)length(unique(z)))==2))
  evidence <- list()
  for(k in c('B1','B2')) {
    z <- tf_read(file.path(old,'tables',paste0('Disease_',k,'_all_SGBs_EN.csv')))
    z <- z[match(candidates$feature_id,z$feature_id), ]
    stopifnot(!anyNA(z$feature_id))
    evidence[[paste0('old_',k)]] <- data.frame(feature_id=z$feature_id,
      analysis=if(k=='B1')'Overall\nModel 3'else'Adults\nModel 4',
      framework='MaAsLin2: family random intercept', q=z$P_BH, effect=z$beta_log2,
      unit='log2 zero-replaced relative abundance', direction=sign(z$beta_log2))
  }
  for(k in c('B1','Adult_B1_same_sample','B2')) {
    z <- e[e$model==k, ]; z <- z[match(candidates$feature_id,z$feature_id), ]
    stopifnot(!anyNA(z$feature_id))
    evidence[[paste0('new_',k)]] <- data.frame(feature_id=z$feature_id,
      analysis=switch(k,B1='Within family\nall-age Model 3',Adult_B1_same_sample='Within family\nadult Model 3',B2='Within family\nadult Model 4'),
      framework='LDM: family fixed effects and within-family permutations', q=z$P_BH,
      effect=z$beta_percentage_points, unit='relative-abundance percentage points', direction=sign(z$beta_relative))
  }
  evidence <- ov_bind(evidence)
  family_data<-fc_family_data(parent,candidates$feature_id,out)
  list(x=x,candidates=candidates,evidence=evidence,family_summary=family_data$summary,protected=protected,
    primary_N=nrow(m),primary_cases=sum(m$fatty_liver),primary_families=nlevels(m$family_id))
}
fc_partial_spearman <- function(x,y,controls) {
  ok <- is.finite(x)&is.finite(y)&apply(controls,1,function(z)all(is.finite(z)))
  x <- x[ok]; y <- y[ok]; controls <- controls[ok,,drop=FALSE]; n <- length(x)
  empty <- function(status) data.frame(N=n,df=NA_real_,rho=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P_raw=NA_real_,status=status)
  if(n<20 || length(unique(x))<2 || length(unique(y))<2)return(empty('NOT_ESTIMABLE: small N or constant measurement'))
  R <- vapply(seq_len(ncol(controls)),function(j)rank(controls[,j],ties.method='average'),numeric(n))
  Z <- cbind(Intercept=1,R); qrZ <- qr(Z,tol=1e-9); df <- n-qrZ$rank-1
  if(df<=3)return(empty('NOT_ESTIMABLE: residual degrees of freedom'))
  rx <- qr.resid(qrZ,rank(x,ties.method='average')); ry <- qr.resid(qrZ,rank(y,ties.method='average'))
  if(sd(rx)<1e-12 || sd(ry)<1e-12)return(empty('NOT_ESTIMABLE: constant residual ranks'))
  rho <- max(-1,min(1,cor(rx,ry))); statistic <- rho*sqrt(df/max(1-rho^2,.Machine$double.eps))
  ci <- tanh(atanh(max(-1+1e-15,min(1-1e-15,rho)))+c(-1,1)*qnorm(.975)/sqrt(n-qrZ$rank-2))
  data.frame(N=n,df=df,rho=rho,CI_low=ci[1],CI_high=ci[2],P_raw=2*pt(-abs(statistic),df),status='ESTIMATED')
}
fc_host_map <- function(x,candidates,out) {
  fields <- c(BMI18='BMI',WAIST18='Waist',syst18='SBP',GLUK18='Glucose',HbA1c18='HbA1c',
    TRIGLY18='Triglycerides',HDLKOL18='HDL-C',ALAT18='ALT',ASAT18='AST')
  m <- x$metadata; stopifnot(all(names(fields)%in%names(m)))
  eligible <- !is.na(m$age)&m$age>=18&v3_cc(m,c('family_id','age','sex','BMI18'))&m$BMI18>0
  mm <- droplevels(m[eligible,]); mm <- mm[order(mm$family_id,mm$sample_id), ]
  set.seed(20261104L)
  ix <- sort(vapply(split(seq_len(nrow(mm)),mm$family_id),function(ii)ii[sample.int(length(ii),1L)],integer(1)))
  mm <- droplevels(mm[ix,]);stopifnot(!anyDuplicated(mm$family_id),!anyDuplicated(mm$sample_id))
  mm$sex_indicator <- as.integer(factor(mm$sex))-1L
  A <- x$relative[,match(mm$sample_id,colnames(x$relative)),drop=FALSE]
  prevalence <- rowMeans(A>0); ids <- rownames(A)[prevalence>.20]
  filter <- data.frame(feature_id=rownames(A),prevalence=prevalence,eligible=prevalence>.20,
    candidate=rownames(A)%in%candidates$feature_id,rule='Detected in >20% of outcome-blind one-adult-per-family population')
  flow <- data.frame(stage=c('Adopted microbiome cohort','Age/sex/BMI-complete adults','One outcome-blind adult per family'),
    N=c(nrow(m),sum(eligible),nrow(mm)))
  saveRDS(list(metadata=mm,seed=20261104L,selection='Before phenotype availability or abundance filtering'),file.path(out,'private/host_map_population.rds'))
  rows <- list(); r <- 0L
  for(mode in c('Source-style','Steatosis-conditional'))for(field in names(fields)) {
    controls <- c('age','sex_indicator',if(field!='BMI18')'BMI18',if(mode=='Steatosis-conditional')'fatty_liver')
    y <- mm[[field]]; lower <- if(field%in%c('ALAT18','ASAT18'))0 else .Machine$double.eps
    y[!is.finite(y)|y<lower] <- NA_real_
    Z <- as.matrix(mm[,controls,drop=FALSE]); storage.mode(Z)<-'double'
    for(id in ids) {
      r <- r+1L
      rows[[r]] <- cbind(feature_id=id,phenotype=field,phenotype_label=unname(fields[field]),
        mode=mode,adjustment=paste(controls,collapse=' + '),fc_partial_spearman(as.numeric(A[id,]),y,Z))
    }
  }
  tests <- ov_bind(rows);tests$P_BH<-NA_real_
  for(mode in unique(tests$mode)) {
    j <- which(tests$mode==mode);good<-is.finite(tests$P_raw[j])&tests$status[j]=='ESTIMATED'
    p <- rep(1,length(j));p[good]<-tests$P_raw[j][good]
    tests$P_BH[j[good]]<-p.adjust(p,'BH')[good]
  }
  tests$FDR_family <- 'All eligible background SGBs x all nine prespecified traits within mode'
  tests$candidate <- tests$feature_id%in%candidates$feature_id
  display <- expand.grid(feature_id=candidates$feature_id,phenotype=names(fields),mode=unique(tests$mode),stringsAsFactors=FALSE)
  display <- merge(display,tests,by=c('feature_id','phenotype','mode'),all.x=TRUE,sort=FALSE)
  display$phenotype_label<-unname(fields[display$phenotype]);display$status[is.na(display$status)]<-'NOT_TESTED: prevalence <=20%'
  display$strict <- is.finite(display$P_BH)&display$P_BH<.05
  display$exploratory <- is.finite(display$P_BH)&display$P_BH>=.05&display$P_BH<.10
  summary <- ov_bind(lapply(split(tests,list(tests$mode,tests$phenotype),drop=TRUE),function(z)
    data.frame(mode=z$mode[1],phenotype=z$phenotype[1],tested=nrow(z),N_min=min(z$N),N_max=max(z$N),
      all_SGB_q005=sum(z$P_BH<.05,na.rm=TRUE),candidate_q005=sum(z$candidate&z$P_BH<.05,na.rm=TRUE),
      candidate_q010=sum(z$candidate&z$P_BH<.10,na.rm=TRUE))))
  list(tests=tests,display=display,summary=summary,filter=filter,flow=flow,fields=fields,
    N=nrow(mm),families=length(unique(mm$family_id)),source_eligible=sum(eligible))
}
fc_validate <- function(r,out) {
  c <- r$candidates; h <- r$host
  stopifnot(!anyDuplicated(c$feature_id),all(c$P_BH<.1),nrow(r$evidence)==5*nrow(c),
    nrow(h$display)==18*nrow(c),h$N==h$families,all(h$tests$rho[is.finite(h$tests$rho)]>=-1),
    all(h$tests$rho[is.finite(h$tests$rho)]<=1))
  for(mode in unique(h$tests$mode)) {
    z<-h$tests[h$tests$mode==mode,];good<-is.finite(z$P_raw)&z$status=='ESTIMATED'
    p<-rep(1,nrow(z));p[good]<-z$P_raw[good]
    stopifnot(max(abs(z$P_BH[good]-p.adjust(p,'BH')[good]))<1e-12)
  }
  # Algebraic consistency: no-controls partial Spearman equals ordinary Spearman.
  a<-c(0,0,seq_len(28));b<-c(5:1,1:25)
  z<-fc_partial_spearman(a,b,matrix(1,length(a),1))
  stopifnot(abs(z$rho-cor(a,b,method='spearman'))<1e-12)
  ade_protect(r$protected)
  v3_table(data.frame(check=c('Frozen complete-discovery q values','Full background BH reproduced','One adult per family',
    'All candidates retained in all displays','Partial-Spearman estimator equivalence','Protected inputs unchanged'),passed=TRUE),
    'Statistical_validation_EN',out)
}
fc_tables <- function(r,out) {
  for(n in c('candidates','evidence','family_summary'))v3_table(r[[n]],switch(n,candidates='All_selected_candidates_EN',
    evidence='Saved_model_evidence_EN',family_summary='Family_state_descriptive_summary_EN'),out)
  h<-r$host
  for(n in c('tests','display','summary','filter','flow'))v3_table(h[[n]],paste0('Host_',n,'_EN'),out)
  v3_table(fc_refs(),'References_and_adaptations_EN',out)
  overlap<-ov_bind(lapply(c(.05,.10),function(q) {
    ids<-r$candidates$feature_id[r$candidates$P_BH<q]
    ov_bind(lapply(c('Overall\nModel 3','Adults\nModel 4'),function(a){z<-r$evidence[r$evidence$analysis==a & r$evidence$feature_id%in%ids, ]
      current<-r$candidates$beta_relative[match(z$feature_id,r$candidates$feature_id)]
      data.frame(within_family_threshold=q,saved_comparison=gsub('\n',' ',a),candidates=length(ids),
        previous_q005=sum(z$q<.05),same_direction=sum(sign(z$effect)==sign(current)))}))
  }))
  v3_table(overlap,'Candidate_overlap_and_direction_EN',out)
}
fc_themed <- function(p) p+tf_plot_theme()+ggplot2::theme(plot.title=ggplot2::element_blank(),
  plot.subtitle=ggplot2::element_text(size=5.7),legend.title=ggplot2::element_text(size=5.7),
  plot.tag=ggplot2::element_text(size=7,face='bold'))
fc_family_estimation <- function(r,out) {
  input<-file.path(out,'private/family_mean_plot_data.rds')
  signature<-digest::digest(list(input=unname(tools::md5sum(input)),method='equal-family-mean-BCa-5000',seed=20261104L),algo='sha256')
  checkpoint<-file.path(out,'private/family_estimation_bootstrap.rds')
  if(file.exists(checkpoint)) {
    old<-readRDS(checkpoint)
    if(identical(old$signature,signature))return(old$value)
  }
  stopifnot(requireNamespace('boot',quietly=TRUE))
  fm<-readRDS(input)$family_means;raw<-list();summary<-list();distributions<-list()
  for(i in seq_len(nrow(r$candidates))) {
    id<-r$candidates$feature_id[i];z<-fm[fm$feature_id==id,]
    a<-z[z$state==0,];b<-z[z$state==1,];b<-b[match(a$family_id,b$family_id),]
    stopifnot(!anyNA(b$family_id),identical(a$family_id,b$family_id),nrow(a)==r$primary_families)
    delta<-b$abundance-a$abundance;set.seed(20261104L+i)
    bt<-boot::boot(delta,function(d,j)mean(d[j]),R=5000)
    ci<-boot::boot.ci(bt,type='bca',L=delta-mean(delta))$bca[1,4:5]
    stopifnot(all(is.finite(ci)),ci[1]<=ci[2],abs(bt$t0-mean(delta))<1e-12)
    raw[[i]]<-data.frame(feature_id=id,family_id=a$family_id,difference_pp=delta)
    summary[[i]]<-data.frame(feature_id=id,families=length(delta),mean_difference_pp=mean(delta),
      CI_low=ci[1],CI_high=ci[2],bootstrap_replicates=5000,seed=20261104L+i,
      saved_adjusted_discovery_q=r$candidates$P_BH[i],
      interpretation='Descriptive equal-family contrast; not the adjusted LDM coefficient or a new discovery test')
    distributions[[i]]<-data.frame(feature_id=id,bootstrap_mean_pp=as.numeric(bt$t))
  }
  value<-list(raw=ov_bind(raw),summary=ov_bind(summary),bootstrap=ov_bind(distributions))
  saveRDS(list(signature=signature,value=value),checkpoint);value
}
fc_figures <- function(r,cfg,out,redraw=c('evidence','family','host')) {
  dir.create(file.path(out,'figures'),showWarnings=FALSE)
  capsule<-file.path(out,'audit/figure_source_capsule.R')
  writeLines(unlist(lapply(c('nf_theme','tf_export','fc_themed','fc_family_estimation','fc_figures'),function(n)
    paste0(n,' <- ',paste(deparse(get(n,envir=.GlobalEnv),width.cutoff=120),collapse='\n')))),capsule)
  invisible(parse(capsule,encoding='UTF-8'))
  system2(cfg$python,c(shQuote(file.path(nf_skill_root(),'scripts/validate_figure.py')),shQuote(capsule)),
    stdout=file.path(out,'audit/figure_static_preflight.txt'),stderr=file.path(out,'audit/figure_static_preflight_stderr.txt'))
  labels<-setNames(vapply(r$candidates$plot_label,function(a)paste(strwrap(a,34),collapse='\n'),character(1)),r$candidates$feature_id)
  ids<-r$candidates$feature_id; e<-r$evidence
  exports<-list()
  if('evidence'%in%redraw) {
  e$row<-factor(e$feature_id,levels=rev(ids),labels=rev(labels))
  e$mark<-ifelse(e$q<.05,'*',ifelse(e$q<.10,'+',''))
  e$framework_display<-ifelse(grepl('MaAsLin2',e$framework),'MaAsLin2','Within-family LDM')
  e$display_value<-ifelse(e$framework_display=='MaAsLin2',2^e$effect,e$effect)
  e$unit_display<-ifelse(e$framework_display=='MaAsLin2','Abundance ratio','Difference (pp)')
  e$cell_label<-vapply(seq_len(nrow(e)),function(i)paste0(format(signif(e$display_value[i],2),trim=TRUE,scientific=FALSE),
    if(nzchar(e$mark[i]))paste0('\n',e$mark[i])else''),character(1))
  plots<-lapply(seq_along(unique(e$framework_display)),function(i){
    method<-unique(e$framework_display)[i];d<-e[e$framework_display==method,]
    d$column<-factor(d$analysis,levels=unique(d$analysis),labels=if(i==1)c('All ages\nModel 3','Adults\nModel 4')else c('All ages\nModel 3','Adults\nModel 3','Adults\nModel 4'))
    lim<-if(i==1)max(1,ceiling(max(abs(d$effect))))else max(.25,ceiling(max(abs(d$effect))/.25)*.25)
    ticks<-c(-lim,0,lim)
    tick_labels<-if(i==1)format(signif(2^ticks,2),trim=TRUE,scientific=FALSE)else format(signif(ticks,2),trim=TRUE)
    g<-ggplot2::ggplot(d,ggplot2::aes(column,row,fill=effect))+
      ggplot2::geom_tile(colour='white',linewidth=.45)+
      ggplot2::geom_text(ggplot2::aes(label=cell_label),size=5.7/ggplot2::.pt,lineheight=.90,colour='black')+
      ggplot2::scale_fill_gradient2(low='#79CAFB',mid='white',high='#F7A24F',limits=c(-lim,lim),
        name=if(i==1)'Abundance ratio'else'Difference (pp)',breaks=ticks,labels=tick_labels,
        guide=ggplot2::guide_colourbar(display='rectangles',barwidth=grid::unit(28,'mm'),barheight=grid::unit(2.7,'mm'),title.position='top'))+
      ggplot2::scale_x_discrete(expand=ggplot2::expansion(add=0))+
      ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=0))+
      ggplot2::labs(x=method,y=NULL,tag=letters[i])
    fc_themed(g)+ggplot2::theme(axis.ticks=ggplot2::element_blank(),axis.line=ggplot2::element_blank(),
      legend.position='bottom',axis.title.x=ggplot2::element_text(size=6.2,margin=ggplot2::margin(t=5,unit='pt')))
  })
  exports[[length(exports)+1L]]<-tf_export(list(plot=patchwork::wrap_plots(plots,nrow=1),width=183,height=135,panels=2,exemptions=list()),
    'Figure50_Candidate_model_evidence_v3',out,cfg)
  v3_table(e,'Model_evidence_heatmap_display_EN',out)
  }
  if('family'%in%redraw) {
    # Never inherit the earlier incomplete stage-24 plotting subset on resume.
    fc_family_data(fc_adopted(cfg$project_root,'45_Family_de_novo_discovery'),ids,out)
    fm<-readRDS(file.path(out,'private/family_mean_plot_data.rds'))$family_means
    plots<-lapply(seq_along(ids),function(i){
      id<-ids[i];d<-fm[fm$feature_id==id,];d$state<-factor(d$state,levels=c(0,1),labels=c('No steatosis','Steatosis'))
      boxes<-ov_bind(lapply(split(d,d$state,drop=TRUE),function(z){
        qs<-unname(quantile(z$abundance,c(.25,.5,.75)));iqr<-qs[3]-qs[1]
        inside<-z$abundance[z$abundance>=qs[1]-1.5*iqr&z$abundance<=qs[3]+1.5*iqr]
        data.frame(state=z$state[1],ymin=min(inside),lower=qs[1],middle=qs[2],upper=qs[3],ymax=max(inside))
      }))
      g<-ggplot2::ggplot(d,ggplot2::aes(state,abundance))+
        ggplot2::geom_line(ggplot2::aes(group=family_id),colour='#999999',alpha=.16,linewidth=.16)+
        ggplot2::geom_boxplot(data=boxes,ggplot2::aes(x=state,ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,fill=state),
          inherit.aes=FALSE,stat='identity',colour='black',linewidth=.3,width=.50)+
        ggplot2::geom_point(position=ggplot2::position_jitter(width=.06,height=0,seed=50L+i),colour='#777777',alpha=.25,size=.24)+
        ggplot2::scale_fill_manual(values=c('No steatosis'='#79CAFB','Steatosis'='#F7A24F'))+
        ggplot2::scale_y_continuous(trans=scales::pseudo_log_trans(sigma=.001,base=10),breaks=c(0,.01,.1,1,10,100),
          labels=function(a)vapply(a,function(v)if(is.na(v))NA_character_ else if(v==0)'0'else format(v,scientific=FALSE,trim=TRUE),character(1)))+
        ggplot2::labs(x=NULL,y='Relative abundance (%)',tag=letters[i],
          subtitle=paste0(labels[id],' | BH q = ',formatC(r$candidates$P_BH[i],digits=3,format='f')))
      fc_themed(g)+ggplot2::theme(legend.position='none')
    })
    exports[[length(exports)+1L]]<-tf_export(list(plot=patchwork::wrap_plots(plots,ncol=2),width=183,height=48*ceiling(length(ids)/2),
      panels=length(ids),uniform_y_axes=TRUE,exemptions=list()),'Figure50_Family_abundance_profiles_restored_v3',out,cfg)
  }
  if('host'%in%redraw) {
  d<-r$host$display; d$row<-factor(d$feature_id,levels=rev(ids),labels=rev(labels))
  d$trait<-factor(d$phenotype,levels=names(r$host$fields),labels=unname(r$host$fields))
  d$mark<-ifelse(d$strict,'*',ifelse(d$exploratory,'+',ifelse(!is.finite(d$rho),'x','')))
  lim<-max(.05,ceiling(max(abs(d$rho),na.rm=TRUE)/.05)*.05)
  heat<-lapply(seq_along(unique(d$mode)),function(i){mode<-unique(d$mode)[i];z<-d[d$mode==mode, ]
    g<-ggplot2::ggplot(z,ggplot2::aes(trait,row,fill=rho))+
      ggplot2::geom_tile(colour='white',linewidth=.35)+
      ggplot2::geom_text(ggplot2::aes(label=mark),size=6/ggplot2::.pt,colour='black')+
      ggplot2::scale_fill_gradient2(low='#2166AC',mid='#FFFFFF',high='#B2182B',limits=c(-lim,lim),
        na.value='#D9D9D9',name='Partial rho',breaks=c(-lim,0,lim),
        guide=ggplot2::guide_colourbar(display='rectangles',barwidth=grid::unit(25,'mm'),barheight=grid::unit(2.7,'mm'),title.position='top'))+
      ggplot2::scale_x_discrete(position='bottom',expand=ggplot2::expansion(add=0))+
      ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=0))+
      ggplot2::labs(x=NULL,y=NULL,tag=letters[i])
    fc_themed(g)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55,hjust=1,vjust=1),
      axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),legend.position='bottom')
  })
  exports[[length(exports)+1L]]<-tf_export(list(plot=patchwork::wrap_plots(heat,nrow=1),width=183,height=105,panels=2,
    oriented_text=TRUE,exemptions=list()),'Figure50_Candidate_host_phenotypes_v2',out,cfg)
  v3_table(data.frame(lower=-lim,midpoint=0,upper=lim,rule='One shared symmetric rho scale; no row/column scaling'),
    'Host_heatmap_colour_scale_EN',out)
  }
  oldmanifest<-file.path(out,'tables/Figure_export_manifest_EN.csv');newmanifest<-ov_bind(exports)
  if(file.exists(oldmanifest)) {
    old<-tf_read(oldmanifest);old<-old[!old$figure%in%newmanifest$figure,]
    newmanifest<-ov_bind(list(old,newmanifest))
  }
  v3_table(newmanifest,'Figure_export_manifest_EN',out)
}
fc_reports <- function(r,out) {
  h<-r$host;d<-h$display
  counts<-ov_bind(lapply(unique(d$mode),function(mode){z<-d[d$mode==mode,];data.frame(mode=mode,
    displayed_tests=nrow(z),estimable=sum(is.finite(z$rho)),q005=sum(z$strict),q010=sum(z$strict|z$exploratory))}))
  v3_table(counts,'Candidate_host_association_counts_EN',out)
  supported<-d[is.finite(d$P_BH)&d$P_BH<.10,]
  supported<-supported[order(supported$mode,supported$P_BH,supported$feature_id,supported$phenotype),]
  v3_table(supported,'Candidate_host_exploratory_supported_links_EN',out)
  writeLines(c('# Post-discovery candidate characterization',paste('Protocol:',fc_protocol),'',
    paste('All',nrow(r$candidates),'adult Model 4 candidates with saved conventional BH q<0.10 are displayed;',
      sum(r$candidates$P_BH<.05),'pass q<0.05. The original complete-discovery BH correction is retained.'),
    paste('Discovery population:',r$primary_N,'participants,',r$primary_cases,'cases,',r$primary_families,'primary families.'),
    'The five-model evidence heatmap follows the coefficient-plus-significance display of Chen et al. Figure 6A. Left: MaAsLin2 LOG effects displayed as 2^coefficient abundance ratios, on a log2 colour scale centered at ratio=1. Right: LDM adjusted frequency contrasts in percentage points, on a separate symmetric linear scale centered at zero. Units and scales are not pooled. Cell numbers are effect values; * q<0.05, + 0.05<=q<0.10, using saved full-feature BH. All source models account for family structure. Same-cohort comparisons are not independent validation.',
    'Family abundance profiles restore connected disease-state means and boxplots. Grey lines connect the affected and unaffected mean relative abundance in each of 289 primary families; boxes and grey points summarize/represent family-state means. The exact discovery matrix includes all 731 participants. Original individual-level adjusted LDM P and BH q remain unchanged; no new paired test or bootstrap is run.',
    '',paste('Host map:',h$N,'adults, one per family, seed 20261104. Selection is before trait filtering and independent of abundance/outcome values.'),
    paste('Background SGBs tested:',sum(h$filter$eligible),'; nine traits. BH covers all background SGB-by-trait tests within each mode, not only selected disease candidates.'),
    'Partial Spearman uses relative abundance including zeros and average ranks for ties; age/sex/BMI are controlled. BMI itself adjusts age/sex only. A distinct sensitivity additionally conditions on ultrasound steatosis. This source-style reduced adjustment is for host characterization, not a replacement for discovery Model 4.',
    'No participant-level imputation. Observed complete cases are used per trait. Non-estimable planned tests count as P=1 in BH and retain NA results. Source populations overlap; host characterization is exploratory, not an independent validation or mechanistic experiment.',
    '',capture.output(print(counts,row.names=FALSE)),
    '', 'SGBs are species-level genome bins, not strain types. No heritability, directional transmission, mediation, diagnostic performance or strain-specific functional attribution is inferred.',
    'The q<0.10 disease threshold was chosen after viewing q<0.05 results and is explicitly exploratory. The saved discovery Monte Carlo P resolution is 1/10001 at the floor; no extra permutation was run.',
    '',apply(fc_refs(),1,function(z)paste0('- ',z['reference'],' (',z['adopted_content'],'; adaptation: ',z['adaptation'],')'))),
    file.path(out,'reports/Methods_results_and_limits_EN.md'),useBytes=TRUE)
  writeLines(c('Figure 50 | Characterization of within-family steatosis-associated SGBs.',
    'Model evidence heatmaps: all ten adult Model 4 candidates with saved BH q<0.10. (a) Previous MaAsLin2 mixed-model abundance ratios (2^coefficient), with colour on the log2-ratio scale centered at ratio=1. Ratios refer to zero-replaced relative abundance. (b) Within-family LDM adjusted frequency differences in percentage points, with a separate symmetric linear colour scale centered at zero. Cell numbers show the effects; blue/orange indicate lower/higher abundance. * BH q<0.05; + 0.05<=BH q<0.10; all q values retain their complete-feature test family. Do not compare effect magnitudes across panels or interpret the same-cohort comparison as independent validation. The display structure follows Chen et al. Cell Reports Medicine 2024, Figure 6A; disease, SGB level, comparison columns, units and BH markers are explicit adaptations.',
    'Family abundance profiles: all ten candidates, 731 participants and 289 primary families. Grey lines connect no-steatosis and steatosis means within the same family; grey points represent those means. Boxes show median/IQR across family means, whiskers extend to observations within 1.5 raw-scale IQR. Blue/orange denote no steatosis/steatosis. The pseudo-log percentage display (sigma=0.001 percentage points) preserves zeros; ticks remain fixed percentage values. Subtitles show the saved individual-level adult Model 4 LDM BH q, not a new test of plotted means. The plot restores the original connected-box design using the corrected complete discovery matrix.',
    'Host phenotype map: all candidates against nine prespecified traits. (a) Partial Spearman adjusted for age, sex and BMI (BMI adjusted for age/sex only). (b) Same adjustment plus steatosis status. One outcome-blind selected adult per primary family; relative abundances include zeros; per-trait complete measurements. Colour shows partial rho on the same symmetric scale; * q<0.05, + 0.05<=q<0.10; x, not tested/not estimable. Conventional BH covers all eligible background SGBs times nine traits separately within each panel, not just the displayed candidates. Rows retain discovery ranking; columns follow clinical domains. These exploratory cross-sectional associations do not establish mechanisms.'),
    file.path(out,'reports/Figure_legends_EN.txt'),useBytes=TRUE)
  first<-counts[counts$mode=='Source-style',];second<-counts[counts$mode=='Steatosis-conditional',]
  writeLines(c('# 家系差异菌的文献式后续分析',
    paste0('纳入成人Model 4的全部',nrow(r$candidates),'个BH q<0.10候选SGB，其中',sum(r$candidates$P_BH<.05),'个属于q<0.05主要发现。原发现模型与BH校正保持不变。'),
    '1. 模型一致性：将既往总体/成人MaAsLin2与本轮三套家系内LDM方向和q分类并列展示。两种方法的系数单位不同，不直接比较系数大小；同一队列证据不称为独立验证。',
    paste0('2. 家系内丰度图恢复原来的连线箱线图；',r$primary_N,'人、',r$primary_families,'个家系。每条灰线连接同一家系患病与未患病成员的平均相对丰度；蓝橙箱体、黑色线条、灰色散点。使用修正后的完整发现矩阵，不恢复旧的缺失样本问题。'),
    paste0('3. 宿主指标：每家随机选一名成人，共',h$N,'人。',sum(h$filter$eligible),'个检出率>20%的背景SGB×9项指标，共',9*sum(h$filter$eligible),'次检验/调整方式；BH不是仅在90个候选组合中计算。'),
    paste0('候选菌中有',sum(h$filter$candidate&h$filter$eligible),'个满足相关分析的文献检出率要求，其余保留空缺标记，不强行检验。'),
    paste0('年龄、性别及BMI调整后：候选菌—表型关联',first$q005,'组q<0.05，',first$q010,'组q<0.10。'),
    paste0('进一步调整脂肪肝状态后：',second$q005,'组q<0.05，',second$q010,'组q<0.10。'),
    '上述是关联组合数，不是菌的数量。BMI自身作为研究指标时不再校正BMI；九项指标为BMI、腰围、收缩压、空腹血糖、HbA1c、甘油三酯、HDL-C、ALT和AST。',
    '该图说明候选菌与宿主表型的横断面关系，不是致病机制、中介作用或独立人群验证。加入脂肪肝状态后的减弱也不能单凭此判断因果方向。',
    '方法依据：Chen等，Cell Reports Medicine，2024；Asnicar等，Nature，2026。具体采用内容和队列适配见References_and_adaptations_EN。'),
    file.path(out,'reports/Results_summary_ZH.md'),useBytes=TRUE)
}
fc_main <- function(args=character(),root=Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')) {
  stopifnot(all(args%in%c('--run','--check','--resume')))
  if(getRversion()<'4.5') {
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts/50_family_candidates_followup.R')),args))
    if(status!=0)stop('Project-R candidate follow-up failed');return(invisible(NULL))
  }
  sys.source(file.path(root,'R/family_four_stage.R'),envir=.GlobalEnv);cfg<-f4_setup(root)
  cfg$output_module<-'50_Family_candidate_characterization';cfg$run_id<-'family_candidates_20261104_v1';cfg$run_id_explicit<-FALSE
  if('--check'%in%args){fc_adopted(root,'45_Family_de_novo_discovery');fc_adopted(root,'25_Taxa_robustness');cat('Candidate inputs and R dependencies available\n');return(invisible(NULL))}
  ctx<-fl_output_start_run(cfg);out<-ctx$run_dir
  for(d in c('tables','reports','audit','private','figures'))dir.create(file.path(out,d),showWarnings=FALSE)
  tryCatch({
    if('--resume'%in%args) {
      paths<-list.files(file.path(ctx$module,'_running'),pattern='^analysis_result[.]rds$',recursive=TRUE,full.names=TRUE)
      if(!length(paths))paths<-list.files(file.path(ctx$module,'current'),pattern='^analysis_result[.]rds$',recursive=TRUE,full.names=TRUE)
      if(!length(paths))paths<-list.files(file.path(ctx$module,'history'),pattern='^analysis_result[.]rds$',recursive=TRUE,full.names=TRUE)
      stopifnot(length(paths)>0L);previous<-paths[which.max(file.info(paths)$mtime)];r<-readRDS(previous)
      stopifnot(identical(r$protocol,fc_protocol));oldout<-dirname(dirname(previous))
      for(f in c('family_mean_plot_data.rds','host_map_population.rds'))stopifnot(file.copy(file.path(oldout,'private',f),file.path(out,'private',f),overwrite=TRUE))
    } else {
      s<-fc_read_sources(root,out);h<-fc_host_map(s$x,s$candidates,out)
      s$x<-NULL;r<-c(s,list(host=h,protocol=fc_protocol))
    }
    saveRDS(r,file.path(out,'private/analysis_result.rds'))
    fc_validate(r,out);fc_tables(r,out);fc_reports(r,out);fc_figures(r,cfg,out)
    ade_protect(r$protected);fl_write_csv(r$protected,file.path(out,'audit/protected_source_fingerprints.csv'))
    writeLines(capture.output(sessionInfo()),file.path(out,'audit/sessionInfo.txt'))
    writeLines(c('COMPLETED',fc_protocol),file.path(out,'run_status.txt'))
    target<-fl_output_complete_run(ctx);writeLines(basename(target),file.path(ctx$module,'LATEST_RUN.txt'))
    cat('COMPLETED:',target,'\n');invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
