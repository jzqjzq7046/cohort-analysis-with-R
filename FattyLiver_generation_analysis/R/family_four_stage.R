# Inert definitions. Fresh, unrestricted within-family discovery; genuine strain
# inference is gated on sequence-derived evidence, never imputed from abundance.
f4_protocol <- 'within_family_discovery_strain_function_20261104_v1'
f4_root <- function() Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
f4_setup <- function(root) {
  sys.source(file.path(root,'R/abc_analysis.R'),envir=.GlobalEnv)
  cfg<-abc_setup(root)
  sys.source(file.path(root,'R/generation_ADE_extensions.R'),envir=.GlobalEnv)
  stopifnot(identical(as.character(Matrix::Matrix.Version()$abi),as.character(lme4:::.Matrix.abi.build.version)))
  # Distinct new checkpoint namespace and signature. Official matched-set LDM
  # code/transform/permutation rules are preserved, with no old discovery filter.
  assign('abc_protocol',f4_protocol,envir=.GlobalEnv)
  cfg
}
f4_refs <- function() data.frame(
  stage=c('Within-family design','Within-family inference','Strain profiling and sharing','Multilevel family microbiome','Functional profiling','Technical covariate rationale','Lifestyle/BMI model hierarchy'),
  reference=c('Chen W, Li Y, Wang W, et al. Enhanced microbiota profiling in patients with quiescent Crohn\'s disease through comparison with paired healthy first-degree relatives. Cell Reports Medicine. 2024;5:101624. doi:10.1016/j.xcrm.2024.101624',
    'Zhu Z, Satten GA, Mitchell C, Hu YJ. Constraining PERMANOVA and LDM to within-set comparisons by projection improves the efficiency of analyses of matched sets of microbiome data. Microbiome. 2021;9:133. doi:10.1186/s40168-021-01034-9',
    'Valles-Colomer M, et al. The person-to-person transmission landscape of the gut and oral microbiomes. Nature. 2023;614:125-135. doi:10.1038/s41586-022-05620-1',
    'Valles-Colomer M, Bacigalupe R, Vieira-Silva S, et al. Variation and transmission of the human gut microbiota across multiple familial generations. Nature Microbiology. 2022;7:87-96. doi:10.1038/s41564-021-01021-8',
    'Beghini F, et al. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. eLife. 2021;10:e65088. doi:10.7554/eLife.65088',
    'Gacesa R, Kurilshikov A, Vich Vila A, et al. Environmental factors shaping the gut microbiome in a Dutch population. Nature. 2022;604:732-739. doi:10.1038/s41586-022-04567-7',
    'Ma W, Wang Y, Nguyen LH, et al. Gut microbiome composition and metabolic activity in women with diverticulitis. Nature Communications. 2024;15:3612. doi:10.1038/s41467-024-47859-4'),
  used_for=c('Disease-discordant relatives as comparators; primary labels do not prove first-degree relation or cohabitation',
    'Unequal matched sets: family indicators, nuisance projection and permutations restricted within family; frequency/AST omnibus',
    'Requires sequence-derived per-SGB marker distances, species-specific thresholds and joint profileability',
    'Generational and familial context; abundance is not strain transmission or host heritability',
    'Existing actual HUMAnN DNA pathway abundance; functional potential, not expression or metabolite concentrations',
    'Age/sex, technical and seasonal covariate rationale; exact generation/pool/family design is this cohort adaptation',
    'Basic versus full lifestyle/BMI hierarchy adapted to measured variables; no substitution for unmeasured Bristol, hormone or race variables'))
f4_models <- function(x) {
  adult<-abc_discordant(x$metadata,v3_B2,adult=TRUE)
  list(B2=list(m=adult,fields=v3_B2,display='Adult Model 4',role='Primary'),
    Adult_B1_same_sample=list(m=adult,fields=v3_B1,display='Adult Model 3',role='Same-sample sensitivity'),
    B1=list(m=abc_discordant(x$metadata,v3_B1),fields=v3_B1,display='All-age Model 3',role='All-age sensitivity'))
}
f4_fit <- function(x,cfg,out,cache,kind) {
  stopifnot(!anyNA(x$metadata$participant_id),!anyDuplicated(x$metadata$participant_id))
  specs<-f4_models(x); fits<-list()
  for(k in names(specs)) {
    a<-specs[[k]]
    cat(format(Sys.time()),kind,a$display,': fresh within-family inference\n');flush.console()
    fits[[k]]<-abc_B_fit(x,a$m,a$fields,k,out,cache)
    saveRDS(fits,file.path(out,'private/incremental_fits.rds'))
  }
  e<-ov_bind(lapply(fits,function(a)a$results));s<-ov_bind(lapply(fits,function(a)a$summary))
  for(k in names(specs)) {
    i<-e$model==k;e$model_display[i]<-specs[[k]]$display;e$analysis_role[i]<-specs[[k]]$role
    j<-s$model==k;s$model_display[j]<-specs[[k]]$display;s$analysis_role[j]<-specs[[k]]$role
  }
  ann<-if(kind=='SGB')ade_labels(x$features)else x$features
  ann<-ann[,intersect(c('feature_id','feature_label','scientific_name','annotation_status','annotation_candidates','taxonomy','pathway_name'),names(ann)),drop=FALSE]
  e<-merge(e,ann,by='feature_id',all.x=TRUE,sort=FALSE)
  if(kind=='Pathway')e$feature_label<-paste(e$feature_id,e$pathway_name,sep=': ')
  e<-e[order(match(e$model,names(specs)),e$P_BH,e$P_raw,e$feature_id),]
  e$direction_interpretation<-'Descriptive conditional mean contrast; omnibus P combines two scales'
  s$global_P_BH_across_three_models<-p.adjust(s$global_P_omnibus,'BH')
  counts<-ov_bind(lapply(names(specs),function(k){m<-specs[[k]]$m;data.frame(model=k,model_display=specs[[k]]$display,
    N=nrow(m),cases=sum(m$fatty_liver),controls=sum(m$fatty_liver==0),families=nlevels(m$family_id),
    G0=sum(m$generation=='G0'),G1=sum(m$generation=='G1'),G2=sum(m$generation=='G2'),minors=sum(m$age<18))}))
  saveRDS(lapply(specs,function(a)a$m),file.path(out,'private/model_populations.rds'))
  list(effects=e,summary=s,counts=counts,annotation=ann,kind=kind,protocol=f4_protocol)
}
f4_verified_function_fit <- function(x,cfg,out,cache) {
  prep<-file.path(cfg$project_root,'outputs/48_Family_functional_explanation/_precomputed',f4_protocol)
  input<-file.path(prep,'private/functional_input.rds')
  if(!file.exists(input))return(f4_fit(x,cfg,out,cache,'Pathway'))
  oldinput<-readRDS(input);current<-readRDS(file.path(out,'private/functional_input.rds'))
  for(n in c('relative','metadata','features','candidates','denominator','source_files','protected'))
    stopifnot(identical(oldinput[[n]],current[[n]]))
  specs<-f4_models(x);approved<-list();checks<-list()
  for(k in names(specs)) {
    cp<-file.path(prep,'cache',paste0(k,'.rds'));design<-file.path(prep,'private',paste0(k,'_analysis_design.rds'))
    stopifnot(file.exists(cp),file.exists(design))
    previous<-readRDS(design);m<-specs[[k]]$m;Z<-abc_design(m,specs[[k]]$fields,as.character(m$family_id))$Z
    er<-drop(m$fatty_liver-Z%*%crossprod(Z,m$fatty_liver))
    stopifnot(identical(previous$metadata,m),identical(previous$Z,Z),identical(previous$disease_residual,er))
    A<-t(x$relative[,match(m$sample_id,colnames(x$relative)),drop=FALSE]);keep<-colMeans(A>0)>=.1&colSums(A>0)>=5
    ids<-colnames(A)[keep];v<-readRDS(cp);ans<-v$value$result
    stopifnot(identical(names(ans$p.otu.omni),ids),ans$n.perm.completed==10000,
      ans$seed==20261029L+match(k,c('B1','G0_same_generation','G2_same_generation','Adult_B1_same_sample','B2')),
      all(v$value$equivalence$relative_error<1e-10),as.character(packageVersion('LDM'))=='6.0.1',
      identical(v$value$official_statistic_hash,digest::digest(body(get('ldm.stat',asNamespace('LDM'))),algo='sha256')),
      identical(v$value$official_ldm_hash,digest::digest(body(get('ldm',asNamespace('LDM'))),algo='sha256')))
    approved[[paste0(k,'.rds')]]<-v
    checks[[k]]<-data.frame(model=k,full_input_values_and_attributes_identical=TRUE,
      participants_and_order_identical=TRUE,nuisance_matrix_identical=TRUE,exposure_residual_identical=TRUE,
      tested_features_and_order_identical=TRUE,official_functions_identical=TRUE,
      permutations=ans$n.perm.completed,seed=ans$seed,checkpoint_MD5=unname(tools::md5sum(cp)))
  }
  # Value-level identity is stronger than a session-dependent serialized hash.
  # Reuse only this protocol's independently completed official results after
  # every input, design, seed, permutation count and function hash check above.
  # No generic checkpoint rule or old analysis module is changed.
  original<-get('tf_checkpoint',envir=.GlobalEnv);logs<-list()
  on.exit(assign('tf_checkpoint',original,envir=.GlobalEnv),add=TRUE)
  verified<-function(path,signature,compute) {
    nm<-basename(path)
    stopifnot(identical(normalizePath(dirname(path),winslash='/',mustWork=TRUE),normalizePath(cache,winslash='/',mustWork=TRUE)),nm%in%names(approved))
    v<-approved[[nm]]
    logs[[nm]]<<-data.frame(checkpoint=nm,stored_signature=v$signature,current_signature=signature,
      serialized_hash_equal=identical(v$signature,signature),action='Reuse numerically identical verified official precomputation; no refit')
    v$value
  }
  assign('tf_checkpoint',verified,envir=.GlobalEnv)
  r<-f4_fit(x,cfg,out,cache,'Pathway')
  fl_write_csv(ov_bind(checks),file.path(out,'audit/precomputed_model_input_equivalence.csv'))
  fl_write_csv(ov_bind(logs),file.path(out,'audit/precomputed_checkpoint_reuse.csv'))
  r
}
f4_validate <- function(r,out) {
  checks<-list()
  for(k in unique(r$effects$model)) {
    e<-r$effects[r$effects$model==k,];d<-readRDS(file.path(out,'private',paste0(k,'_analysis_design.rds')))
    eq<-read.csv(file.path(out,'audit',paste0(k,'_official_statistic_equivalence.csv')))
    stopifnot(!anyDuplicated(e$feature_id),!anyDuplicated(d$metadata$participant_id),max(abs(p.adjust(e$P_raw,'BH')-e$P_BH))<1e-12,
      all(tapply(d$metadata$fatty_liver,d$metadata$family_id,function(z)length(unique(z)))==2),
      all(eq$relative_error<1e-10),max(abs(crossprod(d$Z,d$disease_residual)))<1e-7)
    a<-readRDS(file.path(out,'private',paste0(k,'_official_LDM.rds')))
    pp<-a$p.otu.omni;pp<-if(is.matrix(pp))pp[,ncol(pp)]else pp
    stopifnot(max(abs(sort(as.numeric(pp))-sort(e$P_raw)))<1e-12)
    checks[[k]]<-data.frame(model=k,tested=nrow(e),BH_recomputed=TRUE,official_P_exact=TRUE,
      person_unique=!anyDuplicated(d$metadata$participant_id),sample_unique=!anyDuplicated(d$metadata$sample_id),within_family_permutation=TRUE,
      maximum_equivalence_error=max(eq$relative_error),residual_df=nrow(d$Z)-ncol(d$Z)-1)
  }
  mm<-readRDS(file.path(out,'private/model_populations.rds'))
  stopifnot(identical(mm$B2$sample_id,mm$Adult_B1_same_sample$sample_id))
  fl_write_csv(ov_bind(checks),file.path(out,'audit/statistical_validation.csv'))
  TRUE
}
f4_tables <- function(r,out) {
  v3_table(r$effects,'All_tested_features_EN',out)
  v3_table(r$effects[r$effects$BH_lt_005,],'All_BH_discoveries_EN',out)
  v3_table(r$effects[r$effects$model=='B2'&r$effects$BH_lt_005,],'Primary_adult_Model4_discoveries_EN',out)
  v3_table(r$summary,'Model_and_global_summary_EN',out);v3_table(r$counts,'Analysis_populations_EN',out)
  v3_table(r$annotation,'Complete_annotation_EN',out);v3_table(f4_refs(),'References_and_adaptations_EN',out)
  mm<-readRDS(file.path(out,'private/model_populations.rds'))
  base<-readRDS(r$protected$path[1])$metadata
  flow<-ov_bind(lapply(names(mm),function(k){
    adult<-k!='B1';fields<-if(k=='B2'||k=='Adult_B1_same_sample')v3_B2 else v3_B1
    ageok<-if(adult)!is.na(base$age)&base$age>=18 else rep(TRUE,nrow(base))
    cc<-ageok&v3_cc(base,c('family_id','fatty_liver',fields))
    data.frame(model=k,source_participants=nrow(base),age_eligible=sum(ageok),
      complete_covariates=sum(cc),after_discordant_family_filter=nrow(mm[[k]]),
      excluded_non_discordant_complete=sum(cc)-nrow(mm[[k]]),
      adult_complete_sample_aligned=k=='Adult_B1_same_sample')
  }))
  v3_table(flow,'Sample_flow_EN',out)
  descriptives<-ov_bind(lapply(names(mm),function(k){m<-mm[[k]]
    ov_bind(lapply(split(m,list(m$generation,m$fatty_liver),drop=TRUE),function(d)
      data.frame(model=k,generation=as.character(d$generation[1]),steatosis=d$fatty_liver[1],N=nrow(d),
        age_min=min(d$age),age_Q1=quantile(d$age,.25),age_median=median(d$age),age_Q3=quantile(d$age,.75),age_max=max(d$age))))}))
  v3_table(descriptives,'Age_and_generation_overlap_EN',out)
}
f4_report <- function(r,out) {
  lines<-c('# Within-family de novo discovery and functional follow-up',
    '',paste('Protocol:',f4_protocol),paste('Analyte:',r$kind),'',
    'All eligible features were tested afresh; the previous overall-steatosis discovery list was not used to select inputs.',
    'Primary: adult complete-case Model 4, in families with both ultrasound steatosis states. Adult Model 3 uses exactly the same people; all-age Model 3 is a separate sensitivity analysis.',
    'Model 3: generation, generation-specific age splines, sex, library batch, sequencing pool, cyclic season and recent antibiotics. Model 4 additionally includes smoking, alcohol, diet score, activity, log energy and continuous BMI.',
    'Family indicators absorb shared family-level factors. This does not remove unmeasured individual confounding or establish causation.',
    'Model 4 is conditional on BMI and lifestyle; it is not an estimate of the total causal effect of steatosis. Differences between the adult models cannot be attributed to BMI alone or interpreted as mediation because several variables are added together.',
    'Official LDM 6.0.1: frequency and arcsine-square-root omnibus, up to 10,000 within-family permutations, sequential stopping after 100 exceedances, no between-family permutations; retained rows are not renormalized.',
    'Eligibility: pooled detection >=10% and >=5 carriers in each analysis population. BH is recalculated across every tested feature within each model; the three global omnibus P values are also BH-adjusted.',
    'The reported P values are Monte Carlo/sequential permutation estimates, not exact enumeration. With 10,000 draws the smallest plus-one P is 1/10,001; discoveries at this floor have limited numerical resolution. Official LDM q values are also retained; they can differ from conventional BH applied to reported omnibus P values. The declared discovery criterion here is conventional BH q <0.05.',
    'Effects are descriptive adjusted frequency differences (percentage points) and AST contrasts. A combined omnibus test does not have a unique effect scale. Opposite signs across the scales are explicitly marked.',
    'SGBs are species-level genome bins, not strain types. Primary-family labels do not prove first-degree relationships, shared homes or directional transmission.',
    '',capture.output(print(r$summary,row.names=FALSE)),
    '',if(r$kind=='Pathway')'Pathways are measured HUMAnN DNA functional potential. Community profiles use the fixed sum of all supplied unstratified rows before filtering. Neither microbial expression, metabolite output nor strain-specific mechanism is measured.' else '',
    if(r$kind=='Pathway')'The original HUMAnN version and raw measurement units are not independently verified. These are relative pathway profiles within the supplied community matrix, not total metagenomic functional coverage.'else'',
    '', 'Strain typing, actual cross-generation strain sharing and sharing versus both-member disease status require sequence/marker data. Missing strain inputs are recorded separately and never replaced by species co-detection or abundance concordance.',
    '',apply(f4_refs(),1,function(z)paste0('- ',z['reference'],'; ',z['used_for'])))
  writeLines(lines,file.path(out,'reports/Methods_results_and_limits_EN.md'),useBytes=TRUE)
  writeLines(c(paste0('Figure ',if(r$kind=='SGB')'45'else'48',' | Within-family associations of ultrasound hepatic steatosis with ',if(r$kind=='SGB')'SGB abundance.'else'microbial functional potential.'),
    '(a) All eligible adult Model 4 features. The horizontal axis is the adjusted within-family difference in relative abundance (percentage points); the vertical axis is -log10(BH q) from the official LDM frequency/AST omnibus. Blue, lower; orange, higher; grey, not BH significant; dark grey, scale-dependent direction. The dashed horizontal line marks q=0.05.',
    '(b) Numbers of BH q<0.05 features in all-age Model 3, adult Model 3 and adult Model 4. The two adult models use exactly the same participants. The all-age population is different and is not a nested adjustment comparison.',
    'Primary-discovery heatmaps show every adult Model 4 discovery, paginated only for readability. Cells are adjusted frequency differences in percentage points; * denotes BH q<0.05 within that model and the complete tested feature set. Each figure page uses one symmetric colour scale, with no per-row scaling.',
    if(r$kind=='SGB')'Family-contrast plots show unadjusted mean abundance for affected and unaffected members in each retained primary family. Lines connect the two family-state means; grey points are family-state means, not individual participants. Box outlines are black; boxes show the median and IQR, whiskers extend to observations within 1.5 IQR. These descriptive plots do not replace the adjusted inference.'else'Only primary-significant pathways receive a discovery heatmap. If none pass BH, no selected-pathway discovery plot is produced.',
    'Primary models include family indicators and age, sex, generation, batch, pool, season, antibiotics, smoking, alcohol, diet, activity, energy intake and BMI. See the methods report for exact definitions, sample flow, permutation resolution and references.'),
    file.path(out,'reports/Figure_legends_EN.txt'),useBytes=TRUE)
}
f4_figures <- function(r,x,cfg,out) {
  # Expand shared helpers for static QA: the actual final device/font/alignment
  # contract lives in tf_export, not in each caller. Native R parsing is also
  # required because heuristic scanners do not parse R function-name syntax.
  capsule<-file.path(out,'audit/figure_source_capsule.R')
  funcs<-c('nf_theme','tf_export','f4_figures')
  writeLines(unlist(lapply(funcs,function(n)paste0(n,' <- ',paste(deparse(get(n,envir=.GlobalEnv),width.cutoff=120),collapse='\n')))),capsule)
  invisible(parse(capsule,encoding='UTF-8'))
  system2(cfg$python,c(shQuote(file.path(nf_skill_root(),'scripts/validate_figure.py')),shQuote(capsule)),
    stdout=file.path(out,'audit/figure_static_preflight.txt'),stderr=file.path(out,'audit/figure_static_preflight_stderr.txt'))
  e<-r$effects[r$effects$model=='B2',,drop=FALSE]
  stopifnot(all(is.finite(e$P_BH)),all(e$P_BH>0&e$P_BH<=1))
  e$colour<-ifelse(!e$BH_lt_005,'Not significant',ifelse(!e$scale_direction_agrees,'Scale dependent',ifelse(e$beta_relative<0,'Lower','Higher')))
  pal<-c('Lower'='#79CAFB','Higher'='#F7A24F','Not significant'='#BBBBBB','Scale dependent'='#666666')
  p<-ggplot2::ggplot(e,ggplot2::aes(beta_percentage_points,-log10(P_BH)))+
    ggplot2::geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.3,colour='#888888')+
    ggplot2::geom_vline(xintercept=0,linetype=2,linewidth=.3,colour='#888888')+
    ggplot2::geom_point(ggplot2::aes(colour=colour),size=1.25,alpha=.8)+
    ggplot2::scale_colour_manual(values=pal,drop=FALSE)+
    ggplot2::expand_limits(x=if(r$kind=='SGB')floor(min(e$beta_percentage_points))else min(e$beta_percentage_points))+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.04,.15)))+
    ggplot2::labs(x='Adjusted within-family difference (percentage points)',y='-log10(BH q)',colour=NULL)
  if(sum(e$BH_lt_005)>0&&sum(e$BH_lt_005)<=6) {
    lab<-e[e$BH_lt_005,];near<-abs(lab$beta_percentage_points)<.2
    lab$label_x<-lab$beta_percentage_points-ifelse(near,.03,0);lab$label_hjust<-ifelse(near,1,.5)
    p<-p+ggplot2::geom_text(data=lab,ggplot2::aes(x=label_x,label=sub('^t__','',feature_id),hjust=label_hjust),
      nudge_y=.10,size=6/ggplot2::.pt,colour='black')
  }
  s<-r$summary;s$model_display<-factor(s$model_display,c('All-age Model 3','Adult Model 3','Adult Model 4'))
  q<-ggplot2::ggplot(s,ggplot2::aes(model_display,BH_lt_005))+
    ggplot2::geom_col(fill='#AAAAAA',colour='black',linewidth=.25,width=.55)+
    ggplot2::geom_text(ggplot2::aes(label=BH_lt_005),vjust=-.6,size=6/ggplot2::.pt)+
    ggplot2::scale_x_discrete(labels=c('All-age Model 3'='All ages\nModel 3','Adult Model 3'='Adults\nModel 3','Adult Model 4'='Adults\nModel 4'))+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.16)))+
    ggplot2::labs(x=NULL,y='BH q < 0.05 features')
  axis_spacing<-ggplot2::theme(axis.title.x=ggplot2::element_text(margin=ggplot2::margin(t=7,unit='pt')),
    axis.text.x=ggplot2::element_text(margin=ggplot2::margin(t=3,unit='pt')))
  fig<-patchwork::wrap_plots(list(p+tf_plot_theme()+axis_spacing+ggplot2::labs(tag='a'),
    q+tf_plot_theme()+axis_spacing+ggplot2::labs(tag='b')),nrow=1,guides='collect')+
    patchwork::plot_annotation(theme=ggplot2::theme(legend.position='bottom'))
  tf_export(list(plot=fig,width=183,height=92,panels=2,exemptions=list()),paste0('Figure',if(r$kind=='SGB')'45'else'48','_Within_family_summary'),out,cfg)
  sig<-e[e$BH_lt_005,];sig<-sig[order(sig$P_BH,sig$feature_id),]
  if(!nrow(sig))return(invisible(NULL))
  # All primary discoveries, paginated for readable physical text; no top-N filter.
  for(pg in seq_len(ceiling(nrow(sig)/24))) {
    ids<-sig$feature_id[seq.int((pg-1)*24+1,min(pg*24,nrow(sig)))];z<-r$effects[r$effects$feature_id%in%ids,]
    labels<-sig$feature_label[match(ids,sig$feature_id)];labels<-vapply(labels,function(a)paste(strwrap(a,width=48),collapse='\n'),character(1))
    z$feature<-factor(z$feature_id,rev(ids),rev(labels));z$model_display<-factor(z$model_display,c('All-age Model 3','Adult Model 3','Adult Model 4'))
    lim<-max(abs(z$beta_percentage_points));if(!is.finite(lim)||lim==0)lim<-1
    h<-ggplot2::ggplot(z,ggplot2::aes(model_display,feature,fill=beta_percentage_points))+
      ggplot2::geom_tile(colour='white',linewidth=.5)+
      ggplot2::geom_text(ggplot2::aes(label=if(nrow(sig)<=8)paste0(format(signif(beta_percentage_points,2),trim=TRUE),ifelse(BH_lt_005,'*',''))else ifelse(BH_lt_005,'*','')),size=5.7/ggplot2::.pt)+
      ggplot2::scale_fill_gradient2(low='#79CAFB',mid='white',high='#F7A24F',limits=c(-lim,lim),
        name='Difference (pp)',breaks=c(-lim,0,lim),labels=function(a)format(signif(a,2),trim=TRUE),
        guide=ggplot2::guide_colourbar(display='rectangles',barwidth=grid::unit(3,'mm'),barheight=grid::unit(if(nrow(sig)<=8)12 else 30,'mm'),title.position='top'))+
      ggplot2::scale_x_discrete(labels=c('All-age Model 3'='All ages\nModel 3','Adult Model 3'='Adults\nModel 3','Adult Model 4'='Adults\nModel 4'))+
      ggplot2::labs(x=NULL,y=NULL)+tf_plot_theme()+ggplot2::theme(legend.position='right',legend.title=ggplot2::element_text(size=5.7),axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank())
    tf_export(list(plot=patchwork::wrap_plots(h),width=if(nrow(sig)<=8)89 else 183,height=if(nrow(sig)<=8)max(38,17+length(ids)*4.5)else max(55,13+length(ids)*6),panels=1,
      no_panel_tags=TRUE,exemptions=list()),sprintf('Figure%s_Primary_discoveries_page%02d',if(r$kind=='SGB')'45'else'48',pg),out,cfg)
  }
  if(r$kind=='SGB') {
    mm<-readRDS(file.path(out,'private/model_populations.rds'))$B2
    for(pg in seq_len(ceiling(nrow(sig)/4))) {
      ids<-sig$feature_id[seq.int((pg-1)*4+1,min(pg*4,nrow(sig)))];plots<-list()
      for(id in ids) {
        dat<-data.frame(family=mm$family_id,case=mm$fatty_liver,RA=100*x$relative[id,match(mm$sample_id,colnames(x$relative))])
        d<-aggregate(RA~family+case,dat,mean);d$group<-factor(d$case,0:1,c('No steatosis','Steatosis'))
        plots[[id]]<-ggplot2::ggplot(d,ggplot2::aes(group,RA))+
          ggplot2::geom_line(ggplot2::aes(group=family),colour='#777777',alpha=.16,linewidth=.22)+
          ggplot2::geom_boxplot(ggplot2::aes(fill=group),width=.45,colour='black',outlier.shape=NA,linewidth=.3)+
          ggplot2::geom_point(colour='#777777',size=.3,alpha=.25,position=ggplot2::position_jitter(width=.06,seed=45))+
          ggplot2::scale_fill_manual(values=c('No steatosis'='#79CAFB','Steatosis'='#F7A24F'))+
          ggplot2::labs(title=sub('^t__','',id),x=NULL,y='Family-state mean relative abundance (%)')
      }
      tf_export(list(plot=abc_wrap(plots),width=183,height=94,panels=length(plots),exemptions=list()),
        sprintf('Figure45_Family_contrasts_page%02d',pg),out,cfg)
    }
  }
}
f4_function_context <- function(r,discovery,cfg,out) {
  path<-file.path(cfg$project_root,'data/reference/pathway_redesign_20261102/map_metacyc-pwy_lineage.tsv')
  if(file.exists(path)) {
    a<-read.delim(path,header=FALSE,sep='\t',quote='',comment.char='',col.names=c('feature_id','lineage'))
    z<-merge(r$effects,a,by='feature_id',all.x=TRUE,sort=FALSE)
    z$lineage[is.na(z$lineage)]<-'Unmapped'
    v3_table(z,'All_pathway_results_and_historical_MetaCyc_lineages_EN',out)
  }
  # Reuse only immutable measured source annotations; no old disease estimates.
  mod<-file.path(cfg$project_root,'outputs/44_Pathway_function_sources');p<-file.path(mod,'LATEST_RUN.txt')
  if(!file.exists(p))return(r)
  f<-file.path(mod,'current',readLines(p,1),'tables/HUMAnN_pathway_source_existing_SGB_context_EN.csv')
  if(!file.exists(f))return(r)
  a<-read.csv(f,check.names=FALSE)
  a<-unique(a[,c('taxon_label','source_row','feature_id_pathway','source_label','genus','species',
    'humann_genus','humann_species','exact_named_SGB_candidates','mapping_status','feature_id_SGB')])
  saveRDS(a,file.path(out,'private/existing_source_annotation.rds'))
  # Only immutable annotation is reused; prior clinical coefficients/P/q are removed.
  dmod<-file.path(cfg$project_root,'outputs/45_Family_de_novo_discovery')
  dp<-file.path(dmod,'LATEST_RUN.txt')
  if(file.exists(dp)) {
    dr<-readRDS(file.path(dmod,'current',readLines(dp,1),'private/analysis_result.rds'))
    ds<-dr$effects[dr$effects$model=='B2'&dr$effects$BH_lt_005,]
    exact<-a[a$mapping_status=='Unique exact binomial match'&!is.na(a$feature_id_SGB)&nzchar(a$feature_id_SGB),]
    cov<-ov_bind(lapply(seq_len(nrow(ds)),function(i){z<-exact[exact$feature_id_SGB==ds$feature_id[i],]
      data.frame(feature_id=ds$feature_id[i],feature_label=ds$feature_label[i],
        new_primary_SGB_q=ds$P_BH[i],exact_source_taxa=length(unique(z$taxon_label)),
        measured_source_pathways=length(unique(z$feature_id_pathway)),
        status=if(nrow(z))'Exact named annotation; measured source data present'else'No unique exact species mapping; no functional attribution made')}))
    if(nrow(cov))v3_table(cov,'New_primary_SGB_to_HUMAnN_mapping_coverage_EN',out)
    links<-merge(exact,ds[,c('feature_id','feature_label','P_BH')],by.x='feature_id_SGB',by.y='feature_id',sort=FALSE)
    if(nrow(links)) {
      pe<-r$effects[r$effects$model=='B2',c('feature_id','P_BH','beta_percentage_points','BH_lt_005')]
      names(pe)<-c('feature_id_pathway','new_primary_pathway_q','pathway_beta_percentage_points','pathway_BH_lt_005')
      links<-merge(links,pe,by='feature_id_pathway',all.x=TRUE,sort=FALSE)
      links$interpretation<-'Exact species annotation with observed HUMAnN stratum; not strain-specific attribution or mediation'
      v3_table(links,'New_discovery_SGBs_measured_pathway_annotation_links_EN',out)
    }
    r$discovery_source_mapping<-cov
  }
  v3_table(data.frame(status='Existing HUMAnN taxon-to-SGB annotation retained privately; no unverified strain-specific functional attribution',source=f,md5=unname(tools::md5sum(f))),
    'Functional_source_attribution_scope_EN',out)
  r
}
f4_run <- function(stage,args=character(),root=f4_root()) {
  if(getRversion()<'4.5') {
    entry<-c(discovery='45_family_de_novo_discovery.R',strain='46_family_strain_typing.R',disease='47_strain_sharing_disease.R',functionality='48_family_functional_explanation.R')[[stage]]
    st<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts',entry)),vapply(args,shQuote,character(1))))
    if(st!=0)stop('Project-R stage failed');return(invisible(NULL))
  }
  cfg<-f4_setup(root)
  if(stage%in%c('strain','disease'))return(f4_input_gate(stage,cfg,args))
  key<-if(stage=='discovery')'family_discovery'else'family_function'
  cfg$output_module<-unname(cfg$output_modules[key]);cfg$run_id<-unname(cfg$run_ids[key]);cfg$run_id_explicit<-FALSE
  stopifnot(all(args%in%c('--run','--figures-only','--resume')))
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(d in c('private','audit','tables','figures','reports'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
  writeLines(c('RUNNING',f4_protocol),file.path(out,'run_status.txt'))
  tryCatch({
    src<-abc_source(root);protected<-data.frame(path=src$file,md5=src$md5);x<-src$x
    if(any(args%in%c('--figures-only','--resume'))) {
      paths<-list.files(ctx$module,pattern='analysis_result.rds',recursive=TRUE,full.names=TRUE)
      paths<-paths[!startsWith(normalizePath(paths,winslash='/'),normalizePath(out,winslash='/'))]
      stopifnot(length(paths)>0);old<-paths[which.max(file.info(paths)$mtime)];r<-readRDS(old)
      stopifnot(identical(r$protocol,f4_protocol));prev<-dirname(dirname(old))
      for(d in c('private','audit','tables','reports')) {
        fs<-list.files(file.path(prev,d),full.names=TRUE);if(length(fs))stopifnot(all(file.copy(fs,file.path(out,d),recursive=TRUE,overwrite=TRUE)))
      }
      if(stage=='functionality') {
        r<-f4_function_context(r,NULL,cfg,out)
        saveRDS(r,file.path(out,'private/analysis_result.rds'))
      }
    }else{
      if(stage=='functionality') {
        imp<-ade_import_paths(cfg,x,out);x<-imp$x;protected<-ov_bind(list(protected,imp$protected))
        x$features$official_id<-x$features$feature_id;x$features$taxonomy<-''
        # Exclude only explicit nonbiological rows, without re-closing relative abundance.
        x$relative<-x$relative[x$candidates,,drop=FALSE]
      }
      cache<-file.path(ctx$module,'_checkpoints',f4_protocol,'private');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
      r<-if(stage=='discovery')f4_fit(x,cfg,out,cache,'SGB')else f4_verified_function_fit(x,cfg,out,cache)
      r$protected<-protected
      if(stage=='functionality')r<-f4_function_context(r,NULL,cfg,out)
      saveRDS(r,file.path(out,'private/analysis_result.rds'))
    }
    f4_tables(r,out);f4_report(r,out);f4_validate(r,out);f4_figures(r,x,cfg,out)
    ade_protect(r$protected);writeLines(capture.output(sessionInfo()),file.path(out,'audit/sessionInfo.txt'))
    writeLines(c('COMPLETED',f4_protocol),file.path(out,'run_status.txt'));target<-fl_output_complete_run(ctx)
    writeLines(basename(target),file.path(ctx$module,'LATEST_RUN.txt'));cat('COMPLETED:',target,'\n');invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
f4_input_gate <- function(stage,cfg,args=character()) {
  key<-if(stage=='strain')'family_strain'else'family_strain_disease'
  cfg$output_module<-unname(cfg$output_modules[key]);cfg$run_id<-unname(cfg$run_ids[key]);cfg$run_id_explicit<-FALSE
  ctx<-fl_output_start_run(cfg);out<-ctx$run_dir
  for(d in c('tables','reports','audit','private'))dir.create(file.path(out,d),showWarnings=FALSE)
  status<-data.frame(stage=if(stage=='strain')'2: Strain typing and cross-generation sharing'else'3: Sharing and both-member steatosis states',
    status='BLOCKED_INPUT',reason=if(stage=='strain')'No verified raw metagenomes, sample marker sequences or per-SGB genetic distances supplied'else'No genuine strain-sharing result available from stage 2',
    abundance_substitution=FALSE,transmission_or_heritability_claim=FALSE)
  fl_write_csv(status,file.path(out,'tables/Stage_status_EN.csv'))
  writeLines(c('# Strain evidence required',status$reason,
    'No species abundance, co-detection or family concordance is relabelled as strain sharing.',
    'Required for stage 2: sample-linked FASTQ/alignments or extracted markers; exact reference and profiler versions; marker QC, per-SGB distances and original-paper-compatible species-specific sharing thresholds.',
    'Required for stage 3: actual shared/not-shared calls restricted to jointly profileable sample-SGB pairs, both member IDs and ultrasound states, family IDs, valid relatedness/cohabitation information and individual covariates.',
    'Strain-sharing rate denominator: jointly profileable SGBs, not all detected species. Species-specific missing profiles are missing, not unshared.',
    'Use all eligible cross-generation pairs, including disease-concordant families; do not inherit the discordant-family discovery restriction. Profile background SGBs as well as discovery candidates. A two-species target list alone cannot provide a >=10-SGB person-pair sharing rate.',
    'Pair states planned: neither affected, only older affected, only younger affected, both affected. Primary labels do not prove direction of transmission or inheritance.',
    'Stage 3 prespecification: species-level family-clustered logistic GEE on actual shared-strain calls with both-member age/sex, sequencing technical factors, season and antibiotic use; sparse/separated species not forced. BH across all eligible species. Person-level sharing rates require >=10 jointly profiled SGBs as in Valles-Colomer 2023. Multiple pairs per individual require node/family-aware inference; that implementation remains pending genuine input and its pair structure.',
    'References: Valles-Colomer et al. Nature 2023;614:125-135. doi:10.1038/s41586-022-05620-1; Nature Microbiology 2022;7:87-96. doi:10.1038/s41564-021-01021-8.',
    'This is an input-gate report; neither a completed strain analysis nor a null strain result.'),file.path(out,'reports/Missing_inputs_and_resume_requirements_EN.md'))
  writeLines(c('BLOCKED_INPUT',f4_protocol),file.path(out,'run_status.txt'))
  # Blocked analyses are explicitly not adopted as completed numerical modules.
  report<-file.path(cfg$project_root,'Results_表格与图片/01_当前分析/20_家系内重新发现与功能解释/00_菌株分析缺失资料')
  dir.create(report,recursive=TRUE,showWarnings=FALSE)
  file.copy(file.path(out,'reports/Missing_inputs_and_resume_requirements_EN.md'),file.path(report,paste0(stage,'_missing_inputs_EN.md')),overwrite=TRUE)
  file.copy(file.path(out,'tables/Stage_status_EN.csv'),file.path(report,paste0(stage,'_status_EN.csv')),overwrite=TRUE)
  cat('BLOCKED_INPUT:',stage,'; missing-evidence report saved\n');invisible(out)
}
