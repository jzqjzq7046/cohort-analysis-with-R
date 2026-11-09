# Inert definitions. Genus disease discovery followed by descriptive family mapping.
# Relative abundances remain outcomes; abundance never determines display selection.
sg_protocol <- 'steatosis_genus_Maaslin2_LOG_all_BH005_no_display_cap_20261031_v1'
sg_labels <- c(B1='Model 3: all ages', Adult_B1_same_sample='Model 3: adults',
  B2='Model 4: adults', No_recent_antibiotics='Model 3: no recent antibiotics')

sg_setup <- function(root) {
  sys.source(file.path(root,'R/genus_generation_taxonomy.R'), envir=.GlobalEnv)
  cfg <- gg_setup(root)
  cfg <- lr_setup(root)
  stopifnot(as.character(packageVersion('Maaslin2')) == '1.24.1')
  cfg
}

sg_inputs <- function(cfg) {
  p24 <- lr_adopted(cfg$project_root,'24_Two_generation_reanalysis')
  p36 <- lr_adopted(cfg$project_root,'36_Genus_generation_taxonomy')
  paths <- c(file.path(p24,'private/stage_prepare.rds'),
    file.path(p36,'private/genus_relative_abundance.rds'),
    file.path(p36,'private/analysis_result.rds'),
    file.path(p36,'tables/All_genus_cross_generation_GEE_results_EN.csv'),
    file.path(cfg$project_root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),
    file.path(cfg$project_root,'outputs/36_Genus_generation_taxonomy/LATEST_RUN.txt'))
  list(paths=paths, protected=lr_hashes(paths), p24=p24, p36=p36)
}

sg_design <- function(x, mask, fields) {
  d <- v3_design(x,mask,fields)
  variance <- apply(d$relative,2,var)
  d$filter$variance <- unname(variance)
  d$filter$retained <- d$filter$retained & is.finite(variance) & variance>0
  d$filter$reason <- ifelse(d$filter$prevalence<.10,'FILTERED_PREVALENCE_LT_10_PERCENT',
    ifelse(!is.finite(variance)|variance<=0,'FILTERED_NO_VARIATION','RETAINED'))
  d
}

sg_fit <- function(d, model, x, cfg, cache, out) {
  ids <- d$filter$feature_id[d$filter$retained]
  if (!length(ids)) stop('No eligible genera: ',model)
  sig <- digest::digest(list(protocol=sg_protocol, model=model, A=d$relative[,ids,drop=FALSE],
    design=d$design, extraction=body(ov_extract),
    Maaslin2=as.character(packageVersion('Maaslin2'))),algo='sha256')
  tf_log(sg_labels[model],': ',nrow(d$metadata),' participants; ',sum(d$metadata$fatty_liver),
    ' cases; ',length(ids),' eligible genera')
  fit <- tf_official(d,ids,model,cfg,file.path(cache,model),sig)
  z <- ov_expand(fit$effects,d,x$candidates,model)
  at <- match(z$feature_id,d$filter$feature_id)
  omitted <- !d$filter$retained[at]
  z$status[omitted] <- d$filter$reason[at][omitted]
  t <- data.frame(model=model, model_label=unname(sg_labels[model]), feature_id=z$feature_id,
    genus=sub('^g__','',z$feature_id), N=nrow(d$metadata),cases=sum(d$metadata$fatty_liver),
    families=length(unique(d$metadata$family_id)),prevalence=d$filter$prevalence[at],
    beta_log2=z$estimate,SE=z$SE,df=z$df,CI_low=z$CI_low,CI_high=z$CI_high,P_raw=z$p_raw,
    analysis_valid=z$analysis_valid,status=z$status,singular=z$singular)
  t <- tf_effect_units(tf_adjust(t,paste0('Disease_',model,'_all_eligible_genera')))
  stopifnot(nrow(t)==length(x$candidates),!anyDuplicated(t$feature_id))
  fl_write_csv(d$mapping,file.path(out,'audit',paste0(model,'_nuisance_design.csv')))
  fl_write_csv(d$filter,file.path(out,'audit',paste0(model,'_genus_filter.csv')))
  fl_write_csv(fit$diagnostics,file.path(out,'audit',paste0(model,'_model_diagnostics.csv')))
  replacements <- vapply(ids,function(id) min(d$relative[d$relative[,id]>0,id])/2,numeric(1))
  fl_write_csv(data.frame(feature_id=ids,zero_replacement=replacements),
    file.path(out,'audit',paste0(model,'_zero_replacements.csv')))
  saveRDS(list(sample_id=d$metadata$sample_id,family_id=d$metadata$family_id,
    fields=if(model=='B2')v3_B2 else v3_B1,formula=d$formula),
    file.path(out,'private',paste0(model,'_analysis_sample.rds')))
  t
}

sg_select <- function(t, taxonomy) {
  s <- t[t$FDR_lt_005 & tf_bool(t$analysis_valid),,drop=FALSE]
  s <- s[order(s$P_BH,s$P_raw,s$feature_id),,drop=FALSE]
  s$association_evidence_rank <- seq_len(nrow(s))
  cols <- setdiff(names(taxonomy),c('genus','main_display','feature_id'))
  for (nm in cols) s[[nm]] <- taxonomy[[nm]][match(s$feature_id,taxonomy$feature_id)]
  s$selection_rule <- rep('All valid all-age Model 3 genera with full-family BH q < 0.05; no count cap',nrow(s))
  s$family_key <- paste(s$phylum,s$family,sep='|')
  # Taxonomic display order is independent of mean abundance and pair results.
  s <- s[order(s$phylum,s$family,s$P_BH,s$P_raw,s$feature_id),,drop=FALSE]
  s$display_row <- seq_len(nrow(s))
  s
}

sg_build <- function(cfg,out,cache) {
  inp <- sg_inputs(cfg)
  x <- readRDS(inp$paths[1]);a <- readRDS(inp$paths[2]);g <- readRDS(inp$paths[3])
  stopifnot(identical(x$protocol,v3_protocol),identical(g$protocol,gg_protocol),
    identical(colnames(a$relative),as.character(x$metadata$sample_id)),
    identical(a$metadata$sample_id,x$metadata$sample_id),
    all(rowSums(table(x$metadata$family_id,x$metadata$generation)>0)>=2),
    all(a$relative>=0),all(is.finite(a$relative)),max(colSums(a$relative))<=1+1e-10)
  lr_check_hashes(g$protected)
  saved <- tf_read(inp$paths[4]);stopifnot(isTRUE(all.equal(saved,g$results,
    tolerance=1e-12,check.attributes=FALSE)))
  stopifnot(setequal(rownames(a$relative),g$taxonomy$feature_id))
  x$relative <- a$relative;x$candidates <- rownames(a$relative)
  specs <- v3_micro_samples(x)
  specs$No_recent_antibiotics <- list(mask=specs$B1$mask & !is.na(x$metadata$abx_recent) &
    x$metadata$abx_recent==0,fields=v3_B1)
  tables <- samples <- list();adult_ids <- NULL
  for (model in names(specs)) {
    d <- sg_design(x,specs[[model]]$mask,specs[[model]]$fields)
    if (model=='Adult_B1_same_sample') adult_ids <- d$metadata$sample_id
    if (model=='B2') stopifnot(identical(adult_ids,d$metadata$sample_id))
    tables[[model]] <- sg_fit(d,model,x,cfg,cache,out)
    z <- tables[[model]]
    samples[[model]] <- data.frame(model=model,model_label=unname(sg_labels[model]),N=nrow(d$metadata),
      cases=sum(d$metadata$fatty_liver),families=length(unique(d$metadata$family_id)),
      total_genus_groups=nrow(z),eligible=sum(d$filter$retained),valid=sum(z$analysis_valid),
      singular_valid=sum(z$analysis_valid & tf_bool(z$singular)),
      BH_supported=sum(z$FDR_lt_005),lower=sum(z$FDR_lt_005 & z$beta_log2<0),
      higher=sum(z$FDR_lt_005 & z$beta_log2>0),
      constant_or_aliased_terms=paste(d$dropped,collapse=';'))
    v3_table(z,paste0('Disease_',model,'_all_genera_EN'),out)
  }
  selected <- sg_select(tables$B1,g$taxonomy)
  all <- ov_bind(tables)
  sens <- all[all$feature_id%in%selected$feature_id,,drop=FALSE]
  sens$main_beta <- tables$B1$beta_log2[match(sens$feature_id,tables$B1$feature_id)]
  sens$same_direction_as_main <- ifelse(sens$analysis_valid,
    sign(sens$beta_log2)==sign(sens$main_beta),NA)
  pairs <- saved[saved$feature_id%in%selected$feature_id,,drop=FALSE]
  pairs$model_label <- ifelse(pairs$model=='P2','Model 6','Model 5')
  # Keep original all-genus/all-comparison q values: never re-adjust the subset.
  key <- paste(saved$model,saved$feature_id,saved$comparison)
  stopifnot(!anyDuplicated(key),identical(pairs$P_BH,
    saved$P_BH[match(paste(pairs$model,pairs$feature_id,pairs$comparison),key)]))
  ps <- lapply(gp_pairs,function(k) {
    t <- pairs[pairs$model=='P2' & pairs$comparison==k,,drop=FALSE]
    data.frame(comparison=k,selected_genera=nrow(selected),eligible=sum(t$eligible),
      valid=sum(t$analysis_valid),positive_BH=sum(t$FDR_lt_005 & t$beta_log2_per_log2>0,na.rm=TRUE),
      negative_BH=sum(t$FDR_lt_005 & t$beta_log2_per_log2<0,na.rm=TRUE),
      original_BH_tests=unique(saved$FDR_tests[saved$model=='P2'])[1])
  })
  lr_check_hashes(inp$protected)
  list(protocol=sg_protocol,protected=inp$protected,tables=tables,all=all,
    samples=ov_bind(samples),selected=selected,sensitivity=sens,
    pairs=pairs,pair_summary=ov_bind(ps),all_pairs=saved,
    taxonomy=g$taxonomy,profiles=g$profiles,coverage=g$coverage,counts=g$counts,
    baseline_N=nrow(x$metadata),families=length(unique(x$metadata$family_id)))
}

sg_tables <- function(r,out) {
  v3_table(r$all,'All_personal_steatosis_genus_results_EN',out)
  for (nm in names(r$tables)) v3_table(r$tables[[nm]],paste0('Disease_',nm,'_all_genera_EN'),out)
  v3_table(r$samples,'Analysis_sample_and_discovery_summary_EN',out)
  v3_table(r$selected,'Selected_steatosis_associated_genera_EN',out)
  v3_table(r$sensitivity,'Selected_genera_sensitivity_results_EN',out)
  v3_table(r$pairs,'Selected_genera_cross_generation_results_EN',out)
  v3_table(r$all_pairs,'All_genus_cross_generation_results_original_BH_EN',out)
  v3_table(r$pair_summary,'Selected_genera_cross_generation_summary_EN',out)
  v3_table(r$profiles[r$profiles$feature_id%in%r$selected$feature_id,],
    'Selected_genera_generation_profiles_EN',out)
  v3_table(r$coverage,'Genus_annotation_coverage_EN',out)
  v3_table(r$counts,'Genus_aggregation_and_display_counts_EN',out)
  comparisons <- list(
    tf_compare(r$tables$Adult_B1_same_sample,r$tables$B2,'Adult Model 3 versus Model 4: identical participants'),
    tf_compare(r$tables$B1,r$tables$No_recent_antibiotics,'All-age Model 3 versus antibiotic exclusion: different samples'))
  comparison_summary <- ov_bind(lapply(comparisons,`[[`,'summary'))
  names(comparison_summary)[names(comparison_summary)=='common_valid_SGBs'] <- 'common_valid_genera'
  v3_table(comparison_summary,'Sensitivity_effect_comparison_summary_EN',out)
  v3_table(ov_bind(lapply(comparisons,`[[`,'data')),'Sensitivity_effect_comparisons_all_common_genera_EN',out)
  if (requireNamespace('openxlsx',quietly=TRUE)) {
    openxlsx::write.xlsx(list(Discovery_summary=r$samples,Selected_genera=r$selected,
      All_genus_disease_results=r$all,Selected_sensitivities=r$sensitivity,
      Selected_cross_generation=r$pairs,Cross_generation_summary=r$pair_summary),
      file.path(out,'tables/Steatosis_associated_genera_EN.xlsx'),overwrite=TRUE)
  }
}

sg_reports <- function(r,out) {
  n <- nrow(r$selected);s <- r$samples
  main <- s[s$model=='B1',];adult <- s[s$model=='B2',]
  lower <- r$selected$genus[r$selected$beta_log2<0]
  higher <- r$selected$genus[r$selected$beta_log2>0]
  lines <- c('# 脂肪肝相关菌属与代际特征',
    sprintf('基线%d人；全龄主模型%d人、%d例脂肪肝，检验%d个合格菌属/基因组属分组。',
      r$baseline_N,main$N,main$cases,main$eligible),
    sprintf('Model 3全范围BH q<0.05共%d个：降低%d个、升高%d个；全部展示，无30个上限。',
      n,main$lower,main$higher),
    sprintf('其中具名菌属%d个，未具名GGB基因组属分组%d个。',
      sum(!r$selected$unnamed_genome_genus_bin),sum(r$selected$unnamed_genome_genus_bin)),
    paste('降低：',if(length(lower))paste(lower,collapse=', ')else'无'),
    paste('升高：',if(length(higher))paste(higher,collapse=', ')else'无'),
    sprintf('成人Model 3/4使用同一批%d人、%d例；成人与全龄效应不能仅归因于增加协变量。',adult$N,adult$cases),
    paste0(s$model_label,'：BH支持',s$BH_supported,'个，降低',s$lower,'个、升高',s$higher,'个。'),
    '主模型调整代际、代际内年龄样条、性别、建库批次、测序池、采便季节、近期抗生素；家系随机截距。成人Model 4进一步加入吸烟、饮酒、饮食评分、体力活动、总能量和连续BMI。',
    '脂肪肝为超声判定；全龄保留未成年人。属由明确注释SGB逐人相加，未重新归一化。丰度只用于测量、检出率过滤和描述，不用于挑选前若干个属。',
    '主筛选固定为全龄Model 3。成人和抗生素排除结果只检验稳健性，不要求各模型均显著。BH分母包含全部合格属；无法可靠估计的合格模型按P=1计入分母。',
    paste0(r$pair_summary$comparison,'：所选属中Model 6可估计',r$pair_summary$valid,
      '个，正向BH支持',r$pair_summary$positive_BH,'个、负向',r$pair_summary$negative_BH,'个；沿用原',
      r$pair_summary$original_BH_tests,'项检验的BH。'),
    '同队列筛选后的跨代图为探索性整合展示，不是独立验证；支持丰度一致性，不证明菌株传播、因果或中介。原始数据和36号结果未更改。')
  writeLines(lines,file.path(out,'reports/Results_summary_zh.md'),useBytes=TRUE)
  writeLines(c('# Methods, references and study adaptations',
    'Ultrasound-defined hepatic steatosis is the exposure. The baseline contains primary families with at least two observed generations; minors remain in the all-age discovery model. Labels do not verify direct parent-child relationships or cohabitation.',
    'All SGBs with one unambiguous genus annotation were summed within each participant using the frozen module-36 aggregation. No closure to the classified fraction was applied. A genus must have >=10% prevalence and nonzero variance in the actual model sample. Named genera and unnamed GGB groups are explicit.',
    'Official MaAsLin2 1.24.1 linear mixed models use normalization NONE, transform LOG, no outcome standardization, and a primary-family random intercept. LOG replaces zeros by half the smallest positive abundance of the genus within that model sample. The steatosis contrast is Yes minus No. Model 3 includes generation, generation-specific age splines, sex, preparation batch, sequencing pool, seasonal sine/cosine and recent antibiotic use. Adult Model 4 adds smoking, alcohol, diet score, physical activity, log10 energy and continuous BMI. Adult Model 3 and Model 4 use identical participants. Recent-antibiotic exclusion retains the same basic adjustment with constant terms removed.',
    'Official coefficients and P values are checked against the stored lmerTest fit. Confidence intervals use Satterthwaite degrees of freedom. Nonconverged or mismatched fits are ineligible for discovery; converged singular random effects are retained and flagged. All eligible genera enter a separate BH family for each disease model. Invalid eligible fits count as P=1 but retain unavailable reported q. Per-batch package q values are discarded.',
    'Every valid all-age Model 3 genus with BH q<0.05 is displayed; there is no top-30 cap and no abundance ranking. Evidence ranks use q, raw P and ID. Taxonomy-map row ordering uses phylum, family and q. Multiple pages, when needed, retain every selected genus and common colour limits. Effect magnitude and evidence rank are distinct.',
    'Model 5/6 cross-generation coefficients and q values are reused unchanged from completed genus-level analyses, with no new regression fit and no BH recomputation on selected genera. Model 6 adjusts both endpoints for age splines, sex, batch, pool, seasonal sine/cosine and recent antibiotic use. Gaussian independence GEE has primary-family robust standard errors. Original BH includes all eligible genus-by-three-comparison tests. Selected genera below endpoint prevalence requirements are shown as unavailable.',
    'Disease effects are 2^beta for steatosis versus no steatosis on zero-replaced relative abundance. Pair effects are 2^beta for younger abundance per doubling of older abundance. They have different interpretations and independent legend scales. Neither is an OR, an absolute-load ratio, or a strain-transmission estimate. The classification hierarchy is not an inferred phylogenetic tree. In-sample selection/follow-up is exploratory and not independent validation.',
    '',
    'Mallick H, Rahnavard A, McIver LJ, et al. Multivariable association discovery in population-scale meta-omics studies. PLOS Computational Biology. 2021;17(11):e1009442. doi:10.1371/journal.pcbi.1009442. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009442 (Official multivariable and mixed-model framework; prevalence filtering. Family labels and the exact adjustment set are cohort adaptations.)',
    'Mei Z, Wang F, Bhosle A, et al. Strain-specific gut microbial signatures in type 2 diabetes identified in a cross-cohort analysis of 8,117 metagenomes. Nature Medicine. 2024;30:2265-2276. doi:10.1038/s41591-024-03067-7. https://www.nature.com/articles/s41591-024-03067-7 (Disease-association discovery and sensitivity to additional adjustment; our genus resolution, single-cohort design and q<0.05 threshold are explicit adaptations.)',
    'Valles-Colomer M, Bacigalupe R, Vieira-Silva S, et al. Variation and transmission of the human gut microbiota across multiple familial generations. Nature Microbiology. 2022;7:87-96. doi:10.1038/s41564-021-01021-8. https://www.nature.com/articles/s41564-021-01021-8 (Genus-level familial patterns motivate the question; the present older-to-younger abundance GEE is a cohort adaptation, not the source strain-transmission analysis.)'),
    file.path(out,'reports/Methods_and_references_EN.md'),useBytes=TRUE)
}

sg_validate <- function(r,out) {
  for(nm in names(r$tables)) {
    t <- r$tables[[nm]];eligible <- !grepl('^FILTERED',t$status)
    good <- eligible & tf_bool(t$analysis_valid) & is.finite(t$P_raw)
    p <- rep(1,sum(eligible));p[good[eligible]] <- t$P_raw[good]
    q <- p.adjust(p,'BH')
    stopifnot(max(abs(q[good[eligible]]-t$P_BH[good]))<1e-12,
      all(t$FDR_tests==sum(eligible)),
      max(abs(2^t$beta_log2[good]-t$abundance_ratio[good]))<1e-10)
  }
  main <- r$tables$B1
  stopifnot(setequal(r$selected$feature_id,main$feature_id[main$FDR_lt_005]),
    !anyDuplicated(r$selected$feature_id))
  a <- readRDS(file.path(out,'private/Adult_B1_same_sample_analysis_sample.rds'))
  b <- readRDS(file.path(out,'private/B2_analysis_sample.rds'))
  stopifnot(identical(a$sample_id,b$sample_id),identical(a$family_id,b$family_id))
  if(nrow(r$selected)) {
    paths <- list.files(file.path(out,'tables'),pattern='^Figure38_map_page_[0-9]+_source_EN[.]csv$',full.names=TRUE)
    cells <- ov_bind(lapply(paths,tf_read))
    stopifnot(setequal(cells$feature_id,r$selected$feature_id),all(table(cells$feature_id)==10L))
    paths <- list.files(file.path(out,'tables'),pattern='^Figure38_forest_page_[0-9]+_source_EN[.]csv$',full.names=TRUE)
    forest <- ov_bind(lapply(paths,tf_read))
    stopifnot(setequal(forest$feature_id,r$selected$feature_id),!anyDuplicated(forest$feature_id))
  }
  manifest <- tf_read(file.path(out,'tables/Figure_export_manifest_EN.csv'))
  for(stem in manifest$figure)for(ext in c('pdf','png','tiff','svg'))
    stopifnot(file.exists(file.path(out,'figures',paste0(stem,'.',ext))))
  input <- r$protected$path[basename(r$protected$path)=='genus_relative_abundance.rds']
  z <- readRDS(input);ids <- readRDS(file.path(out,'private/B1_analysis_sample.rds'))$sample_id
  at <- match(ids,z$metadata$sample_id);stopifnot(!anyNA(at))
  mass <- 100*colSums(z$relative)[at];disease <- z$metadata$fatty_liver[at]
  quality <- ov_bind(lapply(0:1,function(k){v<-mass[disease==k]
    data.frame(steatosis=k,N=length(v),mean_classified_percent=mean(v),
      median_classified_percent=median(v),Q1=unname(quantile(v,.25)),Q3=unname(quantile(v,.75)))
  }))
  fl_write_csv(quality,file.path(out,'audit/genus_annotation_coverage_by_steatosis.csv'))
  lr_check_hashes(r$protected)
  writeLines(c('PASS: independent full-family BH and ratio reconciliation.',
    'PASS: discovery set equals every valid primary-model q<0.05 genus.',
    'PASS: adult sample IDs and family IDs match exactly.',
    'PASS: every selected genus appears in ten map tracks and one forest row.',
    'PASS: all requested export formats exist; protected sources unchanged.'),
    file.path(out,'audit/numerical_and_display_validation.txt'))
}

sg_main <- function(args=character(),root=Sys.getenv('FATTYLIVER_ROOT',
  unset='E:/YFS/workspace/FattyLiver_generation_analysis')) {
  if (getRversion()<'4.5') {
    status <- system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',
      shQuote(file.path(root,'scripts/38_steatosis_selected_genera.R')),vapply(args,shQuote,character(1))))
    if (status!=0L) stop('Project-R steatosis-genus analysis failed')
    return(invisible(NULL))
  }
  modes <- c('--run','--check','--figures-only','--resume-figures','--help')
  if (any(!args%in%modes) || sum(args%in%modes)>1L) stop('Choose one execution mode')
  cfg <- sg_setup(root)
  sys.source(file.path(root,'R/steatosis_selected_genera_figures.R'),envir=.GlobalEnv)
  if ('--help'%in%args) {cat('38_steatosis_selected_genera.R [--run | --check | --figures-only | --resume-figures]\n');return(invisible(NULL))}
  if ('--check'%in%args) {sg_inputs(cfg);cat('PASS: syntax, project dependencies, ABI and completed inputs; no models run\n');return(invisible(NULL))}
  key <- 'steatosis_genus_generation';cfg$output_module <- unname(cfg$output_modules[key])
  cfg$run_id <- unname(cfg$run_ids[key]);cfg$run_id_explicit <- FALSE
  module <- file.path(root,'outputs',cfg$output_module)
  cache <- file.path(module,'private','checkpoints',sg_protocol)
  if ('--figures-only'%in%args) {
    out <- lr_adopted(root,cfg$output_module);r <- readRDS(file.path(out,'private/analysis_result.rds'))
    lr_check_hashes(r$protected)
    archive <- file.path(out,'audit',paste0('previous_figures_',format(Sys.time(),'%Y%m%d_%H%M%S')))
    dir.create(archive);files <- list.files(file.path(out,'figures'),full.names=TRUE)
    if(length(files))stopifnot(all(file.copy(files,archive)))
    sys.source(file.path(root,'R/steatosis_selected_genera_figures.R'),envir=.GlobalEnv)
    sg_preflight(cfg,out);sg_figures(r,cfg,out);sg_validate(r,out);lr_check_hashes(r$protected)
    fl_publish_organized_results(cfg);return(invisible(out))
  }
  resume <- NULL
  if ('--resume-figures'%in%args) {
    paths <- list.dirs(file.path(module,'history'),recursive=FALSE,full.names=TRUE)
    paths <- paths[file.exists(file.path(paths,'private/analysis_result.rds'))]
    if(!length(paths))stop('No saved completed inference available for figure recovery')
    scientific_audit <- file.exists(file.path(paths,'audit/B1_model_diagnostics.csv'))
    if(any(scientific_audit))paths <- paths[scientific_audit]
    resume <- paths[which.max(file.info(file.path(paths,'private/analysis_result.rds'))$mtime)]
  }
  ctx <- fl_output_start_run(cfg);cfg <- ctx$cfg;out <- ctx$run_dir
  for(p in c('tables','figures','reports','audit','private'))dir.create(file.path(out,p),recursive=TRUE)
  writeLines(c('RUNNING',sg_protocol),file.path(out,'run_status.txt'))
  tryCatch({
    if(is.null(resume))r <- sg_build(cfg,out,cache)else{
      r <- readRDS(file.path(resume,'private/analysis_result.rds'));lr_check_hashes(r$protected)
      for(p in c('private','audit')){
        dest <- if(p=='private')file.path(out,p)else file.path(out,p,'saved_fit_audit')
        dir.create(dest,showWarnings=FALSE)
        files <- list.files(file.path(resume,p),full.names=TRUE)
        # Keep deep development render trees in their original history location;
        # copying them repeatedly can exceed the Windows path-length limit.
        files <- files[!dir.exists(files)]
        if(length(files))stopifnot(all(file.copy(files,dest)))
      }
      writeLines('Saved inference recovered; no MaAsLin2/GEE refits.',file.path(out,'audit/recovery_mode.txt'))
    }
    r$counts$main_display_genera <- nrow(r$selected)
    r$counts$mean_total_RA_percent_in_display <- sum(r$selected$overall_mean_RA_percent)
    r$taxonomy$main_display <- r$taxonomy$feature_id%in%r$selected$feature_id
    saveRDS(r,file.path(out,'private/analysis_result.rds'));sg_tables(r,out);sg_reports(r,out)
    sys.source(file.path(root,'R/steatosis_selected_genera_figures.R'),envir=.GlobalEnv)
    sg_preflight(cfg,out);sg_figures(r,cfg,out);sg_validate(r,out);lr_check_hashes(r$protected)
    fl_write_csv(r$protected,file.path(out,'audit/protected_sources.csv'))
    writeLines(capture.output(sessionInfo()),file.path(out,'audit/sessionInfo.txt'))
    writeLines('PASS: protected inputs unchanged; all q<0.05 genera displayed; original pair BH retained; strict rendered QA passed.',
      file.path(out,'audit/delivery_validation.txt'))
    writeLines(c('COMPLETED',sg_protocol),file.path(out,'run_status.txt'))
    target <- fl_output_complete_run(ctx);writeLines(basename(target),file.path(module,'LATEST_RUN.txt'))
    tf_log('COMPLETED: ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
