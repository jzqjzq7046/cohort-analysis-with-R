# R-only structural adaptation of the adopted circular tree display to SGBs.
# Existing SGB discoveries define the tips; no genus effects are expanded to SGBs.
sst_protocol <- 'all_original_BH_supported_SGB_red_blue_ring_20261102_v4'
sst_stem <- 'Figure40_Steatosis_associated_SGB_circular_tree'

sst_setup <- function(root) {
  sys.source(file.path(root,'R/genus_main_text_figures.R'),envir=.GlobalEnv)
  sgm_setup(root)
}

sst_species_labels <- function(taxonomy,SGB_id) {
  rows <- lapply(seq_along(SGB_id),function(i) {
    tokens <- unlist(strsplit(taxonomy[i],'[|,]'))
    candidates <- sort(unique(trimws(sub('^s__','',tokens[grepl('^s__',tokens)]))))
    candidates <- candidates[!is.na(candidates)&nzchar(candidates)]
    readable <- gsub('_',' ',candidates,fixed=TRUE)
    named <- length(candidates)==1L &&
      grepl('^(Candidatus )?[A-Z][a-z]+ [a-z][a-z-]+$',readable) &&
      !grepl('\\b(bacterium|archaeon|unclassified|uncultured|unidentified|unknown|sp|spp)\\b',readable,ignore.case=TRUE,perl=TRUE)
    status <- if(length(candidates)>1L)'Conflicting species annotations'else
      if(length(candidates)==0L)'No species annotation'else
      if(named)'Unique named species annotation'else'Unnamed or incomplete species annotation'
    name <- if(named)readable else NA_character_
    data.frame(scientific_name=name,annotation_status=status,
      distinct_species_annotations=length(candidates),
      species_annotation_candidates=paste(readable,collapse='; '),
      figure_display_label=if(named)paste0(name,' (',SGB_id[i],')')else SGB_id[i],
      label_rule='Single distinct unambiguous binomial (optional Candidatus) -> italic name plus roman SGB; placeholder, missing or conflicting names -> SGB only',
      stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}

sst_geometry <- function(tree,taxonomy) {
  g <- gp_coordinates(tree,taxonomy,circular=TRUE,gap_degrees=12)
  # Preserve saved relative reference-tree branch lengths; rotate the fan only.
  scale <- .78/.96;angle <- pi/3
  rotate <- function(x,y)cbind(scale*(x*cos(angle)-y*sin(angle)),
    scale*(x*sin(angle)+y*cos(angle)))
  for(n in c('edges','arcs')) {
    xy <- rotate(g[[n]]$x,g[[n]]$y)
    g[[n]]$x <- xy[,1];g[[n]]$y <- xy[,2]
  }
  xy <- rotate(g$edges$xend,g$edges$yend)
  g$edges$xend <- xy[,1];g$edges$yend <- xy[,2]
  g$tips$theta <- g$tips$theta+angle
  g$tips$radius <- g$tips$radius*scale
  g$scale <- g$scale*scale
  g
}

sst_build <- function(cfg) {
  disease <- lr_adopted(cfg$project_root,'25_Taxa_robustness')
  reference <- lr_adopted(cfg$project_root,'35_Generation_reference_phylogeny')
  files <- c(file.path(disease,'tables/Disease_B1_all_SGBs_EN.csv'),
    file.path(reference,'private/figure_data.rds'),
    file.path(reference,'tables/All_SGB_annotations_and_overall_abundance_EN.csv'),
    file.path(cfg$project_root,'outputs/25_Taxa_robustness/LATEST_RUN.txt'),
    file.path(cfg$project_root,'outputs/35_Generation_reference_phylogeny/LATEST_RUN.txt'))
  protected <- lr_hashes(files)
  t <- tf_read(files[1]);ref <- readRDS(files[2])
  stopifnot(identical(ref$protocol,gp_protocol),!anyDuplicated(t$feature_id),
    !anyDuplicated(ref$taxonomy$feature_id),all(t$model=='B1'))
  lr_check_hashes(ref$protected)
  tested <- !grepl('^FILTERED',t$status)
  valid <- tested&tf_bool(t$analysis_valid)&is.finite(t$P_raw)&is.finite(t$P_BH)
  stopifnot(sum(tested)==unique(t$FDR_tests),all(valid[tested]))
  chosen <- valid&tf_bool(t$FDR_lt_005)&t$P_BH<.05
  s <- t[chosen,,drop=FALSE]
  stopifnot(nrow(s)>1,all(s$P_raw>0&s$P_raw<=1),all(s$abundance_ratio>0),
    setequal(t$feature_id[chosen],s$feature_id),all(s$feature_id%in%ref$tree$tip.label))
  at <- match(s$feature_id,ref$taxonomy$feature_id);stopifnot(!anyNA(at))
  for(k in c('phylum','kingdom','overall_mean_RA_percent','overall_prevalence_percent'))
    s[[k]] <- ref$taxonomy[[k]][at]
  # Recorded reference annotations contain legacy phylum synonyms. Preserve the
  # source string and harmonize only the display name, independently of disease.
  s$phylum_source <- s$phylum
  aliases <- c(Firmicutes='Bacillota',Bacteroidetes='Bacteroidota',
    Actinobacteria='Actinomycetota',Proteobacteria='Pseudomonadota',
    Verrucomicrobia='Verrucomicrobiota',Lentisphaerae='Lentisphaerota')
  legacy <- s$phylum%in%names(aliases)
  s$phylum[legacy] <- unname(aliases[s$phylum[legacy]])
  stopifnot(all(s$overall_mean_RA_percent>0),all(is.finite(s$overall_prevalence_percent)))
  s$model_label <- 'Model 3: all ages'
  s$selection_rule <- 'ALL valid SGB-specific Model 3 original full-family BH q<0.05; no count cap or genus-based selection'
  # Names are a label-only mapping; retain every original annotation unchanged.
  annotation <- sst_species_labels(s$taxonomy,s$SGB_id)
  for(k in names(annotation))s[[k]] <- annotation[[k]]
  s$display_name_source <- s$figure_display_label
  groups <- c('Bacillota','Bacteroidota','Actinomycetota','Pseudomonadota',
    'Verrucomicrobiota','Archaea','Other phyla')
  s$branch_colour_group <- s$phylum
  archaea <- grepl('Archaea',s$kingdom,ignore.case=TRUE)|s$phylum%in%c('Methanobacteriota','Thermoplasmatota')
  s$branch_colour_group[archaea] <- 'Archaea'
  s$branch_colour_group[!s$branch_colour_group%in%groups] <- 'Other phyla'
  colours <- setNames(c('#040676','#397FC7','#F1B656','#79CAFB','#F7A24F','#6A95A3','#B4B4B4'),groups)
  tree <- ape::ladderize(ape::keep.tip(ref$tree,s$feature_id),right=TRUE)
  probe <- sort(s$feature_id)[unique(round(seq(1,nrow(s),length.out=min(40,nrow(s)))))]
  distance_delta <- max(abs(gp_pair_distances(ref$tree,probe)-gp_pair_distances(tree,probe)))
  stopifnot(distance_delta<1e-9,length(tree$tip.label)==nrow(s),setequal(tree$tip.label,s$feature_id))
  g <- sst_geometry(tree,s);s <- s[match(g$tips$feature_id,s$feature_id),,drop=FALSE]
  s$leaf_order <- as.integer(g$tips$leaf_order);s$tip_index <- sprintf('%03d',s$leaf_order)
  s$theta <- g$tips$theta;s$step <- g$tips$step;s$radius <- g$tips$radius
  s$label_shown <- TRUE;s$evidence_bar <- -log10(s$P_raw)
  ra <- range(log10(s$overall_mean_RA_percent));ra <- c(floor(ra[1]),ceiling(ra[2]))
  if(diff(ra)==0)ra <- ra+c(-1,1)
  effect_limit <- max(1,ceiling(max(abs(s$beta_log2))))
  limits <- list(RA=ra,disease=c(-effect_limit,effect_limit))
  s$prevalence_fill_fraction <- s$overall_prevalence_percent/100
  s$RA_fill_fraction <- (log10(s$overall_mean_RA_percent)-ra[1])/diff(ra)
  stopifnot(all(s$prevalence_fill_fraction>=0&s$prevalence_fill_fraction<=1),
    all(s$RA_fill_fraction>=0&s$RA_fill_fraction<=1))
  counts <- data.frame(candidates=nrow(t),tested=sum(tested),filtered=sum(!tested),
    significant=nrow(s),lower=sum(s$beta_log2<0),higher=sum(s$beta_log2>0),
    uniquely_named=sum(!is.na(s$scientific_name)),SGB_only=sum(is.na(s$scientific_name)),
    N=unique(s$N),cases=unique(s$cases),baseline_N=ref$baseline_N,families=ref$families)
  phylum_counts <- data.frame(branch_colour_group=groups,
    SGBs=as.integer(table(factor(s$branch_colour_group,levels=groups))))
  phylum_counts$percentage_of_significant_SGBs <- 100*phylum_counts$SGBs/nrow(s)
  r <- list(protocol=sst_protocol,taxonomy=s,all_results=t,tree=tree,geometry=g,colours=colours,
    limits=limits,counts=counts,phylum_counts=phylum_counts,N=ref$baseline_N,
    association_colors=c(Reference='#2166AC',Comparison='#B2182B'),
    families=ref$families,model_N=unique(s$N),cases=unique(s$cases),
    protected=protected,inherited_protected=ref$protected,distance_delta=distance_delta,
    reference_baseline_tips=length(ref$tree$tip.label),
    evidence_bar_max=ceiling(max(s$evidence_bar)),tip_name_column='SGB_id',
    tip_scientific_name_column='scientific_name',
    use_actual_tip_radius=TRUE,canvas_width_mm=590,canvas_height_mm=448,
    plot_x_limits=c(-1.75,2.859375),plot_y_limits=c(-1.75,1.75),
    legend_font_multiplier=2.75,legend_key_left=1.92,legend_key_right=2.72,
    legend_swatch_half_height=.018)
  # All 379 names need a larger reading master; do not shrink below the text floor.
  name_positions <- sgt_name_positions(r)
  stopifnot(nrow(name_positions)==nrow(s),!anyDuplicated(name_positions$display_name),
    identical(name_positions$display_name,s$figure_display_label))
  r
}

sst_tables <- function(r,out,cfg) {
  v3_table(r$counts,'SGB_tree_analysis_summary_EN',out)
  v3_table(r$phylum_counts,'Significant_SGB_phylum_counts_and_percentages_EN',out)
  v3_table(r$all_results,'All_original_SGB_Model3_results_EN',out)
  s <- r$taxonomy[order(r$taxonomy$leaf_order),,drop=FALSE]
  v3_table(s,'Figure40_All_significant_SGB_tip_key_EN',out)
  v3_table(s[,c('feature_id','SGB_id','scientific_name','annotation_status',
    'distinct_species_annotations','species_annotation_candidates','figure_display_label','taxonomy')],
    'Figure40_Species_annotation_and_display_labels_EN',out)
  statuses <- as.data.frame(table(s$annotation_status),stringsAsFactors=FALSE)
  names(statuses) <- c('annotation_status','SGBs')
  v3_table(statuses,'Figure40_Annotation_status_counts_EN',out)
  v3_table(sgt_cells(r,cfg),'Figure40_Three_ring_source_EN',out)
  label_positions <- sgt_name_positions(r)
  v3_table(label_positions,'Figure40_Full_SGB_radial_label_positions_EN',out)
  v3_table(sgt_mixed_text_parts(label_positions,r$canvas_width_mm,diff(r$plot_x_limits)),
    'Figure40_Italic_species_and_roman_SGB_text_parts_EN',out)
  v3_table(s[,c('feature_id','SGB_id','P_raw','P_BH','evidence_bar')],
    'Figure40_Outer_bar_original_P_evidence_EN',out)
  v3_table(data.frame(scale=c('Mean RA log10 percent lower','Mean RA log10 percent upper',
    'Steatosis log2 ratio lower','Steatosis log2 ratio upper'),value=c(r$limits$RA,r$limits$disease)),
    'Figure40_Colour_limits_EN',out)
  v3_table(data.frame(ring=3L,position=c('Lower','Centre','Higher'),
    abundance_ratio=2^c(r$limits$disease[1],0,r$limits$disease[2]),
    colour=c(unname(r$association_colors['Reference']),'#FFFFFF',unname(r$association_colors['Comparison'])),
    meaning=c('Lower abundance with steatosis','Abundance ratio = 1','Higher abundance with steatosis')),
    'Figure40_Association_ring_palette_EN',out)
  ape::write.tree(r$tree,file=file.path(out,'tables/Significant_SGB_pruned_reference_tree.nwk'))
  fl_write_csv(r$protected,file.path(out,'audit/protected_sources.csv'))
  v3_table(data.frame(item=c('Reference baseline tips','Retained significant tips',
    'Reference-distance audit maximum difference','New phylogenetic inference',
    'Model refits','Multiple-testing recalculation'),
    value=c(r$reference_baseline_tips,length(r$tree$tip.label),format(r$distance_delta,scientific=TRUE),
      'None; existing RDS reference tree pruned only','None','None')),
    'Reference_tree_and_inference_audit_EN',out)
}

sst_reports <- function(r,out) {
  z <- r$counts
  writeLines(c('# Figure 40 | Reference-tree distribution of steatosis-associated SGBs.',
    sprintf('All %d species-level genome bins (SGBs) meeting the original full-family BH q<0.05 in the saved all-age Model 3 are shown; %d have lower and %d have higher abundance with ultrasound-defined hepatic steatosis. Selection is among %d tested SGBs, not by abundance, constituent genera or a fixed display count. Disease models included %s participants, including %s cases.',
      z$significant,z$lower,z$higher,z$tested,format(z$N,big.mark=','),format(z$cases,big.mark=',')),
    'The pre-existing reference tree stored with the original microbiome RDS is pruned to exactly the significant SGBs. Original reference-tree relative branch lengths are retained; branch colours denote phylum, with Archaea and other phyla shown as broad display groups. Legacy phylum synonyms are harmonized for display (for example, Firmicutes to Bacillota); the original phylum string remains in the source table. The exact source database version and inference method are not recorded in the available RDS. This is not a newly inferred study-specific maximum-likelihood tree. An SGB is not a sample-specific strain.',
    sprintf('All %d SGBs have radial labels and grey dotted terminal guides. A unique unambiguous species annotation is shown as an italic scientific name followed by its roman SGB identifier in parentheses (%d SGBs). Unnamed, incomplete or conflicting species annotations are shown by SGB identifier alone (%d SGBs); all original annotations and label decisions are supplied in the corresponding table. A unique annotation requires a single distinct binomial, optionally preceded by Candidatus; placeholder names and multiple species candidates are not resolved to an invented name. Orange-brown labels identify higher abundance with steatosis; dark labels identify lower abundance. From inside out, rings show prevalence (blue, linear percentage), mean relative abundance (orange, common logarithmic positive-mean intensity) and the adjusted steatosis-versus-no-steatosis abundance ratio (blue-white-red on a symmetric log2 scale centred at 1). Descriptive tracks use all %s baseline microbiome participants, include non-detection zeros and retain the original full-profile denominator; selected SGBs are not renormalized. Outer teal bars show -log10(original raw P), not genotype counts.',
      z$significant,z$uniquely_named,z$SGB_only,format(z$baseline_N,big.mark=',')),
    sprintf('Model 3 uses official MaAsLin2 log2 zero-replaced relative abundance with generation, generation-specific age splines, sex, preparation batch, sequencing pool, sampling season and recent antibiotics, with a primary-family random intercept. Zeros use the official half-smallest-positive rule within the original model sample. Ratios refer to this transformed relative-abundance scale, not absolute microbial load or disease odds. All original beta, CI, raw P and BH q values are reused unchanged from the complete %d-test family. No selected-subset BH calculation or model fitting is performed.',z$tested),
    sprintf('The %d x %d mm full-label reading master retains a >=5-pt SGB text floor; it is not a 183-mm journal-production figure. The enlarged right-side legend has 14.6-16.5-pt physical text, 2.75 times its preceding size, with wider colour keys and an expanded right-side layout. Its larger typography is an explicit reading-master contract, not a changed journal-production requirement. Prevalence, relative abundance, association colours, outer bars and label geometry are supplied as source tables. No pTR, genotype count or strain-transmission result is inferred.',r$canvas_width_mm,r$canvas_height_mm),
    'Visual composition: Valles-Colomer M, Bacigalupe R, Vieira-Silva S, et al. Variation and transmission of the human gut microbiota across multiple familial generations. Nature Microbiology. 2022;7:87-96. doi:10.1038/s41564-021-01021-8. Figure 2 supplies the layout; available YFS data determine the ring/bar contents.'),
    file.path(out,'reports/Figure40_legend_EN.md'),useBytes=TRUE)
  writeLines(c('# 显著SGB环形参考树',
    sprintf('复用25号SGB脂肪肝Model 3结果：%d人，%d例；%d个候选中%d个通过原模型过滤并受检，%d个原BH q<0.05，其中降低%d个、升高%d个。',
      z$N,z$cases,z$candidates,z$tested,z$significant,z$lower,z$higher),
    '每个显著SGB恰好一个末端。不是将显著属下面的SGB全部展开，也不画全3949个SGB；不按丰度或前30名取舍。',
    sprintf('标签按用户新规则：%d个SGB具有唯一明确的物种注释，显示斜体物种名及正体括号SGB编号；其余%d个仅显示SGB编号。单独表格保留完整原始注释、不同物种候选、标签及判定原因。不同上级分类但物种名一致的重复注释合并；物种名冲突不擅自选取。',z$uniquely_named,z$SGB_only),
    sprintf('环形布局沿用上一版。第三环及其图例改为蓝-白-红：丰度比低于1为蓝，高于1为红，1为白，原对称log2范围不变。检出率环、平均相对丰度环、树枝、名称和外柱配色均保留。外柱为-log10(原始P)。前两环分母为%d人，模型分母为%d人。',z$baseline_N,z$N),
    '参考树从35号保存对象裁剪，保留原树相对分支长度。确切参考库版本和建树算法在当前RDS中未记录，不标为本研究新建最大似然树；SGB也不等同于菌株。',
    '门名称仅为显示统一同义名称，例如Firmicutes显示为Bacillota；原始门名称保留在phylum_source列，分类与配色不按疾病结果改变。',
    '2026-11-01图例比例修订：图例文字放大2.75倍（约14.6–16.5pt），加宽色标和右侧区域。SGB名称仍为5.1pt；不改变任何树枝、数值、颜色或选择集合。',
    sprintf('%d个完整SGB标识以%d×%d毫米阅读版输出，PDF、PNG、TIFF和SVG。未重拟合、未重算BH；前一版菌属图及其他分析保留。',z$significant,r$canvas_width_mm,r$canvas_height_mm)),
    file.path(out,'reports/Results_summary_zh.md'),useBytes=TRUE)
}

sst_preflight <- function(cfg,out) {
  fun <- c('sst_build','sst_species_labels','sst_geometry','sst_figures','sgt_name_positions',
    'sgt_mixed_text_parts','sgt_plot','sgt_cells',
    'sgt_ring_polygons','gp_coordinates','gp_tint','sg_ramp','tf_export','nf_theme')
  path <- file.path(out,'audit/actual_plot_and_export_source.R')
  writeLines(unlist(lapply(fun,function(n)c(paste0(n,' <- '),deparse(get(n),width.cutoff=120),''))),path)
  invisible(parse(path,encoding='UTF-8'))
  nf_qa_command(cfg,'validate_figure.py',c(shQuote(path),'--json'),file.path(out,'audit/plot_source_preflight.json'))
}

sst_figures <- function(r,cfg,out) {
  p <- sgt_plot(r,cfg)
  manifest <- tf_export(list(plot=patchwork::wrap_plots(list(p)),width=r$canvas_width_mm,
    height=r$canvas_height_mm,panels=1L,exemptions=list(),fixed_aspect_single_panel=TRUE,
    no_panel_tags=TRUE,oriented_text=TRUE,max_body_pt=16.51,
    large_format_font_reason='User-requested larger right-side legend in the 590x448-mm reading master; SGB labels remain 5.1pt, legend text 14.6-16.5pt, not a journal-production file'),sst_stem,out,cfg)
  v3_table(manifest,'Figure_export_manifest_EN',out)
}

sst_validate <- function(r,out,cfg) {
  t <- tf_read(r$protected$path[1]);s <- r$taxonomy
  eligible <- !grepl('^FILTERED',t$status)&tf_bool(t$analysis_valid)&tf_bool(t$FDR_lt_005)
  stopifnot(setequal(t$feature_id[eligible],s$feature_id),setequal(r$tree$tip.label,s$feature_id),
    !anyDuplicated(s$feature_id),nrow(s)==length(r$tree$tip.label),r$distance_delta<1e-9)
  columns <- c('beta_log2','SE','CI_low','CI_high','P_raw','P_BH','abundance_ratio','ratio_CI_low','ratio_CI_high')
  stopifnot(identical(s[,columns],t[match(s$feature_id,t$feature_id),columns]))
  cells <- tf_read(file.path(out,'tables/Figure40_Three_ring_source_EN.csv'))
  stopifnot(setequal(cells$feature_id,s$feature_id),all(table(cells$feature_id)==3L))
  annotations <- sst_species_labels(s$taxonomy,s$SGB_id)
  stopifnot(identical(annotations$figure_display_label,s$figure_display_label),
    all(s$figure_display_label[is.na(s$scientific_name)]==s$SGB_id[is.na(s$scientific_name)]))
  example <- match('SGB15339',s$SGB_id)
  if(!is.na(example))stopifnot(s$figure_display_label[example]=='Faecalibacterium prausnitzii (SGB15339)')
  conflict_example <- match('SGB15322',s$SGB_id)
  if(!is.na(conflict_example))stopifnot(s$figure_display_label[conflict_example]=='SGB15322')
  for(ext in c('pdf','png','tiff','svg'))stopifnot(file.exists(file.path(out,'figures',paste0(sst_stem,'.',ext))))
  lr_check_hashes(r$protected);lr_check_hashes(r$inherited_protected)
  writeLines(c('PASS: every original primary BH-supported SGB is shown exactly once; no abundance selection.',
    'PASS: original beta, SE, CI, raw P, q and the complete testing family are unchanged.',
    'PASS: pruned reference-tree path distances are preserved; no new tree inference.',
    'PASS: three observed data cells and one complete SGB label per tip.',
    'PASS: unique named species use italic name plus roman SGB; unresolved labels use SGB only; complete original annotations retained.',
    'PASS: four export formats; strict rendered text/alignment/collision checks passed.'),
    file.path(out,'audit/delivery_validation.txt'))
}

sst_main <- function(args=character(),root=Sys.getenv('FATTYLIVER_ROOT',
  unset='E:/YFS/workspace/FattyLiver_generation_analysis')) {
  if(getRversion()<'4.5') {
    status <- system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',
      shQuote(file.path(root,'scripts/40_steatosis_SGB_circular_tree.R')),vapply(args,shQuote,character(1))))
    if(status!=0L)stop('Project-R SGB circular-tree export failed')
    return(invisible(NULL))
  }
  if(length(args)>1L||any(!args%in%c('--run','--check','--help')))stop('Choose --run, --check or --help')
  cfg <- sst_setup(root)
  if('--help'%in%args){cat('40_steatosis_SGB_circular_tree.R [--run | --check]\nSaved SGB results and reference tree only; no model fits.\n');return(invisible(NULL))}
  if('--check'%in%args){r<-sst_build(cfg);print(r$counts);cat('PASS: original selection, reference mapping, branch distances and full-label data contract. No fits or figure exports.\n');return(invisible(r$counts))}
  cfg$output_module <- unname(cfg$output_modules['steatosis_sgb_tree'])
  cfg$run_id <- unname(cfg$run_ids['steatosis_sgb_tree']);cfg$run_id_explicit <- FALSE
  ctx <- fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(p in c('tables','figures','reports','audit','private'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
  writeLines(c('RUNNING',sst_protocol),file.path(out,'run_status.txt'))
  tryCatch({
    r <- sst_build(cfg);saveRDS(r,file.path(out,'private/display_result.rds'))
    print(r$counts);sst_tables(r,out,cfg);sst_reports(r,out)
    captions <- file.path(root,'analysis_plan/figure_captions_20261102_EN.md')
    if(file.exists(captions)&&!file.copy(captions,file.path(out,'reports/Updated_figure_captions_EN.md'),overwrite=TRUE))
      stop('Cannot preserve the updated manuscript figure captions')
    sst_preflight(cfg,out)
    sst_figures(r,cfg,out);sst_validate(r,out,cfg)
    writeLines(capture.output(sessionInfo()),file.path(out,'audit/sessionInfo.txt'))
    writeLines(c('COMPLETED',sst_protocol),file.path(out,'run_status.txt'))
    target <- fl_output_complete_run(ctx)
    writeLines(basename(target),file.path(root,'outputs',cfg$output_module,'LATEST_RUN.txt'))
    tf_log('COMPLETED: ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
