# Inert definitions. Published display structure; cohort-specific saved effects.
# Complete discovery and supplementary maps remain in module 38.
sgm_protocol <- 'all_tests_volcano_265_named_genera_reference_figure_style_20261031_v4'

sgm_setup <- function(root) {
  sys.source(file.path(root,'R/steatosis_selected_genera.R'),envir=.GlobalEnv)
  cfg <- sg_setup(root)
  sys.source(file.path(root,'R/steatosis_selected_genera_figures.R'),envir=.GlobalEnv)
  sys.source(file.path(root,'R/significant_genus_circular_taxonomy.R'),envir=.GlobalEnv)
  cfg
}

sgm_build <- function(cfg,out) {
  parent <- lr_adopted(cfg$project_root,'38_Steatosis_selected_genus_generation')
  files <- c(file.path(parent,'private/analysis_result.rds'),
    file.path(parent,'tables/Disease_B1_all_genera_EN.csv'),
    file.path(cfg$project_root,'outputs/38_Steatosis_selected_genus_generation/LATEST_RUN.txt'))
  protected <- lr_hashes(files)
  r <- readRDS(files[1]);stopifnot(identical(r$protocol,sg_protocol))
  lr_check_hashes(r$protected)
  t <- r$tables$B1;stopifnot(!anyDuplicated(t$feature_id))
  csv <- tf_read(files[2]);stopifnot(isTRUE(all.equal(t,csv,tolerance=1e-12,check.attributes=FALSE)))
  eligible <- !grepl('^FILTERED',t$status)
  z <- t[eligible,,drop=FALSE]
  # Current data have a valid estimate for every eligible genus. Do not silently
  # drop a future non-estimable taxon from an all-tested volcano.
  stopifnot(nrow(z)==unique(z$FDR_tests),all(tf_bool(z$analysis_valid)),
    all(is.finite(z$P_raw)&z$P_raw>0&z$P_raw<=1),all(is.finite(z$P_BH)),
    all(is.finite(z$abundance_ratio)&z$abundance_ratio>0))
  sig <- z[z$FDR_lt_005,,drop=FALSE]
  sig <- sig[order(sig$P_BH,sig$P_raw,sig$feature_id),,drop=FALSE]
  chosen <- sig$feature_id
  stopifnot(setequal(sig$feature_id,r$selected$feature_id))
  s <- r$selected[match(chosen,r$selected$feature_id),,drop=FALSE]
  s$main_evidence_rank <- seq_len(nrow(s))
  s$main_selection_rule <- rep('ALL valid primary all-age Model 3 BH q<0.05; no count cap; q/P/ID evidence rank',nrow(s))
  s <- s[order(s$phylum,s$family,s$P_BH,s$P_raw,s$feature_id),,drop=FALSE]
  s$display_row <- seq_len(nrow(s))
  z$main_text_detail <- z$feature_id%in%chosen
  z$minus_log10_raw_P <- -log10(z$P_raw)
  z$plot_class <- ifelse(!z$FDR_lt_005,'Not BH-supported',
    ifelse(z$beta_log2<0,'Lower with steatosis','Higher with steatosis'))
  # This line is a plotting translation of the ORIGINAL full-family BH result.
  # No new p.adjust operation is performed on selected taxa.
  raw_boundary <- if(nrow(sig)) .05*nrow(sig)/unique(z$FDR_tests) else NA_real_
  if(is.finite(raw_boundary))stopifnot(identical(z$FDR_lt_005,z$P_raw<raw_boundary))
  counts <- data.frame(N=z$N[1],cases=z$cases[1],tested=nrow(z),BH_supported=nrow(sig),
    lower=sum(z$FDR_lt_005 & z$beta_log2<0),higher=sum(z$FDR_lt_005 & z$beta_log2>0),
    not_BH_supported=sum(!z$FDR_lt_005),main_detail=nrow(s),
    main_detail_lower=sum(s$beta_log2<0),main_detail_higher=sum(s$beta_log2>0),
    main_detail_named=sum(!s$unnamed_genome_genus_bin),
    main_detail_GGB=sum(s$unnamed_genome_genus_bin),
    raw_P_at_BH_boundary=raw_boundary)
  # Retain the supplement's colour limits so main/supplement cells are comparable.
  limits <- sg_scales(r)
  details <- r;details$selected <- s
  list(protocol=sgm_protocol,parent=parent,protected=protected,
    inherited_protected=r$protected,volcano=z,selected=s,counts=counts,
    limits=limits,detail_source=details,all_supported=sig,
    circular_taxonomy=sgt_build(s,r$baseline_N,r$families,limits,cfg))
}

sgm_tables <- function(a,out,cfg) {
  v3_table(a$counts,'Main_text_display_summary_EN',out)
  v3_table(a$volcano,'Figure39_All_tested_genus_volcano_source_EN',out)
  v3_table(a$selected,'Figure39_All_significant_genus_selection_EN',out)
  v3_table(a$all_supported,'All_primary_BH_supported_genera_EN',out)
  v3_table(a$detail_source$sensitivity[a$detail_source$sensitivity$feature_id%in%a$selected$feature_id,],
    'Main_text_genera_adult_and_antibiotic_sensitivities_EN',out)
  v3_table(a$detail_source$pairs[a$detail_source$pairs$feature_id%in%a$selected$feature_id,],
    'Main_text_genera_cross_generation_original_BH_EN',out)
  sgt_tables(a$circular_taxonomy,out,cfg)
}

sgm_volcano <- function(a,cfg) {
  z <- a$volcano;s <- a$counts
  levels <- c('Lower with steatosis','Higher with steatosis','Not BH-supported')
  z$plot_class <- factor(z$plot_class,levels=levels)
  colours <- c('Lower with steatosis'=unname(cfg$binary_colors['Reference']),
    'Higher with steatosis'=unname(cfg$binary_colors['Comparison']),'Not BH-supported'='#C7C7C7')
  labels <- c(sprintf('Lower (%d)',s$lower),sprintf('Higher (%d)',s$higher),
    sprintf('Not significant (%d)',s$not_BH_supported))
  p <- ggplot2::ggplot(z,ggplot2::aes(abundance_ratio,minus_log10_raw_P))+
    ggplot2::geom_vline(xintercept=1,linetype=2,colour='#888888',linewidth=.25)+
    ggplot2::geom_point(ggplot2::aes(fill=plot_class),shape=21,colour='#555555',
      size=1.6,stroke=.12,alpha=.85)+
    ggplot2::scale_fill_manual(values=colours,breaks=levels,labels=labels,drop=FALSE)+
    ggplot2::scale_x_log10(breaks=c(.0625,.125,.25,.5,1,2,4,8),
      labels=function(v)format(v,trim=TRUE,scientific=FALSE,digits=3),
      expand=ggplot2::expansion(mult=c(.07,.07)))+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.015,.06)))+
    ggplot2::labs(title=NULL,x='Adjusted abundance ratio (log scale)',
      y='-log10(raw P)',fill=NULL)+tf_plot_theme()
  if(is.finite(s$raw_P_at_BH_boundary))p <- p+
    ggplot2::geom_hline(yintercept=-log10(s$raw_P_at_BH_boundary),linetype=3,
      colour='#888888',linewidth=.25)
  p
}

sgm_preflight <- function(cfg,out) {
  fun <- c('sgm_volcano','sgm_figures','sg_map','sg_cells','sg_key','sg_scales','sg_ramp',
    'sg_forest','gg_geometry','gp_tint','gp_ra_fraction','gp_text','gp_rect',
    'sgt_plot','sgt_fan_geometry','sgt_name_positions','sgt_cells','sgt_ring_polygons','gp_coordinates','gp_circular_theme','lr_export','tf_export','nf_theme')
  path <- file.path(out,'audit/actual_plot_and_export_source.R')
  writeLines(unlist(lapply(fun,function(n)c(paste0(n,' <- '),deparse(get(n),width.cutoff=120),''))),path)
  invisible(parse(path,encoding='UTF-8'))
  nf_qa_command(cfg,'validate_figure.py',c(shQuote(path),'--json'),file.path(out,'audit/plot_source_preflight.json'))
}

sgm_reuse_unchanged <- function(a,reuse,out) {
  old <- readRDS(file.path(reuse,'private/display_result.rds'))
  stopifnot(identical(a$volcano,old$volcano),identical(a$selected,old$selected))
  manifest <- tf_read(file.path(reuse,'tables/Figure_export_manifest_EN.csv'))
  manifest <- manifest[!grepl('Cohort_reference_phylogeny|Significant_genus_circular_taxonomy',manifest$figure),,drop=FALSE]
  stopifnot(nrow(manifest)==1L+2L*ceiling(nrow(a$selected)/40L))
  copied <- character()
  for(stem in manifest$figure) {
    for(ext in c('pdf','png','tiff','svg')) {
      f <- file.path(reuse,'figures',paste0(stem,'.',ext));dest <- file.path(out,'figures',basename(f))
      stopifnot(file.copy(f,dest),identical(unname(tools::md5sum(f)),unname(tools::md5sum(dest))))
      copied <- c(copied,f)
    }
    audit <- list.files(file.path(reuse,'audit'),pattern=paste0('^',stem,'[.]'),full.names=TRUE)
    stopifnot(all(file.copy(audit,file.path(out,'audit'))))
  }
  data <- list.files(file.path(reuse,'tables'),pattern='^Figure39_(map_page_|family_display_labels_page_|forest_page_)',full.names=TRUE)
  stopifnot(all(file.copy(data,file.path(out,'tables'))))
  fl_write_csv(lr_hashes(copied),file.path(out,'audit/reused_unchanged_figures_sources.csv'))
  manifest
}

sgm_figures <- function(a,cfg,out,reuse=NULL) {
  manifest <- list()
  if(!is.null(reuse))manifest <- list(sgm_reuse_unchanged(a,reuse,out))
  if(is.null(reuse)) {
  manifest[[1]] <- lr_export(list(sgm_volcano(a,cfg)),
    'Figure39_Main_All_genus_volcano',out,cfg,width=183,height=112)
  if(nrow(a$selected)) {
    # Forty is a physical per-page row limit; ALL discoveries are represented.
    groups <- split(seq_len(nrow(a$selected)),ceiling(seq_len(nrow(a$selected))/40L))
    seen <- character()
    for(k in seq_along(groups)) {
      z <- a$selected[groups[[k]],,drop=FALSE];seen <- c(seen,z$feature_id)
      p <- sg_map(a$detail_source,z,cfg,a$limits,k,out,source_prefix='Figure39')
      manifest[[length(manifest)+1L]] <- lr_export(list(p),
        sprintf('Figure39_All_significant_genus_heatmap_page_%02d',k),out,cfg,
        width=183,height=max(120,(4*nrow(z)+82)*.94))
    }
    stopifnot(identical(seen,a$selected$feature_id),!anyDuplicated(seen))
    evidence <- a$selected[order(a$selected$main_evidence_rank),,drop=FALSE]
    groups <- split(seq_len(nrow(evidence)),ceiling(seq_len(nrow(evidence))/40L))
    for(k in seq_along(groups)) {
      z <- evidence[groups[[k]],,drop=FALSE]
      v3_table(z,sprintf('Figure39_forest_page_%02d_source_EN',k),out)
      p <- sg_forest(z,cfg,'Steatosis-associated genera: Model 3')
      manifest[[length(manifest)+1L]] <- lr_export(list(p),
        sprintf('Figure39_All_significant_genus_forest_page_%02d',k),out,cfg,
        width=120,height=max(85,35+4.1*nrow(z)))
    }
  }
  }
  p <- sgt_plot(a$circular_taxonomy,cfg)
  manifest[[length(manifest)+1L]] <- tf_export(list(plot=patchwork::wrap_plots(list(p)),
    width=360,height=288,panels=1L,exemptions=list(),fixed_aspect_single_panel=TRUE,no_panel_tags=TRUE,oriented_text=TRUE),
    'Figure39_Significant_genus_circular_taxonomy',out,cfg)
  v3_table(ov_bind(manifest),'Figure_export_manifest_EN',out)
  v3_table(data.frame(scale=c('Mean RA log10 percent lower','Mean RA log10 percent upper',
    'Steatosis log2 ratio lower','Steatosis log2 ratio upper','Pair log2 ratio lower','Pair log2 ratio upper'),
    value=c(a$limits$RA,a$limits$disease,a$limits$pair)),
    'Figure39_shared_main_and_supplement_colour_limits_EN',out)
}

sgm_reports <- function(a,out) {
  n <- a$counts
  writeLines(c('# 全部受检分组总览、全部显著分组及显著属环形分类树',
    sprintf('仅复用38号保存结果：%d人、%d例脂肪肝，全部%d个合格分组进入火山图；其中BH q<0.05共%d个，降低%d个、升高%d个。',
      n$N,n$cases,n$tested,n$BH_supported,n$lower,n$higher),
    sprintf('详细图保留全部%d个显著分组：降低%d个、升高%d个，具名菌属%d个、GGB分组%d个。不再选择前30个；每页最多40行为物理分页规则，不是筛选上限。森林图按q、原始P、编号排布；排名体现统计证据，不代表生物学效应大小。',
      n$main_detail,n$main_detail_lower,n$main_detail_higher,n$main_detail_named,n$main_detail_GGB),
    '火山图横轴为调整后丰度比（对数刻度），纵轴为-log10(原始P)；蓝色为脂肪肝相关降低，橙色为升高，灰色为未达到全范围BH q<0.05。已取消前30个分组的额外黑色外圈。',
    sprintf('水平虚线为原%d项BH检验的原始P边界：0.05 × %d / %d = %.6g；不是把原始P<0.05当作FDR显著。',
      n$tested,n$BH_supported,n$tested,n$raw_P_at_BH_boundary),
    '热图保留分类层级、个人脂肪肝关联、三代描述特征及Model 6跨代关联。森林图展示全部显著分组的个人脂肪肝丰度比及95%CI。两图集合完全相同，热图按分类排布，森林图按证据排序。',
    '没有重新拟合MaAsLin2或GEE，没有重新计算BH，也没有按相对丰度挑选。颜色沿用现有方案，主图与完整补充图使用同一热图色标。',
    sprintf('环形分类树仅展示%d个显著属/GGB，每属恰好一个末端。补全原始注释已有的界、纲、目层级，按界—门—纲—目—科—属连接；不代表序列推断的进化距离，不展开为SGB。全部265个完整名称沿圆周放射标注，向外接三个环带与证据柱形。',nrow(a$circular_taxonomy$taxonomy)),
    sprintf('按提供的参考图排列：虚线引导、放射全名、黑色环带边界、右侧紧凑图例。三个外环从内至外为检出率、平均相对丰度、脂肪肝调整丰度比。前两者按基线%d人计算，包含0值、未重新归一化。外柱为已有模型的-log10(原始P)，不是基因型数；橙棕名称标出15个丰度升高的属。完整命名图以360×288mm导出，保持所有文字至少5pt；不是183mm的双栏定稿尺寸。分类节点用后代属数量安排绘图半径，不表示进化距离或时间。',a$circular_taxonomy$N),
    '不构建最大似然树，不计算pTR或基因型数。图式借鉴Valles-Colomer等Figure 2的环形布局；分类层级定量展示参考Foster等Metacoder框架，实际使用自定义R绘图。',
    '上一版全SGB参考树由输出管理器移入历史版本，并从当前39号展示目录移除；原38号完整结果保留。'),
    file.path(out,'reports/Main_text_display_summary_zh.md'),useBytes=TRUE)
  writeLines(c('# Figure 39 | All steatosis-associated genera and their detailed profiles.',
    sprintf('All %d eligible genus/genome-genus groups tested in %d participants (%d with ultrasound-defined hepatic steatosis) are shown in the volcano plot. Colours indicate full-family BH q<0.05: lower (blue, n=%d) or higher (orange, n=%d) abundance with steatosis; grey indicates q>=0.05 (n=%d). The horizontal axis shows the adjusted abundance ratio on a logarithmic scale; the vertical axis shows -log10(raw P). No selected-subset outer rings are drawn. The horizontal dotted line translates the original full-family BH boundary into the raw-P scale (alpha times the number of discoveries divided by the full number of tests); no new multiple-testing correction was performed.',
      n$tested,n$N,n$cases,n$lower,n$higher,n$not_BH_supported),
    sprintf('The detailed classification heatmap and companion forest plot show ALL %d valid primary all-age Model 3 BH discoveries. There is no top-30 count limit, abundance ranking or balanced allocation by direction. Taxonomic row order is used in the heatmaps; q, raw P and identifier determine forest order. At most 40 rows per page preserve final-size legibility; pagination does not change the displayed set.',n$main_detail),
    'Heatmap tracks show phylum-family-genus classification (not phylogeny), personal steatosis abundance ratio, generation-specific detection and mean relative abundance, and Model 6 cross-generation abundance ratios. Disease and pair ratios have separate scales centred at 1. Colour limits match the complete supplementary maps. * indicates the original full-family BH q<0.05; grey denotes not tested or not estimable. Full names for abbreviated provisional families are retained in source tables.',
    'Points and horizontal bars in the forest plot represent adjusted abundance ratios and 95% confidence intervals from the original official MaAsLin2 log2 relative-abundance mixed models. Model 3 adjusts for generation, generation-specific age splines, sex, preparation batch, sequencing pool, sampling season and recent antibiotics, with a primary-family random intercept. Confidence intervals use Satterthwaite degrees of freedom. Ratios refer to zero-replaced relative abundance, not absolute microbial load or disease odds.',
    'Model 6 adjusts both pair endpoints for age splines, sex, batch, pool, season and recent antibiotics using primary-family-clustered Gaussian GEE. Its original all-genus/all-comparison q values remain unchanged. These same-cohort cross-generation follow-ups are exploratory and do not establish strain transmission or independent validation. All original tests remain in the discovery result module.'),
    file.path(out,'reports/Figure39_legends_EN.md'),useBytes=TRUE)
  writeLines(c('# Literature basis and display adaptations',
    'Fackelmann G, Manghi P, Carlino N, et al. Gut microbiome signatures of vegan, vegetarian and omnivore diets and associated health outcomes across 21,561 individuals. Nature Microbiology. 2025;10:41-52. doi:10.1038/s41564-024-01870-z. Figures 3-4 display the 30 SGBs with the largest absolute SMD, with complete results in supplementary tables. https://www.nature.com/articles/s41564-024-01870-z',
    'Asnicar F, Manghi P, Fackelmann G, et al. Gut micro-organisms associated with health, nutrition and dietary interventions. Nature. 2026;650:450-458. doi:10.1038/s41586-025-09854-7. Figure 2 selects the 15 top and bottom health-ranked SGBs for visualization. https://www.nature.com/articles/s41586-025-09854-7',
    'Mei Z, Wang F, Bhosle A, et al. Strain-specific gut microbial signatures in type 2 diabetes identified in a cross-cohort analysis of 8,117 metagenomes. Nature Medicine. 2024;30:2265-2276. doi:10.1038/s41591-024-03067-7. Figure 2 combines an overview of significant taxa with selected detailed profiles. https://www.nature.com/articles/s41591-024-03067-7',
    'The current user-approved revision displays ALL primary BH discoveries, without a fixed count cap. The preceding top-30 and all-SGB tree versions are archived. The cited selected-detail examples provide context, not a requirement to limit the number displayed. The new circular taxonomic hierarchy has exactly one terminal node for each of the 265 significant genus/GGB groups. Its ring composition is visually adapted from Valles-Colomer et al. Nature Microbiology 2022 Figure 2; it is not a sequence phylogeny. No disease/family models or multiple-testing families are changed.'),
    file.path(out,'reports/Display_literature_and_adaptations_EN.md'),useBytes=TRUE)
  sgt_legend(a$circular_taxonomy,out)
}

sgm_validate <- function(a,out) {
  z <- a$volcano
  stopifnot(!anyDuplicated(z$feature_id),nrow(z)==a$counts$tested,
    sum(z$FDR_lt_005)==a$counts$BH_supported,
    identical(a$selected$feature_id[order(a$selected$main_evidence_rank)],
      a$all_supported$feature_id))
  cell_files <- list.files(file.path(out,'tables'),pattern='^Figure39_map_page_[0-9]+_source_EN[.]csv$',full.names=TRUE)
  cells <- ov_bind(lapply(cell_files,tf_read))
  stopifnot(setequal(cells$feature_id,a$selected$feature_id),all(table(cells$feature_id)==10L))
  forest_files <- list.files(file.path(out,'tables'),pattern='^Figure39_forest_page_[0-9]+_source_EN[.]csv$',full.names=TRUE)
  forests <- ov_bind(lapply(forest_files,tf_read))
  stopifnot(!anyDuplicated(forests$feature_id),setequal(forests$feature_id,a$selected$feature_id))
  r <- readRDS(a$protected$path[1]);main <- r$tables$B1
  cols <- c('beta_log2','SE','CI_low','CI_high','P_raw','P_BH','abundance_ratio','ratio_CI_low','ratio_CI_high')
  stopifnot(identical(z[,cols],main[match(z$feature_id,main$feature_id),cols]))
  manifest <- tf_read(file.path(out,'tables/Figure_export_manifest_EN.csv'))
  stopifnot(nrow(manifest)==2L+2L*ceiling(nrow(a$selected)/40L))
  for(stem in manifest$figure)for(ext in c('pdf','png','tiff','svg'))
    stopifnot(file.exists(file.path(out,'figures',paste0(stem,'.',ext))))
  lr_check_hashes(a$protected);lr_check_hashes(a$inherited_protected)
  tr <- a$circular_taxonomy
  lr_check_hashes(tr$protected)
  stopifnot(!anyDuplicated(tr$taxonomy$feature_id),setequal(tr$taxonomy$feature_id,a$selected$feature_id),
    identical(tr$tree$tip.label,tr$taxonomy$feature_id),length(tr$tree$tip.label)==nrow(a$selected))
  tree_cells <- tf_read(file.path(out,'tables/Figure39_Significant_genus_circular_tree_cells_EN.csv'))
  stopifnot(setequal(tree_cells$feature_id,a$selected$feature_id),all(table(tree_cells$feature_id)==3L),
    !any(grepl('Cohort_reference_phylogeny',manifest$figure)))
  writeLines(c('PASS: every eligible tested group is represented in the volcano.',
    'PASS: detail set equals ALL primary BH-supported groups; no count cap.',
    'PASS: all heatmap/forest pages cover identical complete IDs once; ten heatmap tracks each.',
    'PASS: circular taxonomy contains exactly all significant genus IDs, one tip and three observed data cells each.',
    'PASS: all-SGB tree removed from current bundle; no inferred genetic distances, pTR or genotype proxies.',
    'PASS: coefficients, confidence intervals, raw P and original q are unchanged.',
    'PASS: all four export formats and strict rendered checks completed.'),
    file.path(out,'audit/delivery_validation.txt'))
}

sgm_main <- function(args=character(),root=Sys.getenv('FATTYLIVER_ROOT',
  unset='E:/YFS/workspace/FattyLiver_generation_analysis')) {
  if(getRversion()<'4.5') {
    status <- system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',
      shQuote(file.path(root,'scripts/39_genus_main_text_figures.R')),vapply(args,shQuote,character(1))))
    if(status!=0L)stop('Project-R main-text genus visualization failed')
    return(invisible(NULL))
  }
  if(any(!args%in%c('--run','--check','--help','--taxonomy-only')) || length(args)>1L)stop('Choose --run, --check, --help or --taxonomy-only')
  cfg <- sgm_setup(root)
  if('--help'%in%args){cat('39_genus_main_text_figures.R [--run | --check | --taxonomy-only]\nSaved results only; never refits models.\n');return(invisible(NULL))}
  if('--check'%in%args){lr_adopted(root,'38_Steatosis_selected_genus_generation');cat('PASS: saved parent, syntax, dependencies and ABI; no statistical fitting.\n');return(invisible(NULL))}
  reuse <- if('--taxonomy-only'%in%args)lr_adopted(root,'39_Genus_main_text_figures')else NULL
  cfg$output_module <- unname(cfg$output_modules['genus_main_text'])
  cfg$run_id <- unname(cfg$run_ids['genus_main_text']);cfg$run_id_explicit <- FALSE
  ctx <- fl_output_start_run(cfg);cfg <- ctx$cfg;out <- ctx$run_dir
  for(p in c('tables','figures','reports','audit','private'))dir.create(file.path(out,p),recursive=TRUE)
  writeLines(c('RUNNING',sgm_protocol),file.path(out,'run_status.txt'))
  tryCatch({
    a <- sgm_build(cfg,out);saveRDS(a,file.path(out,'private/display_result.rds'))
    sgm_tables(a,out,cfg);sgm_reports(a,out);sgm_preflight(cfg,out);sgm_figures(a,cfg,out,reuse)
    sgm_validate(a,out);fl_write_csv(a$protected,file.path(out,'audit/protected_sources.csv'))
    writeLines(capture.output(sessionInfo()),file.path(out,'audit/sessionInfo.txt'))
    writeLines(c('COMPLETED',sgm_protocol),file.path(out,'run_status.txt'))
    target <- fl_output_complete_run(ctx)
    writeLines(basename(target),file.path(root,'outputs',cfg$output_module,'LATEST_RUN.txt'))
    tf_log('COMPLETED: ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
