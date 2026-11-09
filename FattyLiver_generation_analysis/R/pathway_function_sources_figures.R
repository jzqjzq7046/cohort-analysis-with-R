# R-only, saved-result displays. Figure selection never changes inference.
pfs_theme <- function() ade_theme()+ggplot2::theme(
  plot.title=ggplot2::element_text(size=6.7,margin=ggplot2::margin(b=6)),
  axis.title=ggplot2::element_text(size=6.3),
  legend.title=ggplot2::element_text(size=6.1,margin=ggplot2::margin(b=5)),
  legend.text=ggplot2::element_text(size=5.6),
  legend.margin=ggplot2::margin(4,4,4,4),legend.box.spacing=grid::unit(3,'mm'),
  legend.position='bottom')
pfs_series <- function(r,ids) {
  specs<-data.frame(series=c('All ages','Adults: core','Adults: full','G0','G1','G2'),
    model=c('B1','Adult_B1_same_sample','B2','B1','B1','B1'),
    generation=c('Overall','Overall','Overall','G0','G1','G2'),stringsAsFactors=FALSE)
  rows<-lapply(seq_len(nrow(specs)),function(j){
    q<-specs[j,];z<-r$effects[r$effects$model==q$model&r$effects$generation==q$generation,,drop=FALSE]
    k<-match(ids,z$feature_id);z<-z[k,,drop=FALSE];z$feature_id<-ids;z$series<-q$series;z
  })
  z<-ov_bind(rows);z$series<-factor(z$series,levels=specs$series)
  z$feature_id<-factor(z$feature_id,levels=rev(ids))
  z$fill<-ifelse(tf_bool(z$analysis_valid),z$beta_log2,NA_real_)
  z$star<-ifelse(tf_bool(z$FDR_lt_005),'*','');z
}
pfs_effect_limit <- function(z) {
  v<-z$fill[is.finite(z$fill)];stopifnot(length(v)>0)
  max(.25,ceiling(max(abs(v))*4)/4)
}
pfs_effect_plot <- function(z,lim,cfg,title=NULL,labels=ggplot2::waiver()) {
  ggplot2::ggplot(z,ggplot2::aes(series,feature_id,fill=fill))+
    ggplot2::geom_tile(width=.96,height=.92,colour='white',linewidth=.15)+
    ggplot2::geom_text(ggplot2::aes(label=star),size=5.8/(72.27/25.4),colour='black')+
    ggplot2::scale_fill_gradient2(low=unname(cfg$binary_colors[1]),mid='white',high=unname(cfg$binary_colors[2]),
      midpoint=0,limits=c(-lim,lim),breaks=c(-lim,0,lim),labels=ade_ratio_text(2^c(-lim,0,lim)),
      name='Adjusted abundance ratio',na.value='#E5E5E5',
      guide=ggplot2::guide_colourbar(direction='horizontal',title.position='top',display='rectangles',
        nbin=48,ticks=FALSE,barwidth=grid::unit(45,'mm'),barheight=grid::unit(2.8,'mm')))+
    ggplot2::scale_y_discrete(drop=FALSE,labels=labels,expand=ggplot2::expansion(add=.55))+
    ggplot2::scale_x_discrete(drop=FALSE,expand=ggplot2::expansion(add=.55))+
    ggplot2::labs(title=title,x=NULL,y=NULL)+pfs_theme()+
    ggplot2::theme(axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank())
}
pfs_tagged <- function(plots,ncol=1,heights=NULL) {
  patchwork::wrap_plots(plots,ncol=ncol,heights=heights)+patchwork::plot_annotation(tag_levels='a',
    theme=ggplot2::theme(plot.tag=ggplot2::element_text(size=8,face='bold',family='Arial')))
}
pfs_export <- function(p,stem,out,cfg,width=183,height=180,panels=1L,uniform_y_axes=TRUE) {
  tf_export(list(plot=p,width=width,height=height,panels=panels,exemptions=list(),uniform_y_axes=uniform_y_axes),stem,out,cfg)
}
pfs_legacy_main_figure <- function(r,cfg,out,preview_only=FALSE) {
  z<-pfs_primary(r);z<-z[tf_bool(z$analysis_valid),,drop=FALSE]
  z$state<-factor(ifelse(!tf_bool(z$FDR_lt_005),'Not significant',ifelse(z$beta_log2<0,'Lower','Higher')),
    levels=c('Lower','Higher','Not significant'))
  z$y<--log10(pmax(z$P_raw,.Machine$double.xmin))
  # Literal shared biological ratio axis; failure is explicit instead of clipping.
  axis<-list(limits=c(.5,2),breaks=c(.5,.75,1,1.5,2))
  stopifnot(all(z$abundance_ratio>=axis$limits[1]&z$abundance_ratio<=axis$limits[2]))
  bhline<-if(any(z$FDR_lt_005))max(z$P_raw[z$FDR_lt_005])else NA_real_
  p1<-ggplot2::ggplot(z,ggplot2::aes(abundance_ratio,y,colour=state))+
    ggplot2::geom_vline(xintercept=1,linetype=2,colour='#888888',linewidth=.25)+
    ggplot2::geom_point(size=1.25,alpha=.85,shape=16)+
    ggplot2::scale_colour_manual(values=c(Lower=unname(cfg$binary_colors[1]),Higher=unname(cfg$binary_colors[2]),
      'Not significant'='#BDBDBD'),drop=FALSE)+
    ggplot2::scale_x_log10(limits=axis$limits,breaks=axis$breaks,labels=c('0.5','0.75','1','1.5','2'))+
    ggplot2::labs(title='Overall pathway associations',x='Adjusted abundance ratio (log scale)',
      y='−log10(raw P)',colour=NULL)+pfs_theme()+
    # Equal physical A-B and B-C plot-area gutters despite the different
    # discrete-point and colourbar legend heights (checked on final device).
    ggplot2::theme(plot.margin=ggplot2::margin(6,7,15.15,5,'pt'))
  if(is.finite(bhline))p1<-p1+ggplot2::geom_hline(yintercept=-log10(bhline),linetype=3,colour='#999999',linewidth=.25)
  ids<-r$representatives$feature_id
  eff<-pfs_series(r,ids);lim<-pfs_effect_limit(eff)
  labelmap<-setNames(paste0(gsub('-',' ',r$representatives$display_parent,fixed=TRUE),'\n',ids),ids)
  p2<-pfs_effect_plot(eff,lim,cfg,'Representative functional categories',labels=labelmap)
  src<-r$sources$display
  genera<-c(r$sources$top_genera,'Other named genera','Unclassified')
  src<-src[src$group%in%c('No_steatosis','Steatosis'),,drop=FALSE]
  src$feature_id<-factor(src$feature_id,levels=rev(ids))
  ng<-length(genera)
  src$x<-match(src$display_genus,genera)+ifelse(src$group=='Steatosis',ng+1,0)
  ticks<-c(seq_len(ng),seq_len(ng)+ng+1)
  # One physical panel with two explicitly separated source blocks. This avoids
  # an artificial extra inferred endpoint or a misleading species-level effect.
  p3<-ggplot2::ggplot(src,ggplot2::aes(x,feature_id,fill=source_share_percent))+
    ggplot2::geom_tile(width=.94,height=.92,colour='white',linewidth=.15)+
    ggplot2::scale_fill_gradient(low='white',high='#397FC7',limits=c(0,100),breaks=c(0,50,100),
      name='Share of summed stratified abundance (%)',
      guide=ggplot2::guide_colourbar(direction='horizontal',title.position='top',display='rectangles',
        nbin=48,ticks=FALSE,barwidth=grid::unit(45,'mm'),barheight=grid::unit(2.8,'mm')))+
    ggplot2::scale_x_continuous(breaks=ticks,labels=rep(genera,2),expand=ggplot2::expansion(add=.55))+
    ggplot2::scale_y_discrete(drop=FALSE,expand=ggplot2::expansion(add=c(.55,1.5)))+
    ggplot2::annotate('text',x=c((ng+1)/2,(ng+1)/2+ng+1),y=length(ids)+1,
      label=c('No steatosis','Steatosis'),size=6/(72.27/25.4),family='Arial')+
    ggplot2::labs(title='Observed genus-stratified sources',x=NULL,y=NULL)+pfs_theme()+
    ggplot2::theme(axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=.5,size=5.4))
  p<-pfs_tagged(list(p1,p2,p3),heights=c(1.1,1,1.2))
  if(preview_only)return(p)
  ade_table(z,'Figure44a_all_pathway_points_EN',out)
  ade_table(eff,'Figure44b_representative_adjusted_effects_EN',out)
  ade_table(src,'Figure44c_observed_source_matrix_EN',out)
  ade_table(data.frame(panel='a',x_lower=.5,x_upper=2,ticks='0.5;0.75;1;1.5;2',
    BH_boundary_raw_P=bhline,effect_colour_log2_limit=lim), 'Figure44_axis_and_colour_contract_EN',out)
  pfs_export(p,'Figure44_Pathway_function_and_sources',out,cfg,height=268,panels=3L)
}
pfs_main_figure <- function(r,cfg,out,preview_only=FALSE) {
  if(is.null(r$adult81_display))r<-pfs_prepare_adult81(r,cfg,out)
  # Reuse the exact original overview plot: all 420 all-age tests, original
  # point colours, values, BH boundary and fixed literal ratio ticks.
  a<-pfs_legacy_main_figure(r,cfg,out,preview_only=TRUE)[[1]]+
    ggplot2::labs(title=NULL)+ggplot2::theme(plot.title=ggplot2::element_blank(),
      plot.margin=ggplot2::margin(2,2,4,2,'pt'),
      legend.margin=ggplot2::margin(1,1,1,1,'pt'),legend.box.spacing=grid::unit(1,'mm'))
  heat<-pfs_adult81_panels(r,cfg)
  # Six effect columns and twenty source columns intentionally have different
  # grid spans. One pathway-label axis serves both adjacent, row-aligned maps.
  p<-patchwork::wrap_plots(a,heat$b,heat$c,design='AAAAA\nBBCCC',heights=c(1,4.0))+
    patchwork::plot_annotation(tag_levels='a',
      theme=ggplot2::theme(plot.tag=ggplot2::element_text(size=8,face='bold',family='Arial')))
  if(preview_only)return(p)
  z<-pfs_primary(r);z<-z[tf_bool(z$analysis_valid),,drop=FALSE]
  z$state<-factor(ifelse(!tf_bool(z$FDR_lt_005),'Not significant',ifelse(z$beta_log2<0,'Lower','Higher')),
    levels=c('Lower','Higher','Not significant'));z$y<--log10(pmax(z$P_raw,.Machine$double.xmin))
  ade_table(z,'Figure44a_all_pathway_points_EN',out)
  ade_table(heat$source,'Figure44c_observed_source_matrix_EN',out)
  ade_table(data.frame(panel='a',x_lower=.5,x_upper=2,ticks='0.5;0.75;1;1.5;2',
    BH_boundary_raw_P=max(z$P_raw[z$FDR_lt_005]),effect_colour_log2_limit=heat$colour_limit,
    display_rule='B/C: all original adult Model-4 BH discoveries, matched functional row order'),
    'Figure44_axis_and_colour_contract_EN',out)
  pfs_export(p,'Figure44_Pathway_function_and_sources_compact',out,cfg,width=210,height=310,panels=3L,
    uniform_y_axes=FALSE)
}
pfs_complete_maps <- function(r,cfg,out,ids,prefix) {
  if(!length(ids))return(data.frame())
  full<-pfs_series(r,ids);lim<-pfs_effect_limit(full)
  pages<-split(ids,ceiling(seq_along(ids)/32L));exports<-rows<-list()
  for(j in seq_along(pages)){
    keep<-pages[[j]];z<-full[as.character(full$feature_id)%in%keep,,drop=FALSE]
    z$feature_id<-factor(as.character(z$feature_id),levels=rev(keep))
    p<-pfs_effect_plot(z,lim,cfg,NULL)
    p<-pfs_tagged(list(p));stem<-sprintf('%s_%02d',prefix,j)
    exports[[j]]<-pfs_export(p,stem,out,cfg,width=183,height=45+length(keep)*3.7)
    rows[[j]]<-data.frame(figure=stem,page=j,display_order=seq_along(keep)+(j-1)*32L,
      feature_id=keep,colour_limit_log2=lim,selection='All discoveries in the original complete BH family; no count cap')
  }
  ade_table(ov_bind(rows),paste0(prefix,'_complete_page_manifest_EN'),out)
  ade_table(full,paste0(prefix,'_complete_matrix_EN'),out)
  ov_bind(exports)
}
pfs_figures <- function(r,cfg,out) {
  # Saved-result render recovery is exact: both plotting inputs/functions and
  # all vector/raster files plus their strict QA records must match their hashes.
  signature<-digest::digest(list(effects=r$effects,heterogeneity=r$heterogeneity,
    adult81=r$adult81_display,representatives=r$representatives,
    support_md5=unname(tools::md5sum(file.path(cfg$project_root,'R/pathway_adult81_heatmaps.R'))),
    binary=cfg$binary_colors,generation=cfg$generation_colors,
    functions=lapply(list(pfs_main_figure,pfs_legacy_main_figure,pfs_series,pfs_effect_plot,pfs_effect_limit,
      pfs_complete_maps,pfs_theme,pfs_tagged,pfs_export,tf_export),body)),algo='sha256')
  cachepath<-file.path(out,'audit/pfs_figure_cache.rds')
  if(file.exists(cachepath)){
    cache<-readRDS(cachepath);paths<-file.path(out,cache$files$relative_path)
    if(identical(signature,cache$signature)&&all(file.exists(paths))&&
       identical(unname(tools::md5sum(paths)),cache$files$md5)){
      ade_log('Reused byte-verified four-format figures and strict QA; plotting inputs/functions unchanged')
      return(invisible(NULL))
    }
  }
  # Unchanged supplementary figures retain their original four-format bytes
  # and strict QA. Verify data equality and every reused file before reuse.
  reuse<-NULL
  if(file.exists(cachepath)&&!is.null(r$recovered_from)){
    old<-readRDS(file.path(r$recovered_from,'private/analysis_result.rds'))
    cache<-readRDS(cachepath)
    keep<-grepl('(^|/)FigureS44_All_',cache$files$relative_path)
    paths<-file.path(out,cache$files$relative_path[keep])
    if(any(keep)&&identical(old$effects,r$effects)&&identical(old$heterogeneity,r$heterogeneity)&&
       all(file.exists(paths))&&identical(unname(tools::md5sum(paths)),cache$files$md5[keep])){
      previous<-tf_read(file.path(out,'tables/Figure_export_manifest_EN.csv'))
      reuse<-previous[grepl('^FigureS44_All_',previous$figure),,drop=FALSE]
      ade_log('Preserving ',nrow(reuse),' unaffected complete-discovery/interaction maps byte-for-byte')
    }
  }
  # Retire only the copied, superseded B/C-only display files in this NEW run.
  # Previous adopted versions remain in history via the output manager.
  retired<-c(list.files(file.path(out,'figures'),pattern='^Figure(44bc_All81_|S44_Adult81_)',full.names=TRUE),
    list.files(file.path(out,'audit'),pattern='^Figure(44bc_All81_|S44_Adult81_)',full.names=TRUE),
    list.files(file.path(out,'figures'),pattern='^Figure44_Pathway_function_and_sources\\.',full.names=TRUE),
    list.files(file.path(out,'audit'),pattern='^Figure44_Pathway_function_and_sources\\.',full.names=TRUE),
    list.files(file.path(out,'tables'),pattern='^FigureS44_Adult81_complete_page_manifest_EN\\.',full.names=TRUE))
  if(length(retired)){
    root<-normalizePath(out,winslash='/',mustWork=TRUE)
    resolved<-normalizePath(retired,winslash='/',mustWork=TRUE)
    stopifnot(all(startsWith(resolved,paste0(root,'/'))))
    unlink(retired,recursive=FALSE)
  }
  exports<-list(pfs_main_figure(r,cfg,out))
  primary<-pfs_primary(r);p<-primary[tf_bool(primary$FDR_lt_005),,drop=FALSE]
  a<-r$classification$annotation;p$display_parent<-a$display_parent[match(p$feature_id,a$feature_id)]
  p<-p[order(p$display_parent,p$P_BH,p$P_raw,p$feature_id),,drop=FALSE]
  h<-r$heterogeneity[r$heterogeneity$model=='B1'&tf_bool(r$heterogeneity$FDR_lt_005),,drop=FALSE]
  h<-h[order(h$P_BH,h$P_raw,h$feature_id),,drop=FALSE]
  if(is.null(reuse)){
    exports[[2]]<-pfs_complete_maps(r,cfg,out,p$feature_id,'FigureS44_All_primary_pathway_discoveries')
    exports[[3]]<-pfs_complete_maps(r,cfg,out,h$feature_id,'FigureS44_All_generation_interactions')
  }else exports[[2]]<-reuse
  ade_table(ov_bind(exports),'Figure_export_manifest_EN',out)
  paths<-c(list.files(file.path(out,'figures'),pattern='^Figure.*\\.(pdf|png|tiff|svg)$',full.names=TRUE),
    list.files(file.path(out,'audit'),pattern='^Figure.*\\.(json|csv)$',full.names=TRUE),
    file.path(out,'tables/Figure_export_manifest_EN.csv'))
  relative<-substring(gsub('\\\\','/',paths),nchar(gsub('\\\\','/',out))+2L)
  saveRDS(list(signature=signature,files=data.frame(relative_path=relative,
    md5=unname(tools::md5sum(paths)))),cachepath)
  invisible(NULL)
}
pfs_reports <- function(r,cfg,out) {
  p<-pfs_primary(r);sig<-p[tf_bool(p$FDR_lt_005),,drop=FALSE]
  adult<-r$effects[r$effects$model=='B2'&r$effects$generation=='Overall',,drop=FALSE]
  h<-r$heterogeneity[r$heterogeneity$model=='B1',,drop=FALSE]
  cls<-r$classification$annotation
  diagnostics<-ov_bind(lapply(names(r$fits),function(key){
    z<-r$fits[[key]]$diagnostics;z$fit_set<-key;z$mode<-r$fits[[key]]$mode;z
  }))
  ade_table(diagnostics,'All_model_diagnostics_EN',out)
  diagnostic_summary<-ov_bind(lapply(split(diagnostics,diagnostics$fit_set),function(z)
    data.frame(fit_set=z$fit_set[1],mode=z$mode[1],features=nrow(z),
      estimated=sum(z$status=='ESTIMATED'),singular=sum(z$status=='ESTIMATED_SINGULAR'),
      unreliable=sum(!z$status%in%c('ESTIMATED','ESTIMATED_SINGULAR')))))
  ade_table(diagnostic_summary,'Convergence_and_singularity_summary_EN',out)
  summary<-data.frame(metric=c('Primary participants','Primary steatosis cases','Supplied community pathways',
    'Primary tested pathways','Primary BH discoveries','Lower with steatosis','Higher with steatosis',
    'Adult fully adjusted BH discoveries','Primary discoveries also significant in adult full model',
    'Joint generation-interaction BH discoveries','Historical lineage mapping: all pathways',
    'Historical lineage mapping: primary discoveries','Main displayed adult Model-4 pathways',
    'Actual stratified rows','Distinct source taxon labels','Source taxa uniquely matched to named SGBs',
    'Existing primary SGB discoveries matched exactly'),
    value=c(unique(p$N),unique(p$cases),nrow(r$annotation),nrow(p),nrow(sig),sum(sig$beta_log2<0),sum(sig$beta_log2>0),
      sum(tf_bool(adult$FDR_lt_005)),length(intersect(sig$feature_id,adult$feature_id[tf_bool(adult$FDR_lt_005)])),
      sum(tf_bool(h$FDR_lt_005)),sum(cls$mapped),sum(cls$mapped[match(sig$feature_id,cls$feature_id)]),
      nrow(r$adult81_display$selection),r$sources$source_rows,r$sources$taxon_strata,
      r$source_taxonomy_context$unique_matches,r$source_taxonomy_context$matched_primary_SGBs))
  ade_table(summary,'Pathway_redesign_result_summary_EN',out)
  ade_table(sig[order(sig$P_BH,sig$P_raw,sig$feature_id),,drop=FALSE],'All_primary_pathway_discoveries_EN',out)
  ade_table(h[tf_bool(h$FDR_lt_005),,drop=FALSE],'All_generation_interaction_discoveries_EN',out)
  ade_searchable_table(r,out)
  comp<-r$composition$effects
  lines<-c('# Pathway function and observed microbial sources',
    '',sprintf('Primary sample: %s participants, including %s with ultrasound-defined hepatic steatosis. MASLD is not analysed.',unique(p$N),unique(p$cases)),
    sprintf('Of %d eligible supplied community pathways, %d passed BH q < 0.05 (%d lower, %d higher).',nrow(p),nrow(sig),sum(sig$beta_log2<0),sum(sig$beta_log2>0)),
    sprintf('Adult fully adjusted discoveries: %d; overlapping the all-age discoveries: %d. Joint generation-interaction discoveries: %d.',
      sum(tf_bool(adult$FDR_lt_005)),length(intersect(sig$feature_id,adult$feature_id[tf_bool(adult$FDR_lt_005)])),sum(tf_bool(h$FDR_lt_005))),
    '', '## Functional community composition',
    vapply(seq_len(nrow(comp)),function(i)sprintf('%s; %s: N=%d, R2 of total corrected variation=%.6f%%, partial R2=%.6f%%, pseudo-F=%.5f, permutation P=%.4g (%d draws).',
      ade_model_label(comp$model[i]),comp$scheme[i],comp$N[i],100*comp$R2_total[i],100*comp$R2_partial[i],comp$pseudo_F[i],comp$P_raw[i],comp$permutations[i]),character(1)),
    'These are family-conditioned dbRDA associations. They do not estimate the same between/within-family mixture as the random-intercept feature models. The all-age within-family test is primary; five additional tests are sensitivity checks. The separate six-test BH column is provided in the table.',
    '', '## Main figure',
    'Panel a shows every valid primary pathway test. Colour identifies BH-supported lower/higher abundance; grey denotes q>=0.05. The horizontal dotted line is the largest raw P among the BH discoveries, not a nominal P=0.05 line.',
    'Panels b/c display ALL 81 original valid adult Model-4 overall BH q<0.05 pathways. No additional top-count cap, subset BH or model fitting is used. Rows are organised by a verified historical MetaCyc major class, subclass and deeper lineage; within the deepest class, complete-linkage clustering uses the original available six log2 coefficients without row standardisation or zero imputation. Both panels use the identical row order and group separators. Five unmatched pathways remain in a separate final block.',
    'A canonical lineage is chosen only for display using a fixed major-class priority and lexical remaining lineage. All true multiple parents are retained in the selection table. Functional classification and clustering are descriptive layout choices, not category enrichment or new significance tests.',
    'All ages, adult core and adult full columns are Model 3, Model 3 and Model 4 respectively. G0/G1/G2 use generation-specific disease coefficients from a single Model-3 model. Stars use the relevant original complete effect BH family; a generation-specific star is not a joint interaction test.',
    'Pathway selection is based on the 2,040-person adult Model-4 analysis. The G0/G1/G2 coefficients remain the original all-age Model-3 estimates. Panel c retains the original 3,193-person source sample and is not an adult Model-4-adjusted source comparison.',
    'Panel c uses actual HUMAnN species strata aggregated to their reported genera. Each source share is the group mean stratified quantity divided by the sum of all group mean strata for that pathway. It is descriptive, unadjusted and not a share of the community pathway total. The top eight named genera are chosen using pooled source shares across all 81 selected pathways, without disease-group contrasts. Other named genera and unclassified sources retain all remaining mass.',
    'The complete three-panel figure is a compact 210 x 310 mm reading master with no panel titles. Panels b/c share a single pathway-label column and aligned rows. Visible major-class headings and fine subclass separators preserve functional organisation; all detailed class names remain in the selection/order table. B/C-only masters and adult-81 reading pages are superseded. Original complete primary-discovery and interaction supplemental maps remain unchanged.',
    'The complete primary-discovery and generation-interaction heatmaps contain every discovery, with pagination only for legibility. Full adjusted ratios, confidence intervals, original raw P, BH q and model diagnostics are in the complete searchable tables.',
    '', '## Methods and interpretation',
    'The original RDS/HDF5 are read-only. All 620 supplied unstratified community rows form a fixed sample denominator before the >=10% nonzero prevalence and nonzero-variance filter. No selected-pathway re-closure is used. The original export does not establish the exact HUMAnN version/raw units and has no UNMAPPED/UNINTEGRATED rows.',
    'Fresh official MaAsLin2 v1.24.1 mixed models use NONE normalization after the fixed-denominator conversion, official LOG zero replacement (half the smallest positive value within a feature/model), an unstandardized disease indicator and a primary-family random intercept. The ratio 2^beta refers to zero-replaced relative abundance, not absolute microbial load or measured metabolic flux.',
    'Model 3 adjusts for generation, generation-specific natural age splines (3 df), sex, preparation batch, sequencing pool, cyclic sampling season and antibiotic use in the preceding three months. Model 4 additionally includes smoking, alcohol, diet score, physical activity, log total energy and continuous BMI. Adult Model 3/4 use identical complete cases. Exactly aliased/constant nuisance terms are removed by QR, never by association P values.',
    'All eligible pathways are fitted; complete BH families are retained for each overall effect, the pooled three-generation effects and the joint two-df Satterthwaite generation test. Invalid fits remain flagged and count conservatively as P=1 internally; no unreliable estimate is declared significant. Singular random-intercept fits are flagged separately from failed convergence.',
    'Bray-Curtis functional dbRDA conditions on family and the corresponding nuisance design, with a Lingoes correction. Restricted reduced-model permutations co-permute the conditional response and nuisance terms with re-residualization of exposure. Every accelerated run is checked against ten compact and one original vegan v2.7-6 calculations, tolerance 1e-7.',
    'The historical lineage file was publicly shared by Eric Franzosa in 2021. It is not the current MetaCyc release. Its URL, date, hash, multiple lineages and unmatched pathways are documented. Classification counts overlap and are not independent mechanisms or enrichment results.',
    'HUMAnN community pathway quantities need not equal the sum of species-stratified quantities. The community-versus-strata table audits this empirically; no source percentage uses community totals as its denominator. Source taxon names are retained as HUMAnN annotations, without forced SGB or strain assignment. Source contrasts have no inferential P values.',
    'Pathway names are supplied database annotations. For example, PWY-7383 is labelled "anaerobic energy metabolism (invertebrates, cytosol)" and its observed stratified quantity is almost entirely unclassified. This label alone does not establish an invertebrate organism, a host tissue of origin or measured metabolic activity in these samples.',
    sprintf('An exact binomial-only crosswalk links %d source taxon labels uniquely to existing named SGBs; %d distinct linked SGBs are existing primary steatosis discoveries. Multiple SGBs with the same binomial, unresolved names and absent exact matches remain explicit, without nearest-name or genus-to-species imputation.',
      r$source_taxonomy_context$unique_matches,r$source_taxonomy_context$matched_primary_SGBs),
    '', '## Literature and adaptations',
    'The structure follows functional discovery -> biological organisation -> microbial sources in Mei et al. (Nature Medicine 2024, Fig. 3; doi:10.1038/s41591-024-03067-7) and Kim et al. (Nature Metabolism 2025, Fig. 3; doi:10.1038/s42255-025-01318-6). The source papers use different populations, filters, thresholds and/or additional data layers; this is an explicitly adapted design, not an exact replication.',
    'MaAsLin2: Mallick et al. PLoS Computational Biology 2021;17:e1009442. doi:10.1371/journal.pcbi.1009442.',
    'HUMAnN framework: Beghini et al. eLife 2021;10:e65088. doi:10.7554/eLife.65088. Pathway stratification/additivity: https://github.com/biobakery/humann#3-pathway-abundance-file',
    'Functional integration context: Lloyd-Price et al. Nature 2019;569:655-662. doi:10.1038/s41586-019-1237-9; Vatanen et al. Nature 2018;562:589-594. doi:10.1038/s41586-018-0620-2; Asnicar et al. Nature Medicine 2021;27:321-332. doi:10.1038/s41591-020-01183-8.',
    'This is cross-sectional functional potential, not transcription, metabolite production, strain sharing or causal mediation.')
  writeLines(lines,file.path(out,'reports/Methods_results_and_figure_legend_EN.md'),useBytes=TRUE)
  invisible(NULL)
}

pfs_validate <- function(r,cfg,out) {
  checks<-list();add<-function(name,ok,detail=''){
    checks[[length(checks)+1L]]<<-data.frame(check=name,pass=isTRUE(ok),detail=as.character(detail))
  }
  p<-pfs_primary(r);add('Primary sample preserved',unique(p$N)==3193&&unique(p$cases)==640)
  add('Adult paired complete samples identical',identical(r$fits$Adult_B1_same_sample$sample_ids,r$fits$B2$sample_ids))
  for(nm in names(r$fits))for(part in c('effects','heterogeneity')){
    z<-r$fits[[nm]][[part]];if(!nrow(z))next
    good<-tf_bool(z$analysis_valid)&is.finite(z$P_raw)
    pv<-rep(1,nrow(z));pv[good]<-z$P_raw[good];expected<-p.adjust(pv,'BH',n=nrow(z))
    err<-if(any(good))max(abs(z$P_BH[good]-expected[good]))else 0
    add(paste('Complete BH recomputation',nm,part),err<1e-12,err)
    add(paste('Invalid estimates not discoveries',nm,part),!any(tf_bool(z$FDR_lt_005)&!good))
    if(part=='effects'){
      valid<-tf_bool(z$analysis_valid)
      add(paste('Ratios and CI ordering',nm),all(abs(log2(z$abundance_ratio[valid])-z$beta_log2[valid])<1e-10)&
        all(z$ratio_CI_low[valid]<z$abundance_ratio[valid]&z$abundance_ratio[valid]<z$ratio_CI_high[valid]))
      add(paste('Official package match',nm),all(tf_bool(z$official_result_match[valid])))
    }
  }
  pred<-p[tf_bool(p$FDR_lt_005),,drop=FALSE]
  a<-r$classification$annotation;pred$parent<-a$display_parent[match(pred$feature_id,a$feature_id)]
  pred<-pred[order(pred$P_BH,pred$P_raw,pred$feature_id),,drop=FALSE]
  expected<-head(pred$feature_id[!duplicated(pred$parent)],10)
  add('Legacy ten-representative selection retained for provenance',identical(r$representatives$feature_id,expected))
  adult<-r$effects[r$effects$model=='B2'&r$effects$generation=='Overall'&
    tf_bool(r$effects$analysis_valid)&tf_bool(r$effects$FDR_lt_005),,drop=FALSE]
  display<-r$adult81_display;selection<-display$selection
  add('Main B/C contain all 81 original Model-4 discoveries',nrow(selection)==81L&&
    !anyDuplicated(selection$feature_id)&&setequal(selection$feature_id,adult$feature_id))
  add('Main B/C functional rows share a monotone order',all(diff(selection$y)<0)&&
    all(diff(selection$group_index)>=0)&&!anyDuplicated(selection$y))
  add('All unmatched pathways retained last',sum(!selection$mapped)==5L&&
    all(selection$mapped[seq_len(sum(selection$mapped))]))
  canonical<-r$classification$lineages
  add('Displayed categories use recorded true lineages only',all(vapply(which(selection$mapped),function(i)
    any(canonical$feature_id==selection$feature_id[i]&canonical$lineage==selection$display_lineage[i],na.rm=TRUE),logical(1))))
  key<-paste(r$effects$model,r$effects$generation,r$effects$feature_id)
  e<-display$effects;at<-match(paste(e$model,e$generation,as.character(e$feature_id)),key)
  available<-!is.na(at)&tf_bool(e$analysis_valid)
  add('Displayed effects and original P/q unchanged',all(e$beta_log2[available]==r$effects$beta_log2[at[available]])&&
    all(e$P_raw[available]==r$effects$P_raw[at[available]])&&all(e$P_BH[available]==r$effects$P_BH[at[available]]))
  add('B/C y coordinates identical for every repeated row',
    all(e$y==selection$y[match(as.character(e$feature_id),selection$feature_id)])&&
    all(display$source$y==selection$y[match(display$source$feature_id,selection$feature_id)]))
  mass<-stats::aggregate(source_share_percent~feature_id+group,display$source,sum)
  add('Adult-81 source display retains all source mass',all(abs(mass$source_share_percent-100)<1e-8))
  add('Source descriptive sample unchanged',display$source_sample_N==3193L)
  lookup<-pfs_adult81_panels(r,cfg)$source
  add('Finite source percentages remain finite in colour lookup',
    all(is.finite(lookup$display_share_percent[is.finite(lookup$source_share_percent)])))
  add('Colour lookup changes machine roundoff only',
    max(abs(lookup$display_share_percent-lookup$source_share_percent),na.rm=TRUE)<1e-10)
  add('All supplied annotations retained',nrow(a)==620&&!anyDuplicated(a$feature_id))
  for(v in c('vegan_kernel_error','full_original_vegan_error'))add(paste('Permutation equivalence',v),
    all(r$composition$diagnostics[[v]]<1e-7),max(r$composition$diagnostics[[v]]))
  add('Permutation six scopes complete',nrow(r$composition$effects)==6L&&all(r$composition$effects$permutations==9999L))
  add('All actual source rows retained',r$sources$source_rows==33723L)
  add('Source group size preserved',sum(r$sources$groups$N[r$sources$groups$group!='All'])==3193L)
  sums<-stats::aggregate(source_share_percent~feature_id+group,r$sources$genus,FUN=sum)
  add('All-genus source fractions sum to 100',all(abs(sums$source_share_percent-100)<1e-8))
  sums<-stats::aggregate(source_share_percent~feature_id+group,r$sources$display,FUN=sum)
  add('Display retained Other/unclassified source mass',all(abs(sums$source_share_percent-100)<1e-8))
  add('Source chunk and genus aggregation equivalent',r$sources$max_aggregation_error<1e-9,r$sources$max_aggregation_error)
  cross<-r$source_taxonomy_context$crosswalk
  add('All source taxon labels retained in SGB crosswalk',nrow(cross)==r$sources$taxon_strata&&!anyDuplicated(cross$taxon_label))
  add('No forced ambiguous taxon match',all(is.na(cross$feature_id[cross$mapping_status!='Unique exact binomial match'])))
  # The previous module is an independent run. Comparison is a numerical audit,
  # never reuse of its fitted output or an enforced discovery-count target.
  oldmod<-file.path(cfg$project_root,'outputs/42_Functional_pathways/current')
  oldpaths<-list.files(oldmod,pattern='analysis_result.rds',recursive=TRUE,full.names=TRUE)
  if(length(oldpaths)){
    old<-readRDS(oldpaths[which.max(file.info(oldpaths)$mtime)])
    k<-match(paste(r$effects$model,r$effects$generation,r$effects$feature_id),
      paste(old$effects$model,old$effects$generation,old$effects$feature_id))
    valid<-tf_bool(r$effects$analysis_valid)&tf_bool(old$effects$analysis_valid[k])&!is.na(k)
    delta<-max(abs(r$effects$beta_log2[valid]-old$effects$beta_log2[k[valid]]))
    ade_table(data.frame(comparison='Fresh module-44 vs independently completed module-42 identical model',
      matched_valid_effects=sum(valid),maximum_beta_difference=delta), 'Independent_fresh_fit_comparison_EN',out)
    add('Fresh repeated inference numerically reproducible',delta<1e-6,delta)
  }
  export<-read.csv(file.path(out,'tables/Figure_export_manifest_EN.csv'),check.names=FALSE)
  for(stem in export$figure){
    fs<-file.path(out,'figures',paste0(stem,c('.pdf','.png','.tiff','.svg')))
    add(paste('Four formats',stem),all(file.exists(fs))&&all(file.info(fs)$size>100))
    font<-jsonlite::fromJSON(file.path(out,'audit',paste0(stem,'.effective-text-audit.json')))
    collision<-jsonlite::fromJSON(file.path(out,'audit',paste0(stem,'.collision-audit.json')))
    alignment<-jsonlite::fromJSON(file.path(out,'audit',paste0(stem,'.alignment.json')))
    add(paste('Rendered font gate',stem),identical(font$verdict,'PASS')&&font$minimum_pt>=5)
    add(paste('Rendered collision gate',stem),identical(collision$verdict,'PASS')&&
      collision$summary$fail==0L&&collision$summary$warn==0L)
    add(paste('Rendered alignment gate',stem),alignment$summary$fail==0L&&alignment$summary$warn==0L)
  }
  for(nm in c('primary_pathway_discoveries','generation_interactions')){
    manifest<-tf_read(file.path(out,'tables',paste0('FigureS44_All_',nm,'_complete_page_manifest_EN.csv')))
    expected<-if(nm=='primary_pathway_discoveries')p$feature_id[tf_bool(p$FDR_lt_005)]else
      r$heterogeneity$feature_id[r$heterogeneity$model=='B1'&tf_bool(r$heterogeneity$FDR_lt_005)]
    add(paste('All discoveries plotted exactly once',nm),!anyDuplicated(manifest$feature_id)&&setequal(manifest$feature_id,expected))
  }
  display_rows<-tf_read(file.path(out,'tables/Figure44_Adult_Model4_all81_selection_order_EN.csv'))
  add('All 81 pathways appear exactly once in combined heatmaps',nrow(display_rows)==81L&&
    identical(as.character(display_rows$feature_id),as.character(selection$feature_id)))
  add('Superseded B/C-only files absent from current display',
    !any(grepl('^Figure(44bc_All81_|S44_Adult81_)',export$figure)))
  status<-system2(cfg$python,c(shQuote(file.path(cfg$project_root,'tools/verify_pathway_export_metadata.py')),shQuote(out)))
  add('Independent raster DPI/dimensions and vector page headers',identical(status,0L))
  ade_protect(r$protected);add('Original source fingerprints unchanged',TRUE)
  t<-ov_bind(checks);fl_write_csv(t,file.path(out,'audit/independent_validation.csv'))
  if(!all(t$pass))stop('Independent validation failed: ',paste(t$check[!t$pass],collapse='; '))
  ade_log('Independent validation PASS: ',nrow(t),' checks')
  invisible(t)
}
