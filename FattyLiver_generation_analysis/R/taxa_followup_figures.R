tf_export <- function(spec,stem,out,cfg) {
  spec$plot<-fl_model_plot_labels(spec$plot)
  audit<-file.path(out,'audit');width<-spec$width;height<-spec$height
  probe<-file.path(audit,paste0(stem,'.layout-probe.pdf'))
  grDevices::cairo_pdf(probe,width=width/25.4,height=height/25.4,family='Arial',bg='white')
  opened<-TRUE;on.exit(if(opened)grDevices::dev.off(),add=TRUE)
  grob<-patchwork::patchworkGrob(spec$plot)
  grob<-fl_model_label_grob(grob)
  if(isTRUE(spec$uniform_y_axes)){
    # Free y scales can produce different tick-label widths in each column.
    # Equalize only these fixed axis cells; keep the data, scales, font sizes
    # and flexible panel widths unchanged. Audit the final rendered geometry.
    axes<-grob$layout[grepl('^axis-l-[0-9]+$',grob$layout$name),,drop=FALSE]
    stopifnot(nrow(axes)>0L,all(axes$l==axes$r))
    columns<-sort(unique(axes$l))
    axis_width<-max(vapply(columns,function(j)grid::convertWidth(grob$widths[j],'bigpts',valueOnly=TRUE),numeric(1)))
    grob$widths[columns]<-grid::unit(rep(axis_width,length(columns)),'bigpts')
  }
  grid::grid.newpage();grid::grid.draw(grob);grid::grid.force()
  rows<-grob$layout[grepl('^panel-[0-9]+$',grob$layout$name),,drop=FALSE]
  # patchwork nests a single coord_fixed/equal panel in a combined panel cell.
  # Opt-in for the axis-free reference-tree canvas; other modules are unchanged.
  fixed_single<-isTRUE(spec$fixed_aspect_single_panel)&&spec$panels==1L&&nrow(rows)==0L
  if(fixed_single)rows<-grob$layout[grepl('^panel;',grob$layout$name),,drop=FALSE]
  rows<-rows[order(rows$t,rows$l),,drop=FALSE]
  stopifnot(nrow(rows)==spec$panels)
  vps<-grid::grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)$name
  measure<-function(name,anchor=FALSE){hit<-vps[startsWith(vps,paste0(name,'.'))]
    if(length(hit)!=1L)stop('Cannot locate rendered viewport: ',name)
    grid::seekViewport(hit);on.exit(grid::upViewport(0))
    lo<-grid::deviceLoc(grid::unit(0,'npc'),grid::unit(0,'npc'),valueOnly=TRUE)
    hi<-grid::deviceLoc(grid::unit(1,'npc'),grid::unit(1,'npc'),valueOnly=TRUE)
    if(anchor)return(72*c(lo$x,hi$y));72*c(lo$x,lo$y,hi$x,hi$y)
  }
  panels<-lapply(seq_len(nrow(rows)),function(i){r<-rows[i,];k<-if(fixed_single)1L else as.integer(sub('panel-','',r$name))
    list(id=letters[k],bbox_pt=measure(r$name),grid_id='followup-patchwork',
      row_start=r$t-1L,row_stop=r$b,col_start=r$l-1L,col_stop=r$r,
      panel_label=if(isTRUE(spec$no_panel_tags))'' else letters[k],
      panel_label_anchor_pt=measure(if(isTRUE(spec$no_panel_tags))r$name else paste0('tag-',k),TRUE))})
  jsonlite::write_json(list(schema_version=1,backend='R-patchwork-Cairo',measurement='Final device rendered panel viewports',
    figure=list(width_pt=width/25.4*72,height_pt=height/25.4*72),panels=panels,exemptions=spec$exemptions),
    file.path(audit,paste0(stem,'.alignment-layout.json')),auto_unbox=TRUE,pretty=TRUE,digits=NA)
  grDevices::dev.off();opened<-FALSE
  status<-system2(cfg$python,c(shQuote(file.path(nf_skill_root(),'scripts/audit_panel_alignment.py')),
    shQuote(file.path(audit,paste0(stem,'.alignment-layout.json'))),'--json-out',shQuote(file.path(audit,paste0(stem,'.alignment.json'))),'--strict'))
  if(status!=0L)stop('Rendered panel alignment blocked export: ',stem)
  draw<-function(){grid::grid.newpage();grid::grid.draw(grob)}
  grDevices::png(file.path(audit,paste0(stem,'.preview.png')),width=width,height=height,units='mm',res=300,type='cairo',bg='white')
  tryCatch(draw(),finally=grDevices::dev.off())
  # Export editable vector PDF first, then require physical text/collision PASS.
  pdf<-file.path(out,'figures',paste0(stem,'.pdf'))
  grDevices::cairo_pdf(pdf,width=width/25.4,height=height/25.4,family='Arial',bg='white')
  tryCatch(draw(),finally=grDevices::dev.off())
  text_contract_args <- character()
  if(!is.null(spec$max_body_pt)) {
    if(width<=183||is.null(spec$large_format_font_reason)||!nzchar(spec$large_format_font_reason))
      stop('A larger physical text ceiling requires an explicit large-format artifact contract')
    stopifnot(is.finite(spec$max_body_pt),spec$max_body_pt>=7.01)
    text_contract_args <- c('--max-body-pt',as.character(spec$max_body_pt))
    jsonlite::write_json(list(reason=spec$large_format_font_reason,width_mm=width,height_mm=height,
      glyph_floor_pt=5,maximum_text_pt=spec$max_body_pt,other_figures_default_maximum_pt=7.01),
      file.path(audit,paste0(stem,'.artifact-font-contract.json')),auto_unbox=TRUE,pretty=TRUE)
  }
  status<-system2(cfg$python,c(shQuote(file.path(cfg$project_root,'scripts/audit_nature_cairo_pdf.py')),shQuote(pdf),
    '--skill',shQuote(nf_skill_root()),'--prefix',shQuote(file.path(audit,stem)),
    if(isTRUE(spec$oriented_text))'--oriented-text'else character(),text_contract_args))
  if(status!=0L)stop('Rendered text/collision blocked export: ',stem)
  svglite::svglite(file.path(out,'figures',paste0(stem,'.svg')),width=width/25.4,height=height/25.4,bg='white')
  tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::png(file.path(out,'figures',paste0(stem,'.png')),width=width,height=height,units='mm',res=600,type='cairo',bg='white')
  tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::tiff(file.path(out,'figures',paste0(stem,'.tiff')),width=width,height=height,units='mm',res=600,type='cairo',compression='lzw',bg='white')
  tryCatch(draw(),finally=grDevices::dev.off())
  data.frame(figure=stem,width_mm=width,height_mm=height,panels=spec$panels,formats='PDF;PNG;TIFF;SVG',raster_dpi=600)
}

tf_plot_theme <- function() nf_theme()+ggplot2::theme(
  plot.title=ggplot2::element_text(size=6.8,margin=ggplot2::margin(l=12,b=6,unit='pt')),
  legend.background=ggplot2::element_blank(),legend.key=ggplot2::element_blank(),
  legend.box.background=ggplot2::element_blank(),plot.margin=ggplot2::margin(6,7,5,5,'pt'))

tf_scatter <- function(z,title,xlabel,ylabel,tag,cfg) {
  z$direction<-factor(ifelse(z$beta_log2_x<0,'Lower with steatosis','Higher with steatosis'),levels=c('Lower with steatosis','Higher with steatosis'))
  rho<-cor(z$beta_log2_x,z$beta_log2_y,method='spearman');lims<-range(c(z$beta_log2_x,z$beta_log2_y),finite=TRUE);pad<-diff(lims)*.08;lims<-lims+c(-pad,pad)
  ggplot2::ggplot(z,ggplot2::aes(beta_log2_x,beta_log2_y))+
    ggplot2::geom_hline(yintercept=0,colour='#BBBBBB',linewidth=.25)+ggplot2::geom_vline(xintercept=0,colour='#BBBBBB',linewidth=.25)+
    ggplot2::annotate('segment',x=lims[1],y=lims[1],xend=lims[2],yend=lims[2],colour='#777777',linetype=2,linewidth=.25)+
    ggplot2::geom_point(ggplot2::aes(fill=direction),shape=21,colour='#444444',stroke=.15,size=1.3,alpha=.8)+
    ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),levels(z$direction)))+
    ggplot2::scale_x_continuous(limits=lims,breaks=scales::breaks_pretty(4))+
    ggplot2::scale_y_continuous(limits=lims,breaks=scales::breaks_pretty(4))+
    ggplot2::annotate('text',x=lims[1]+diff(lims)*.04,y=lims[2]-diff(lims)*.04,hjust=0,vjust=1,
      label=sprintf('n = %d SGBs\nSpearman rho = %.2f',nrow(z),rho),size=5.7/(72.27/25.4))+
    ggplot2::labs(title=title,x=xlabel,y=ylabel,tag=tag)+tf_plot_theme()
}

tf_figures25 <- function(result,cfg,out) {
  a<-result$comparisons[result$comparisons$comparison=='Adult B1 versus adult B2: same participants',]
  b<-result$comparisons[result$comparisons$comparison=='Overall B1 versus antibiotic exclusion: different samples',]
  p1<-tf_scatter(a,'Adult lifestyle and BMI adjustment','Adult B1 coefficient (log2)','Adult B2 coefficient (log2)','a',cfg)
  p2<-tf_scatter(b,'Exclusion of recent antibiotic use','Overall B1 coefficient (log2)','Antibiotic-excluded coefficient (log2)','b',cfg)
  fig<-patchwork::wrap_plots(p1,p2,nrow=1,guides='keep')
  one<-tf_export(list(plot=fig,width=183,height=100,panels=2,exemptions=list()),'Figure25_Sensitivity_effect_comparisons',out,cfg)
  t<-result$tables$B1;t<-t[t$FDR_lt_005&tf_bool(t$analysis_valid),];t<-head(t[order(t$P_BH,t$P_raw,t$feature_id),],16L)
  t$display_order<-seq_len(nrow(t));t$SGB_id<-factor(t$SGB_id,levels=rev(t$SGB_id))
  t$direction<-factor(ifelse(t$beta_log2<0,'Lower with steatosis','Higher with steatosis'),levels=c('Lower with steatosis','Higher with steatosis'))
  p<-ggplot2::ggplot(t,ggplot2::aes(abundance_ratio,SGB_id))+
    ggplot2::geom_vline(xintercept=1,linetype=2,colour='#888888',linewidth=.25)+
    ggplot2::geom_segment(ggplot2::aes(x=ratio_CI_low,xend=ratio_CI_high,yend=SGB_id),linewidth=.4,colour='black')+
    ggplot2::geom_point(ggplot2::aes(fill=direction),shape=21,size=2,colour='black',stroke=.25)+
    ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),levels(t$direction)))+
    ggplot2::scale_x_log10(breaks=c(.125,.25,.5,1,2,4,8,16))+
    ggplot2::labs(title=NULL,x='Adjusted abundance ratio (95% CI)',y=NULL,tag='a')+tf_plot_theme()
  two<-tf_export(list(plot=patchwork::wrap_plots(p),width=89,height=115,panels=1,exemptions=list()),'Figure25_Overall_abundance_ratios',out,cfg)
  fl_write_csv(t,file.path(out,'tables/Figure25_displayed_overall_SGBs_EN.csv'))
  writeLines(c('# Figure legends','',
    '**Figure 25. Sensitivity comparisons of steatosis-associated microbial abundance.** Each point is an SGB with a valid estimate in both compared models. (a) Adult B1 and lifestyle/BMI-adjusted B2 use the same participants. (b) Overall B1 and the antibiotic-excluded analysis use different samples. Dashed lines indicate equal coefficients; colours indicate the overall/reference-model direction. Spearman correlations are descriptive, not independent validation.',
    '**Companion figure. Overall steatosis-associated abundance ratios.** Up to 16 overall B1 BH-supported SGBs ranked by q, raw P and ID are displayed. Points and bars show 2^beta and 95% confidence intervals from official MaAsLin2 log2 relative-abundance mixed models. Ratios refer to zero-replaced relative abundance. B1 adjusts for generation, generation-specific age splines, sex, library batch, sequencing pool, sampling season and recent antibiotics, with a primary-family random intercept. Full results, including nonsignificant and filtered SGBs, are retained in the tables.'),file.path(out,'reports/Figure_legends_EN.md'))
  ov_bind(list(one,two))
}

tf_figures26 <- function(result,cfg,out) {
  t<-result$results[result$results$model=='P2'&result$results$analysis_valid,]
  ranked<-aggregate(P_BH~feature_id,t,min);ranked<-ranked[order(ranked$P_BH,ranked$feature_id),];ids<-head(ranked$feature_id,16L)
  d<-t[t$feature_id%in%ids,];d$position<-match(d$feature_id,rev(ids))
  levels_pair<-c('G0-G1','G1-G2','G0-G2');d$comparison<-factor(d$comparison,levels=levels_pair)
  d$offset<-c(-.22,0,.22)[as.integer(d$comparison)];d$y<-d$position+d$offset
  colors<-setNames(unname(cfg$generation_colors),levels_pair)
  p1<-ggplot2::ggplot(d,ggplot2::aes(younger_ratio_per_older_doubling,y))+
    ggplot2::geom_vline(xintercept=1,colour='#888888',linetype=2,linewidth=.25)+
    ggplot2::geom_segment(ggplot2::aes(x=ratio_CI_low,xend=ratio_CI_high,yend=y,colour=comparison),linewidth=.4)+
    ggplot2::geom_point(ggplot2::aes(fill=comparison),shape=21,colour='black',stroke=.2,size=1.5)+
    ggplot2::scale_fill_manual(values=colors)+ggplot2::scale_colour_manual(values=colors,guide='none')+
    ggplot2::scale_x_log10(breaks=c(.8,1,1.25,1.5,2,3,4))+
    ggplot2::scale_y_continuous(breaks=seq_along(ids),labels=sub('^t__','',rev(ids)),expand=ggplot2::expansion(add=.65))+
    ggplot2::labs(title='Cross-generation abundance concordance',x='Younger abundance ratio per older doubling',y=NULL,tag='a')+tf_plot_theme()
  s<-result$summary[result$summary$model=='P2',]
  counts<-rbind(data.frame(comparison=s$comparison,direction='Positive concordance',N=s$BH_positive),data.frame(comparison=s$comparison,direction='Negative association',N=s$BH_negative))
  counts$comparison<-factor(counts$comparison,levels=levels_pair);counts$direction<-factor(counts$direction,levels=c('Negative association','Positive concordance'))
  counts$N_offset<-ifelse(counts$direction=='Negative association',-counts$N,counts$N)
  p2<-ggplot2::ggplot(counts,ggplot2::aes(N_offset,comparison,fill=direction))+
    ggplot2::geom_vline(xintercept=0,colour='#888888',linewidth=.25)+
    ggplot2::geom_col(width=.56,colour='black',linewidth=.2)+
    ggplot2::geom_text(ggplot2::aes(label=ifelse(N>0,N,''),hjust=ifelse(N_offset<0,1.15,-.15)),size=5.7/(72.27/25.4))+
    ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),levels(counts$direction)))+
    ggplot2::scale_x_continuous(labels=abs,expand=ggplot2::expansion(mult=c(.2,.25)))+
    ggplot2::scale_y_discrete(limits=rev(levels_pair))+
    ggplot2::labs(title='All tested SGBs: BH q < 0.05',x='Number of SGBs',y=NULL,tag='b')+tf_plot_theme()
  spec<-list(plot=patchwork::wrap_plots(p1,p2,nrow=1),width=183,height=125,panels=2,exemptions=list())
  manifest<-tf_export(spec,'Figure26_Taxon_cross_generation_concordance',out,cfg)
  fl_write_csv(d,file.path(out,'tables/Figure26_displayed_SGBs_EN.csv'))
  writeLines(c('# Figure legend','',
    '**Figure 26. Single-taxon abundance concordance across observed generations.** (a) Family-cluster GEE estimates for 16 SGBs selected only for display by minimum P2 BH q across the three comparisons. Points and bars show the younger abundance ratio per doubling of older abundance and 95% confidence intervals. (b) Counts of positive and negative BH-supported associations among all eligible tested SGBs. BH correction pools all eligible SGB-by-comparison tests within P2. Both members age splines, sex, library batch, sequencing pool, season and recent antibiotics are adjusted; P1 results appear in the tables. Dyads are clustered by primary-family label. Positive effects indicate abundance concordance on a zero-replaced relative-abundance scale, not strain transmission or verified parent-child inheritance.'),file.path(out,'reports/Figure_legends_EN.md'))
  manifest
}

tf_figures <- function(result,analysis,cfg,out) {
  manifest<-if(analysis=='taxa_robustness')tf_figures25(result,cfg,out)else tf_figures26(result,cfg,out)
  fl_write_csv(manifest,file.path(out,'tables/Figure_export_manifest_EN.csv'));invisible(manifest)
}

