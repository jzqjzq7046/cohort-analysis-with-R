# R-only structural adaptation of the adopted module-25/26 plotting helpers.
# Outcome evidence: observed older steatosis background and unaffected younger
# SGB abundances. Binary colours are fixed; neutral points are nonsignificant.
fp_figures <- function(r,cfg,out) {
  t<-r$results[r$results$model=='B1'&r$results$analysis_valid,,drop=FALSE]
  keys<-c('G0-G1','G1-G2','G0-G2');stopifnot(all(keys%in%t$comparison))
  t$plot_group<-factor(ifelse(!t$FDR_lt_005,'BH >= 0.05',ifelse(t$beta_log2<0,'Lower with older case','Higher with older case')),
    levels=c('BH >= 0.05','Lower with older case','Higher with older case'))
  pal<-c('BH >= 0.05'='#C6C6C6','Lower with older case'=unname(cfg$binary_colors['Reference']),
    'Higher with older case'=unname(cfg$binary_colors['Comparison']))
  maxx<-max(abs(t$beta_log2))*1.08;maxy<-max(-log10(pmax(t$P_raw,1e-300)))*1.15
  plots<-lapply(seq_along(keys),function(i){z<-t[t$comparison==keys[i],,drop=FALSE]
    ggplot2::ggplot(z,ggplot2::aes(beta_log2,-log10(pmax(P_raw,1e-300))))+
      ggplot2::geom_vline(xintercept=0,colour='#999999',linetype=2,linewidth=.25)+
      ggplot2::geom_point(ggplot2::aes(fill=plot_group),shape=21,colour='#555555',stroke=.1,size=1.1,alpha=.8)+
      ggplot2::scale_fill_manual(values=pal,drop=TRUE)+
      ggplot2::scale_x_continuous(limits=c(-maxx,maxx),breaks=scales::breaks_pretty(4))+
      ggplot2::scale_y_continuous(limits=c(0,maxy),breaks=scales::breaks_pretty(4),expand=ggplot2::expansion(mult=c(0,.02)))+
      ggplot2::labs(title=keys[i],x='Adjusted coefficient (log2)',y='-log10(raw P)',tag=letters[i])+tf_plot_theme()+
      ggplot2::theme(legend.position='none',axis.title.x=ggplot2::element_text(margin=ggplot2::margin(t=8,unit='pt')),
        axis.text.x=ggplot2::element_text(margin=ggplot2::margin(t=3,unit='pt')))
  })
  fig<-patchwork::wrap_plots(plots,nrow=1,guides='keep')
  manifest<-list(tf_export(list(plot=fig,width=183,height=90,panels=3,exemptions=list()),
    'Figure27_Family_background_all_SGBs',out,cfg))
  displayed<-list();forests<-lapply(seq_along(keys),function(i){
    z<-t[t$comparison==keys[i],,drop=FALSE];z<-head(z[order(z$P_BH,z$P_raw,z$feature_id),,drop=FALSE],8L)
    z$SGB_id<-factor(z$SGB_id,levels=rev(z$SGB_id));displayed[[i]]<<-z
    ggplot2::ggplot(z,ggplot2::aes(abundance_ratio,SGB_id))+
      ggplot2::geom_vline(xintercept=1,colour='#999999',linetype=2,linewidth=.25)+
      ggplot2::geom_segment(ggplot2::aes(x=ratio_CI_low,xend=ratio_CI_high,yend=SGB_id),colour='black',linewidth=.35)+
      ggplot2::geom_point(ggplot2::aes(fill=plot_group),shape=21,colour='black',size=1.8,stroke=.2)+
      ggplot2::scale_fill_manual(values=pal,drop=TRUE)+
      ggplot2::scale_x_log10(breaks=c(.0625,.125,.25,.5,1,2,4,8,16,32))+
      ggplot2::labs(title=keys[i],x='Abundance ratio (95% CI)',y=NULL,tag=letters[i])+tf_plot_theme()+
      ggplot2::theme(legend.position='none',axis.title.x=ggplot2::element_text(margin=ggplot2::margin(t=8,unit='pt')),
        axis.text.x=ggplot2::element_text(margin=ggplot2::margin(t=3,unit='pt')))
  })
  forest_figure<-patchwork::wrap_plots(forests,nrow=1,guides='keep')
  manifest[[2]]<-tf_export(list(plot=forest_figure,width=183,height=105,panels=3,exemptions=list()),
    'Figure27_Representative_abundance_ratios',out,cfg)
  v3_table(ov_bind(displayed),'Figure27_Displayed_SGBs_EN',out)
  # Two representative raw-distribution panels. Selection is for display only;
  # all raw participant-level plotting data remain private.
  z<-head(t[t$comparison=='G1-G2',][order(t$P_BH[t$comparison=='G1-G2'],t$P_raw[t$comparison=='G1-G2'],t$feature_id[t$comparison=='G1-G2']),],2L)
  x<-readRDS(file.path(out,'private/analysis_population.rds'));md<-readRDS(file.path(out,'private/G1-G2_B1_design.rds'))$metadata
  raw<-ov_bind(lapply(z$feature_id,function(id)data.frame(feature_id=id,SGB_id=sub('^t__','',id),sample_id=md$sample_id,
    older_case=factor(md$older_any_case,levels=0:1,labels=c('None','At least one')),relative_percent=100*x$relative[id,match(md$sample_id,colnames(x$relative))])))
  saveRDS(raw,file.path(out,'private/representative_distribution_source.rds'))
  pp<-lapply(seq_along(z$feature_id),function(i){d<-raw[raw$feature_id==z$feature_id[i],,drop=FALSE]
    positive<-d$relative_percent[d$relative_percent>0];sigma<-min(positive)/2
    lo<-ceiling(log10(sigma*10));hi<-floor(log10(max(positive)))
    ticks<-c(0,if(lo<=hi)10^seq(lo,hi,by=1)else max(positive))
    ticks<-ticks[ticks<=max(positive)*1.1]
    ggplot2::ggplot(d,ggplot2::aes(older_case,relative_percent,fill=older_case))+
      ggplot2::geom_boxplot(width=.5,outlier.shape=NA,colour='black',linewidth=.3)+
      ggplot2::geom_point(position=ggplot2::position_jitter(width=.09,seed=20261027),colour='#888888',alpha=.35,size=.45)+
      ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),levels(d$older_case)))+
      ggplot2::scale_y_continuous(breaks=ticks,labels=function(v)vapply(v,function(a)if(a==0)'0'else format(a,scientific=FALSE,trim=TRUE),character(1)),
        expand=ggplot2::expansion(mult=c(0,.08)))+
      ggplot2::coord_transform(y=scales::pseudo_log_trans(sigma=sigma,base=10))+
      ggplot2::labs(title=paste('G1-G2:',sub('^t__','',z$feature_id[i])),x='Observed older steatosis cases',y='Relative abundance (%)',tag=letters[i])+tf_plot_theme()+
      ggplot2::theme(legend.position='none',axis.title.x=ggplot2::element_text(margin=ggplot2::margin(t=8,unit='pt')),
        axis.text.x=ggplot2::element_text(margin=ggplot2::margin(t=3,unit='pt')))
  })
  manifest[[3]]<-tf_export(list(plot=patchwork::wrap_plots(pp,nrow=1),width=183,height=95,panels=2,exemptions=list()),
    'Figure27_Representative_raw_distributions',out,cfg)
  v3_table(ov_bind(lapply(split(raw,list(raw$feature_id,raw$older_case),drop=TRUE),function(d){q<-quantile(d$relative_percent,c(.25,.5,.75));
    data.frame(feature_id=d$feature_id[1],exposure=d$older_case[1],N=nrow(d),prevalence=mean(d$relative_percent>0),
      median_percent=unname(q[2]),Q1_percent=unname(q[1]),Q3_percent=unname(q[3]))})), 'Figure27_Raw_distribution_summaries_EN',out)
  fl_write_csv(ov_bind(manifest),file.path(out,'tables/Figure_export_manifest_EN.csv'))
  writeLines(c('# Figure legends','',
    '**Figure 27. Observed older-generation hepatic steatosis and microbial abundances in unaffected younger members.** All younger participants have no ultrasound-defined steatosis. Each point represents an SGB tested in G0-G1 (a), G1-G2 (b) or G0-G2 (c). The horizontal axis is the official MaAsLin2 log2 abundance coefficient for at least one observed older case versus no case among observed older members; the vertical axis is -log10(raw P). Blue/orange denote negative/positive associations with joint BH q<0.05; grey indicates other valid estimates. BH correction pools all eligible SGB-by-comparison tests within all-age B1. Models adjust younger age splines, sex, library batch, pool, sampling season, recent antibiotics and nonconstant observed older-member count. A family random intercept is used when younger members recur within families; otherwise the official ordinary Gaussian model is used. Primary labels do not verify direct parent-child relationships.',
    '', '**Representative abundance ratios.** Eight SGBs per comparison are selected only for display by joint BH q, raw P and ID. Points and bars show 2^beta and t-based 95% confidence intervals; one is the null. Ratios refer to geometric means of zero-replaced relative abundance, not absolute bacterial load or ORs. Colour rules and B1 covariates match Figure 27; nonsignificant estimates are retained when fewer than eight BH-supported taxa exist. All eligible SGBs and adult extensions are in the full tables.',
    '', '**Representative raw relative-abundance distributions in unaffected G2 members.** The two G1-G2 SGBs ranked first by joint BH q, raw P and ID are displayed regardless of significance. Blue indicates no observed older steatosis case; orange indicates at least one. Boxes show median and interquartile range, whiskers extend to 1.5 times the interquartile range, and grey points are participants. A monotone pseudo-log coordinate axis retains zeros and labels actual relative-abundance percentages; its sigma is half the minimum positive percentage for each displayed SGB. Box statistics are computed on untransformed percentages. A zero-height box indicates coincident quartiles at zero. These are unadjusted distributions; inference comes from the adjusted models, not a separate boxplot test. Sample counts and abundance summaries are in Figure27_Raw_distribution_summaries_EN.csv.'),file.path(out,'reports/Figure_legends_EN.md'),useBytes=TRUE)
  invisible(manifest)
}
