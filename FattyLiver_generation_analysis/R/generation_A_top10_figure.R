# A saved-result display revision, user-selected on 2026-11-01.
# No fitting, multiplicity correction or data-dependent axis selection.
ade_A_fixed_axes <- function()list(
  overview_x_limits=c(1,32),overview_x_breaks=c(1,2,4,8,16,32),
  overview_y_limits=c(0,3),overview_y_breaks=0:3,
  forest_x_limits=c(.01,10),forest_x_breaks=c(.01,.1,1,10))

ade_A_top10_combined <- function(r,cfg,out) {
  stopifnot(identical(r$kind,'A'))
  axes<-ade_A_fixed_axes()
  h<-r$heterogeneity[r$heterogeneity$model=='B1'&tf_bool(r$heterogeneity$analysis_valid),,drop=FALSE]
  e<-r$effects[r$effects$model=='B1'&tf_bool(r$effects$analysis_valid),,drop=FALSE]
  stopifnot(nrow(h)>=10,!anyDuplicated(h$feature_id))
  h<-h[order(h$P_raw,h$feature_id),,drop=FALSE];h$interaction_rank<-seq_len(nrow(h))
  h$effect_spread_ratio<-vapply(h$feature_id,function(id){
    b<-e$beta_log2[e$feature_id==id];stopifnot(length(b)==3L,all(is.finite(b)))
    2^(max(b)-min(b))
  },numeric(1))
  h$minus_log10_BH_q<--log10(pmax(h$P_BH,.Machine$double.xmin))
  h$FDR_display<-factor(ifelse(h$FDR_lt_005,'BH q < 0.05','BH q ≥ 0.05'),levels=c('BH q ≥ 0.05','BH q < 0.05'))
  h$overview_label<-ifelse(h$interaction_rank<=10L,sub('^t__','',h$feature_id),'')
  top<-head(h,10);ids<-top$feature_id
  # The frozen result has seven labels in a dense low-q-height strip. Separate
  # their text into two rows while retaining every original point coordinate.
  # All leader ends sit below the corresponding label's physical text box.
  label_layout<-data.frame(feature_id=top$feature_id,label=top$overview_label,
    interaction_rank=top$interaction_rank,point_x=top$effect_spread_ratio,
    point_y=top$minus_log10_BH_q,joint_interaction_P=top$P_raw,joint_interaction_BH_q=top$P_BH)
  label_layout$label_x<-label_layout$point_x
  label_layout$label_y<-label_layout$point_y+.12
  crowded<-which(label_layout$point_y<1)
  crowded<-crowded[order(label_layout$point_x[crowded],label_layout$feature_id[crowded])]
  stopifnot(length(crowded)==7L)
  label_layout$label_x[crowded]<-2^c(.55,1.25,1.95,2.65,3.35,4.5,4.65)
  label_layout$label_y[crowded]<-c(.85,.98,.85,.98,.85,.98,.85)
  label_layout$leader_y_end<-label_layout$label_y-.06
  forest<-e[e$feature_id%in%ids,,drop=FALSE]
  forest$interaction_rank<-match(forest$feature_id,ids)
  forest$joint_interaction_P<-top$P_raw[forest$interaction_rank]
  forest$joint_interaction_BH_q<-top$P_BH[forest$interaction_rank]
  forest$joint_FDR_lt_005<-top$FDR_lt_005[forest$interaction_rank]
  forest$generation<-factor(forest$generation,levels=c('G0','G1','G2'))
  offset<-c(G0=.21,G1=0,G2=-.21)
  forest$y<-11-forest$interaction_rank+unname(offset[as.character(forest$generation)])
  stopifnot(nrow(forest)==30L,all(table(forest$feature_id)==3L),
    all(h$effect_spread_ratio>=axes$overview_x_limits[1]),all(h$effect_spread_ratio<=axes$overview_x_limits[2]),
    all(h$minus_log10_BH_q>=axes$overview_y_limits[1]),all(h$minus_log10_BH_q<=axes$overview_y_limits[2]),
    all(is.finite(forest$ratio_CI_low)),all(is.finite(forest$ratio_CI_high)),
    all(forest$ratio_CI_low>=axes$forest_x_limits[1]),all(forest$ratio_CI_high<=axes$forest_x_limits[2]))

  theme<-ade_theme()+ggplot2::theme(legend.position='bottom',legend.title=ggplot2::element_blank(),
    legend.text=ggplot2::element_text(size=5.6),legend.key.width=grid::unit(2.4,'mm'),
    legend.spacing.x=grid::unit(.8,'mm'),legend.box.margin=ggplot2::margin(0,0,0,0,'pt'))
  p1<-ggplot2::ggplot(h,ggplot2::aes(effect_spread_ratio,minus_log10_BH_q,colour=FDR_display))+
    ggplot2::geom_hline(yintercept=-log10(.05),linetype=2,colour='#888888',linewidth=.25)+
    ggplot2::geom_point(size=1.25,alpha=.78,shape=16)+
    ggplot2::geom_segment(data=label_layout,ggplot2::aes(x=point_x,y=point_y,
      xend=label_x,yend=leader_y_end),inherit.aes=FALSE,colour='#777777',linewidth=.2)+
    ggplot2::geom_text(data=label_layout,ggplot2::aes(x=label_x,y=label_y,label=label),
      inherit.aes=FALSE,colour='black',family='Arial',size=5.6/(72.27/25.4))+
    ggplot2::scale_colour_manual(values=setNames(unname(cfg$binary_colors),levels(h$FDR_display)),drop=FALSE)+
    ggplot2::scale_x_log10(limits=axes$overview_x_limits,breaks=axes$overview_x_breaks,
      labels=c('1','2','4','8','16','32'),expand=ggplot2::expansion(mult=c(.02,.02)))+
    ggplot2::scale_y_continuous(limits=axes$overview_y_limits,breaks=axes$overview_y_breaks,
      expand=ggplot2::expansion(mult=c(.01,.02)))+
    ggplot2::labs(title=NULL,x='Max/min generation abundance ratio',y='−log10(BH q)')+theme

  labels<-ade_ylabels(r$annotation,ids)
  labels<-as.expression(lapply(seq_along(ids),function(i)
    if(isTRUE(top$FDR_lt_005[i]))bquote(.(labels[[i]])~'*')else labels[[i]]))
  p2<-ggplot2::ggplot(forest,ggplot2::aes(abundance_ratio,y))+
    ggplot2::geom_vline(xintercept=1,linetype=2,colour='#888888',linewidth=.25)+
    ggplot2::geom_segment(ggplot2::aes(x=ratio_CI_low,xend=ratio_CI_high,yend=y),colour='black',linewidth=.28)+
    ggplot2::geom_point(ggplot2::aes(fill=generation),shape=21,size=1.45,colour='black',stroke=.2)+
    ggplot2::scale_fill_manual(values=cfg$generation_colors,drop=FALSE)+
    ggplot2::scale_x_log10(limits=axes$forest_x_limits,breaks=axes$forest_x_breaks,
      labels=c('0.01','0.1','1','10'),expand=ggplot2::expansion(mult=c(.02,.02)))+
    ggplot2::scale_y_continuous(limits=c(.4,10.6),breaks=10:1,labels=labels,expand=ggplot2::expansion(mult=0))+
    ggplot2::labs(title=NULL,x='Adjusted abundance ratio (95% CI)',y=NULL)+theme
  p<-patchwork::wrap_plots(p1,p2,nrow=1,widths=c(1,1))+
    patchwork::plot_annotation(tag_levels='a',caption='* Joint generation-interaction BH q < 0.05',
      theme=ggplot2::theme(plot.tag=ggplot2::element_text(size=8,face='bold',family='Arial'),
        plot.caption=ggplot2::element_text(size=5.6,family='Arial',hjust=1)))
  stem<-'Figure41_Generation_heterogeneity_overview_top10_no_titles'
  top$display_selection<-'Ten smallest primary raw joint-interaction P values, ties by SGB ID; independent of q threshold'
  ade_table(top,'Top10_joint_generation_tests_EN',out)
  ade_table(forest,'Top10_generation_effects_EN',out)
  ade_table(h,'All_SGB_heterogeneity_overview_source_EN',out)
  ade_table(label_layout,'Figure41_Overview_label_manifest_EN',out)
  ax<-data.frame(panel=c('a','b'),scale=c('log10','log10'),
    x_lower=c(axes$overview_x_limits[1],axes$forest_x_limits[1]),
    x_upper=c(axes$overview_x_limits[2],axes$forest_x_limits[2]),
    x_ticks=c(paste(axes$overview_x_breaks,collapse=';'),paste(axes$forest_x_breaks,collapse=';')),
    axis_rule='Fixed literal domain and labelled ticks; no pretty(), data-derived limits, clipping or rescaling')
  ade_table(ax,'Figure41_Fixed_axis_contract_EN',out)
  fl_write_csv(data.frame(figure=stem,display_rank=seq_along(ids),feature_id=ids,
    primary_raw_P=top$P_raw,primary_BH_q=top$P_BH,FDR_lt_005=top$FDR_lt_005,
    selection=top$display_selection),file.path(out,'tables/Figure41_Top10_display_manifest_EN.csv'))
  spec<-list(plot=p,width=183,height=132,panels=2,exemptions=list())
  tf_export(spec,stem,out,cfg)
}
