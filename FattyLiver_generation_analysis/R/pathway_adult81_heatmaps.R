# Inert Figure44 display extension. No model fitting, P values or FDR changes.
pfs_adult81_revision <- 'adult_Model4_all81_compact_combined_20261102_v2'
pfs_class_text <- function(x) {
  labels<-c('Biosynthesis'='Biosynthesis','Degradation'='Degradation',
    'Energy-Metabolism'='Energy metabolism','Metabolic-Clusters'='Metabolic clusters',
    'Super-Pathways'='Superpathways','Unmapped'='Unmapped',
    'Amino-Acid-Biosynthesis'='Amino acid biosynthesis','Nucleotide-Biosynthesis'='Nucleotide biosynthesis',
    'Cofactor-Biosynthesis'='Cofactor biosynthesis','Cell-Structure-Biosynthesis'='Cell structure biosynthesis',
    'Carbohydrates-Biosynthesis'='Carbohydrate biosynthesis','Lipid-Biosynthesis'='Lipid biosynthesis',
    'AROMATIC-COMPOUNDS-BIOSYN'='Aromatic compound biosynthesis',
    'SECONDARY-METABOLITE-BIOSYNTHESIS'='Secondary metabolite biosynthesis',
    'Aminoacyl-tRNAs-Charging'='Aminoacyl-tRNA charging',
    'Carbohydrates-Degradation'='Carbohydrate degradation','Amino-Acid-Degradation'='Amino acid degradation',
    'NUCLEO-DEG'='Nucleotide/nucleoside degradation','CARBOXYLATES-DEG'='Carboxylate degradation',
    'SECONDARY-METABOLITE-DEGRADATION'='Secondary metabolite degradation',
    'Alcohol-Degradation'='Alcohol degradation','Fermentation'='Fermentation',
    'TCA-VARIANTS'='TCA variants','C1-COMPOUNDS'='One-carbon compounds',
    'Pentose-Phosphate-Cycle'='Pentose phosphate cycle')
  out<-unname(labels[x]);missing<-is.na(out)
  out[missing]<-gsub('-',' ',x[missing],fixed=TRUE);out
}
pfs_prepare_adult81 <- function(r,cfg,out) {
  all<-r$effects[r$effects$model=='B2'&r$effects$generation=='Overall',,drop=FALSE]
  selected<-all[tf_bool(all$analysis_valid)&tf_bool(all$FDR_lt_005),,drop=FALSE]
  stopifnot(nrow(selected)==81L,!anyDuplicated(selected$feature_id),all(selected$P_BH<.05))
  lineages<-r$classification$lineages
  priority<-c('Biosynthesis','Degradation','Energy-Metabolism','Activation-Inactivation-Interconversion',
    'Detoxification','Metabolic-Clusters','Super-Pathways','Unmapped')
  annotation<-ov_bind(lapply(sort(selected$feature_id),function(id){
    z<-lineages[lineages$feature_id==id&!is.na(lineages$lineage),,drop=FALSE]
    if(nrow(z)){
      rank<-match(z$level_1,priority);rank[is.na(rank)]<-length(priority)
      z<-z[order(rank,z$level_1,z$level_2,z$lineage),,drop=FALSE]
      canonical<-z$lineage[1];parts<-strsplit(canonical,'|',fixed=TRUE)[[1]]
      major<-parts[1];subclass<-if(length(parts)>=2)parts[2]else parts[1]
      # The deepest parent shared in the canonical lineage is a functional
      # annotation, never inferred from a disease coefficient or name keyword.
      deeper<-if(length(parts)>=3)paste(parts[3:length(parts)],collapse='|')else subclass
    }else{canonical<-'';major<-subclass<-deeper<-'Unmapped'}
    data.frame(feature_id=id,display_major=major,display_subclass=subclass,
      display_deeper_class=deeper,display_lineage=canonical,
      all_true_lineages=paste(sort(unique(z$lineage)),collapse='; '),
      multiple_lineages=nrow(z)>1L,mapped=nrow(z)>0L,stringsAsFactors=FALSE)
  }))
  annotation$major_rank<-match(annotation$display_major,priority)
  annotation$major_rank[is.na(annotation$major_rank)]<-length(priority)-1L
  annotation$major_rank[annotation$display_major=='Unmapped']<-length(priority)+1L
  annotation<-annotation[order(annotation$major_rank,annotation$display_major,
    annotation$display_subclass,annotation$display_deeper_class,annotation$feature_id),,drop=FALSE]
  z<-pfs_series(r,annotation$feature_id)
  profile<-matrix(NA_real_,nrow(annotation),6L,dimnames=list(annotation$feature_id,levels(z$series)))
  for(j in seq_len(nrow(z)))profile[as.character(z$feature_id[j]),as.character(z$series[j])]<-z$fill[j]
  key<-paste(annotation$display_major,annotation$display_subclass,annotation$display_deeper_class,sep=' / ')
  ordered<-unlist(lapply(unique(key),function(k){
    ids<-annotation$feature_id[key==k];a<-profile[ids,,drop=FALSE]
    usable<-colSums(is.finite(a))==nrow(a)
    if(length(ids)>=3L&&any(usable)&&all(annotation$mapped[key==k])){
      # Display only: no row standardisation, significance-based distance,
      # zero imputation or clustering across different functional subclasses.
      hc<-stats::hclust(stats::dist(a[,usable,drop=FALSE]),method='complete')
      ids[hc$order]
    }else if(!all(annotation$mapped[key==k])){
      ids[order(selected$pathway_name[match(ids,selected$feature_id)],ids)]
    }else ids
  }),use.names=FALSE)
  annotation<-annotation[match(ordered,annotation$feature_id),,drop=FALSE]
  groups<-paste(annotation$display_major,annotation$display_subclass,sep='|')
  group_order<-unique(groups);coordinates<-headers<-major_headers<-list();cursor<-0
  previous_major<-NULL
  for(k in seq_along(group_order)){
    ix<-which(groups==group_order[k]);n<-length(ix)
    a<-annotation[ix[1],]
    # One compact heading per major class; subclass boundaries retain their
    # original order without reserving an extra text row for every subclass.
    new_major<-!identical(previous_major,a$display_major)
    if(new_major){
      major_headers[[length(major_headers)+1L]]<-data.frame(
        y=-cursor-.2,label=pfs_class_text(a$display_major),major=a$display_major)
      cursor<-cursor+1.4
    }
    row_y<- -cursor-seq_len(n)+1;head_y<-row_y[1]+.5
    coordinates[[k]]<-data.frame(feature_id=annotation$feature_id[ix],y=row_y,group_index=k)
    label<-if(a$display_major==a$display_subclass)pfs_class_text(a$display_major)else
      paste(pfs_class_text(a$display_major),pfs_class_text(a$display_subclass),sep=' / ')
    headers[[k]]<-data.frame(group_index=k,y=head_y,label=label,
      major=a$display_major,subclass=a$display_subclass,major_start=new_major)
    cursor<-cursor+n;previous_major<-a$display_major
  }
  coord<-ov_bind(coordinates);headers<-ov_bind(headers)
  k<-match(annotation$feature_id,coord$feature_id)
  annotation$y<-coord$y[k];annotation$group_index<-coord$group_index[k]
  annotation$display_rank<-seq_len(nrow(annotation))
  s<-selected[match(annotation$feature_id,selected$feature_id),,drop=FALSE]
  for(v in setdiff(names(annotation),'feature_id'))s[[v]]<-annotation[[v]]
  s$selection_rule<-'All original valid adult Model 4 overall BH q<0.05 pathways; no new tests or subset BH'
  s$order_rule<-'Canonical verified lineage (major/subclass/deeper class); complete-link Euclidean clustering of original available six log2 effects within the deepest class; unmapped last'
  effects<-pfs_series(r,s$feature_id)
  effects$y<-s$y[match(as.character(effects$feature_id),s$feature_id)]
  # Use the original descriptive source sample and estimator. Only the pathway
  # set, pooled source ranking and row order change; no HDF5 import or refit.
  g<-r$sources$genus
  pool<-g[g$group=='All'&g$feature_id%in%s$feature_id&g$genus!='Unclassified',,drop=FALSE]
  score<-stats::aggregate(source_share_percent~genus,pool,sum)
  score<-score[order(-score$source_share_percent,score$genus),,drop=FALSE]
  top<-head(score$genus,8L);g<-g[g$feature_id%in%s$feature_id,,drop=FALSE]
  g$display_genus<-ifelse(g$genus=='Unclassified','Unclassified',
    ifelse(g$genus%in%top,g$genus,'Other named genera'))
  source<-stats::aggregate(mean_RA_percent~feature_id+group+display_genus,g,sum)
  total<-r$sources$totals
  source<-merge(source,total,by=c('feature_id','group'),sort=FALSE,all.x=TRUE)
  grid<-expand.grid(feature_id=s$feature_id,group=c('All','No_steatosis','Steatosis'),
    display_genus=c(top,'Other named genera','Unclassified'),stringsAsFactors=FALSE)
  source<-merge(grid,source,by=c('feature_id','group','display_genus'),sort=FALSE,all.x=TRUE)
  at<-match(paste(source$feature_id,source$group),paste(total$feature_id,total$group))
  denominator<-total$summed_stratified_mean_RA_percent[at]
  source$mean_RA_percent[is.na(source$mean_RA_percent)&denominator>0]<-0
  source$source_share_percent<-ifelse(denominator>0,100*source$mean_RA_percent/denominator,NA_real_)
  source$y<-s$y[match(source$feature_id,s$feature_id)]
  source$descriptive_only<-TRUE;source$source_sample_N<-r$sources$samples
  values<-source$source_share_percent[is.finite(source$source_share_percent)]
  fl_write_csv(data.frame(raw_minimum=format(min(values),digits=22),raw_maximum=format(max(values),digits=22),
    roundoff_excursions=sum(values<0|values>100),maximum_permitted_roundoff_percentage_points=1e-10,
    rule='Original numerical estimates retained; only <=1e-10 percentage-point floating roundoff may snap to the mathematical 0..100 colour domain'),
    file.path(out,'audit/source_colour_domain_roundoff.csv'))
  r$adult81_display<-list(revision=pfs_adult81_revision,selection=s,effects=effects,source=source,
    headers=headers,major_headers=ov_bind(major_headers),top_genera=top,y_limits=c(min(s$y)-.6,.9),
    source_sample_N=r$sources$samples,source_groups=r$sources$groups)
  ade_table(s,'Main_figure_representatives_EN',out)
  ade_table(s,'Figure44_Adult_Model4_all81_selection_order_EN',out)
  ade_table(headers,'Figure44_Functional_group_headers_EN',out)
  ade_table(ov_bind(major_headers),'Figure44_Visible_major_class_headers_EN',out)
  ade_table(effects,'Figure44b_Adult_Model4_all81_effects_EN',out)
  ade_table(source,'Figure44c_observed_source_matrix_EN',out)
  ade_table(source,'Main_figure_descriptive_source_shares_EN',out)
  ade_table(data.frame(rank=seq_along(top),genus=top,
    rule='Pooled source shares across all 81 adult Model-4 discoveries; no disease-group contrast used for genus ranking'),
    'Displayed_source_genera_EN',out)
  # Retire obsolete ten-row source aliases in this NEW run only. The previous
  # adopted run and its complete original files are retained under history.
  obsolete<-file.path(out,'tables',paste0('Figure44b_representative_adjusted_effects_EN',c('.csv','.html')))
  unlink(obsolete[file.exists(obsolete)],recursive=FALSE)
  r
}
pfs_adult81_panels <- function(r,cfg,ids=r$adult81_display$selection$feature_id) {
  d<-r$adult81_display;s<-d$selection[d$selection$feature_id%in%ids,,drop=FALSE]
  ids<-s$feature_id;effects<-d$effects[as.character(d$effects$feature_id)%in%ids,,drop=FALSE]
  source<-d$source[d$source$feature_id%in%ids&d$source$group!='All',,drop=FALSE]
  headers<-d$headers[d$headers$group_index%in%s$group_index,,drop=FALSE]
  stopifnot(identical(as.character(ids),as.character(d$selection$feature_id)))
  ylim<-d$y_limits
  lim<-pfs_effect_limit(d$effects);generas<-c(d$top_genera,'Other named genera','Unclassified')
  ng<-length(generas);source$x<-match(source$display_genus,generas)+ifelse(source$group=='Steatosis',ng+1,0)
  # Preserve exported source estimates. Only machine-roundoff excursions from
  # the mathematical 0..100 domain may be snapped for colour lookup; a real
  # out-of-range value blocks rendering rather than being clipped/hidden.
  observed<-source$source_share_percent[is.finite(source$source_share_percent)]
  stopifnot(all(observed>=-1e-10&observed<=100+1e-10))
  source$display_share_percent<-pmin(100,pmax(0,source$source_share_percent))
  base<-pfs_theme()+ggplot2::theme(axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),
    plot.title=ggplot2::element_blank(),
    axis.text.y=ggplot2::element_text(size=5.4),legend.position='bottom',
    plot.margin=ggplot2::margin(2,2,2,2,'pt'),
    legend.margin=ggplot2::margin(1,1,1,1,'pt'),legend.box.spacing=grid::unit(1,'mm'))
  b<-ggplot2::ggplot(effects,ggplot2::aes(series,y,fill=fill))+
    ggplot2::geom_hline(yintercept=headers$y,colour='#CCCCCC',linewidth=.18)+
    ggplot2::geom_tile(width=.96,height=.98,colour='white',linewidth=.12)+
    # Keep the complete rendered 5.1-pt glyph box inside a compact row. The
    # optical centre of Arial's asterisk differs from its typographic box.
    ggplot2::geom_text(ggplot2::aes(label=star),colour='black',size=5.1/(72.27/25.4),vjust=.4)+
    ggplot2::geom_text(data=d$major_headers,ggplot2::aes(x=.52,y=y,label=label),inherit.aes=FALSE,
      hjust=0,size=5.6/(72.27/25.4),family='Arial',fontface='bold')+
    ggplot2::scale_y_continuous(limits=ylim,breaks=s$y,labels=ids,expand=ggplot2::expansion(mult=0))+
    ggplot2::scale_x_discrete(drop=FALSE,labels=c('Model 3\nAll ages','Model 3\nAdults','Model 4\nAdults','G0','G1','G2'),
      expand=ggplot2::expansion(add=.55))+
    ggplot2::scale_fill_gradient2(low=unname(cfg$binary_colors[1]),mid='white',high=unname(cfg$binary_colors[2]),
      limits=c(-lim,lim),midpoint=0,na.value='#E5E5E5',breaks=c(-lim,0,lim),
      labels=ade_ratio_text(2^c(-lim,0,lim)),name='Abundance ratio',
      guide=ggplot2::guide_colourbar(direction='horizontal',title.position='top',display='rectangles',
        nbin=48,ticks=FALSE,barwidth=grid::unit(45,'mm'),barheight=grid::unit(2.8,'mm')))+
    ggplot2::labs(title=NULL,x=NULL,y=NULL)+base
  c<-ggplot2::ggplot(source,ggplot2::aes(x,y,fill=display_share_percent))+
    ggplot2::geom_hline(yintercept=headers$y,colour='#CCCCCC',linewidth=.18)+
    ggplot2::geom_tile(width=.94,height=.92,colour='white',linewidth=.12)+
    ggplot2::scale_y_continuous(limits=ylim,breaks=NULL,expand=ggplot2::expansion(mult=0))+
    ggplot2::scale_x_continuous(breaks=c(seq_len(ng),seq_len(ng)+ng+1),labels=rep(generas,2),
      expand=ggplot2::expansion(add=.55))+
    ggplot2::scale_fill_gradient(low='white',high='#397FC7',limits=c(0,100),breaks=c(0,50,100),
      name='Source share (%)',
      guide=ggplot2::guide_colourbar(direction='horizontal',title.position='top',display='rectangles',
        nbin=48,ticks=FALSE,barwidth=grid::unit(45,'mm'),barheight=grid::unit(2.8,'mm')))+
    ggplot2::annotate('text',x=c((ng+1)/2,(ng+1)/2+ng+1),y=ylim[2]-.35,
      label=c('No steatosis','Steatosis'),family='Arial',size=6/(72.27/25.4))+
    ggplot2::labs(title=NULL,x=NULL,y=NULL)+base+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=.5,size=5.4),
      axis.text.y=ggplot2::element_blank(),axis.title.y=ggplot2::element_blank())
  list(b=b,c=c,source=source,limits=ylim,colour_limit=lim)
}
