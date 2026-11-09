# Inert R-only definitions. Figure design adapted from Valles-Colomer et al.
# Nat Microbiol 2022 Fig. 2; reference tree and YFS abundance tracks, NOT pTR.
gp_protocol <- 'reference_phylogeny_generation_profiles_saved_P2_20261028_v1'
gp_pairs <- c('G0-G1','G1-G2','G0-G2')
gp_states <- c('Positive; BH q < 0.05','Negative; BH q < 0.05',
               'Tested; BH q >= 0.05','Not tested / unavailable')

gp_setup <- function(root) {
  .libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),file.path(R.home(),'library')))
  for (f in c('config.R','R/common.R','R/entry_helpers.R','R/nature_export_approved.R',
    'R/reanalysis_20261024.R','R/taxa_followup_figures.R','R/literature_ABC_20261028.R'))
    sys.source(file.path(root,f),envir=.GlobalEnv)
  for (p in c('ape','TreeSummarizedExperiment','ggplot2','patchwork','svglite','jsonlite'))
    if (!requireNamespace(p,quietly=TRUE)) stop('Missing dependency: ',p)
  Sys.setenv(PYTHONPATH=paste(c(file.path(root,'.figure_qa_python'),Sys.getenv('PYTHONPATH')),collapse=.Platform$path.sep))
  fl_config(root)
}

gp_pair_distances <- function(tr,ids) {
  # Small deterministic path audit without an N-tip squared distance matrix.
  nt <- length(tr$tip.label); parent <- integer(nt+tr$Nnode)
  parent[tr$edge[,2]] <- tr$edge[,1]; depth <- ape::node.depth.edgelength(tr)
  nodes <- match(ids,tr$tip.label); stopifnot(!anyNA(nodes))
  paths <- lapply(nodes,function(v) {a<-v;while(parent[v]!=0L){v<-parent[v];a<-c(a,v)};a})
  ans <- matrix(0,length(nodes),length(nodes),dimnames=list(ids,ids))
  for(i in seq_along(nodes))for(j in seq_len(i)) {
    lca <- paths[[i]][match(TRUE,paths[[i]]%in%paths[[j]])]
    ans[i,j] <- ans[j,i] <- depth[nodes[i]]+depth[nodes[j]]-2*depth[lca]
  }
  ans
}

gp_build <- function(cfg,out) {
  p24 <- lr_adopted(cfg$project_root,'24_Two_generation_reanalysis')
  p26 <- lr_adopted(cfg$project_root,'26_Taxon_family_sharing')
  stage <- file.path(p24,'private/stage_prepare.rds')
  tests <- file.path(p26,'tables/Taxon_cross_generation_GEE_all_results_EN.csv')
  originals <- c(cfg$microbiome_file,stage,tests,
    file.path(cfg$project_root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),
    file.path(cfg$project_root,'outputs/26_Taxon_family_sharing/LATEST_RUN.txt'))
  protected <- lr_hashes(originals)
  x <- readRDS(stage);m<-x$metadata;rel<-x$relative
  stopifnot(identical(x$protocol,v3_protocol),identical(colnames(rel),as.character(m$sample_id)),
    !anyDuplicated(m$sample_id),!anyDuplicated(rownames(rel)),all(is.finite(rel)),all(rel>=0),
    max(abs(colSums(rel)-1))<1e-10,all(rowSums(rel)>0),!anyNA(m$generation))
  fam <- table(m$family_id,factor(m$generation,levels=c('G0','G1','G2')))
  stopifnot(all(rowSums(fam>0)>=2L))
  tse <- suppressWarnings(readRDS(cfg$microbiome_file))
  tr <- TreeSummarizedExperiment::rowTree(tse)
  stopifnot(inherits(tr,'phylo'),!anyDuplicated(tr$tip.label),!is.null(tr$edge.length),
    all(is.finite(tr$edge.length)),all(tr$edge.length>=0))
  rd <- as.data.frame(SummarizedExperiment::rowData(tse))
  ids <- rownames(rel);at<-match(ids,rownames(rd))
  stopifnot(!anyNA(at),all(ids%in%tr$tip.label))
  tree <- ape::ladderize(ape::keep.tip(tr,ids),right=TRUE)
  probe <- sort(ids)[unique(round(seq(1,length(ids),length.out=min(40,length(ids)))))]
  delta <- max(abs(gp_pair_distances(tr,probe)-gp_pair_distances(tree,probe)))
  stopifnot(delta<1e-9)
  taxonomy <- data.frame(feature_id=ids,SGB_id=sub('^t__','',ids),
    kingdom=as.character(rd$kingdom[at]),phylum=sub('^p__','',as.character(rd$phylum[at])),
    species_annotation=as.character(rd$species[at]),stringsAsFactors=FALSE)
  taxonomy$phylum[is.na(taxonomy$phylum)|!nzchar(taxonomy$phylum)] <- 'Unclassified'
  taxonomy$overall_mean_RA_percent <- 100*rowMeans(rel)
  taxonomy$overall_prevalence_percent <- 100*rowMeans(rel>0)
  profiles <- do.call(rbind,lapply(c('G0','G1','G2'),function(g){
    z<-rel[,as.character(m$generation)==g,drop=FALSE]
    data.frame(feature_id=ids,generation=g,N=ncol(z),detected_n=rowSums(z>0),
      prevalence_percent=100*rowMeans(z>0),mean_RA_percent=100*rowMeans(z),stringsAsFactors=FALSE)
  }))
  stopifnot(all(profiles$prevalence_percent>=0 & profiles$prevalence_percent<=100),
    max(abs(tapply(profiles$mean_RA_percent,profiles$generation,sum)-100))<1e-8)
  t <- read.csv(tests,fileEncoding='UTF-8-BOM',check.names=FALSE,stringsAsFactors=FALSE)
  stopifnot(all(c('comparison','model','feature_id','eligible','analysis_valid','beta_log2_per_log2','P_raw','P_BH')%in%names(t)))
  t <- t[t$model=='P2',,drop=FALSE]
  stopifnot(!anyDuplicated(paste(t$comparison,t$feature_id)),all(t$feature_id%in%ids))
  valid <- t$eligible & t$analysis_valid & is.finite(t$P_raw) & is.finite(t$P_BH) & is.finite(t$beta_log2_per_log2)
  stopifnot(all(!t$eligible|valid),max(abs(p.adjust(t$P_raw[valid],'BH')-t$P_BH[valid]))<1e-10)
  marks <- expand.grid(feature_id=ids,comparison=gp_pairs,stringsAsFactors=FALSE)
  k <- match(paste(marks$comparison,marks$feature_id),paste(t$comparison,t$feature_id))
  marks$eligible <- !is.na(k) & t$eligible[k]
  marks$analysis_valid <- !is.na(k) & t$analysis_valid[k]
  for(v in c('beta_log2_per_log2','CI_low','CI_high','P_raw','P_BH',
    'younger_ratio_per_older_doubling','ratio_CI_low','ratio_CI_high'))marks[[v]]<-t[[v]][k]
  marks$state <- gp_states[4]
  ok <- marks$eligible & marks$analysis_valid & is.finite(marks$P_BH)
  marks$state[ok] <- gp_states[3]
  marks$state[ok & marks$P_BH<.05 & marks$beta_log2_per_log2>0] <- gp_states[1]
  marks$state[ok & marks$P_BH<.05 & marks$beta_log2_per_log2<0] <- gp_states[2]
  # Unavailable tests are not set to beta=0, P=1, or non-significant.
  stopifnot(all(is.na(marks$P_BH[!marks$eligible])))
  detail <- head(taxonomy$feature_id[order(-taxonomy$overall_mean_RA_percent,taxonomy$feature_id)],30L)
  detail_tree <- ape::ladderize(ape::keep.tip(tree,detail),right=TRUE)
  stopifnot(max(abs(gp_pair_distances(tree,detail)-gp_pair_distances(detail_tree,detail)))<1e-9)
  positive <- profiles$mean_RA_percent[profiles$mean_RA_percent>0]
  ra_limits <- c(floor(log10(min(positive))),ceiling(log10(max(positive))))
  counts <- do.call(rbind,lapply(c('G0','G1','G2'),function(g){
    z<-m[as.character(m$generation)==g,,drop=FALSE]
    data.frame(generation=g,participants=nrow(z),families=length(unique(z$family_id)),
      minors=sum(z$age<18,na.rm=TRUE),age_missing=sum(!is.finite(z$age)))
  }))
  sm <- do.call(rbind,lapply(gp_pairs,function(p){
    z<-marks[marks$comparison==p,,drop=FALSE]
    data.frame(comparison=p,SGBs_in_reference=nrow(z),eligible=sum(z$eligible),
      positive_BH_q_lt_005=sum(z$state==gp_states[1]),negative_BH_q_lt_005=sum(z$state==gp_states[2]),
      tested_not_significant=sum(z$state==gp_states[3]),not_tested=sum(z$state==gp_states[4]),BH_family_tests=sum(valid))
  }))
  taxcounts <- as.data.frame(table(taxonomy$phylum),stringsAsFactors=FALSE);names(taxcounts)<-c('phylum','SGBs')
  taxcounts <- taxcounts[order(-taxcounts$SGBs),]
  top_phyla <- head(taxcounts$phylum,4)
  phylum_colours <- setNames(c('#040676','#397FC7','#F1B656','#79CAFB','#A8A8A8'),c(top_phyla,'Other phyla'))
  taxonomy$branch_colour_group <- ifelse(taxonomy$phylum%in%top_phyla,taxonomy$phylum,'Other phyla')
  producer <- 'E:/YFS/documents-export-2025-07-07/documents-export-2025-07-07/make_MAE.R'
  provenance <- data.frame(item=c('Source tree tips','Retained baseline tips','Unmapped baseline SGBs',
    'Samples','Families','Reference tree origin','Exact reference version / inference algorithm',
    'Earlier export script clue','Distance preservation audit','New sequence inference','Statistical model refits'),
    value=c(length(tr$tip.label),length(tree$tip.label),0,ncol(rel),length(unique(m$family_id)),
      'rowTree in the original 2025 TreeSummarizedExperiment RDS',
      'Not recorded in the current RDS; not independently verified',
      'Earlier make_MAE.R imports mpa_vJun23_CHOCOPhlAnSGB_202403.nwk; exact provenance of current 2025 RDS remains unverified',
      format(delta,scientific=TRUE),'None','None'))
  r <- list(protocol=gp_protocol,tree=tree,detail_tree=detail_tree,taxonomy=taxonomy,profiles=profiles,
    concordance=marks,counts=counts,summary=sm,phylum_counts=taxcounts,
    phylum_colours=phylum_colours,ra_limits=ra_limits,detail_ids=detail,provenance=provenance,
    protected=protected,baseline_N=ncol(rel),families=length(unique(m$family_id)))
  gp_tables(r,out)
  ape::write.tree(tree,file=file.path(out,'tables/Reference_tree_all_SGBs.nwk'))
  fl_write_csv(protected,file.path(out,'audit/protected_source_files.csv'))
  if(file.exists(producer))writeLines(readLines(producer,warn=FALSE)[75:95],file.path(out,'audit/earlier_export_script_excerpt.txt'))
  lr_check_hashes(protected);r
}

gp_tables <- function(r,out) {
  for(n in c('counts','summary','provenance','phylum_counts'))
    v3_table(r[[n]],paste0('Reference_phylogeny_',n,'_EN'),out)
  v3_table(r$taxonomy,'All_SGB_annotations_and_overall_abundance_EN',out)
  v3_table(r$profiles,'All_SGB_generation_profiles_EN',out)
  v3_table(r$concordance,'All_SGB_saved_P2_concordance_EN',out)
  d<-r$taxonomy[match(r$detail_ids,r$taxonomy$feature_id),];d$display_rank<-seq_len(nrow(d))
  d$selection_rule<-'First 30 by overall mean relative abundance (all 3217 baseline participants); ID breaks ties'
  v3_table(d,'Detail_30_SGB_display_selection_EN',out)
}

gp_coordinates <- function(tr,taxonomy,circular=TRUE,gap_degrees=48) {
  tr<-ape::reorder.phylo(tr,'postorder');nt<-length(tr$tip.label);nn<-nt+tr$Nnode
  parent<-integer(nn);parent[tr$edge[,2]]<-tr$edge[,1]
  children<-split(tr$edge[,2],tr$edge[,1]);root<-setdiff(tr$edge[,1],tr$edge[,2])[1]
  # Depth-first planar order, not numerical tip IDs (which may be reordered).
  walk<-function(v)if(v<=nt)v else unlist(lapply(children[[as.character(v)]],walk),use.names=FALSE)
  ordertips<-walk(root);pos<-numeric(nn);pos[ordertips]<-seq_len(nt)
  group<-rep(NA_character_,nn)
  group[seq_len(nt)]<-taxonomy$branch_colour_group[match(tr$tip.label,taxonomy$feature_id)]
  descend<-function(v){if(v<=nt)return(invisible(NULL));ch<-children[[as.character(v)]]
    for(a in ch)descend(a)
    pos[v]<<-mean(pos[ch]);groups<-unique(group[ch]);group[v]<<-if(length(groups)==1L)groups else 'Other phyla'
  };descend(root)
  depth<-ape::node.depth.edgelength(tr);scale<-.96/max(depth);rad<-depth*scale
  xy<-function(rr,aa)cbind(x=rr*cos(aa),y=rr*sin(aa))
  if(circular){
    stopifnot(gap_degrees>=0,gap_degrees<360)
    gap<-gap_degrees*pi/180;step<-(2*pi-gap)/nt
    theta<-pi/2+gap/2+(pos-.5)*step
    begin<-xy(rad[tr$edge[,1]],theta[tr$edge[,2]])
    end<-xy(rad[tr$edge[,2]],theta[tr$edge[,2]])
    edges<-data.frame(x=begin[,1],y=begin[,2],xend=end[,1],yend=end[,2],group=group[tr$edge[,2]])
    arcs<-lapply(as.integer(names(children)),function(v){
      a<-range(theta[children[[as.character(v)]]]);n<-max(3,ceiling(diff(a)/(.6*pi/180)))
      q<-xy(rad[v],seq(a[1],a[2],length.out=n));data.frame(x=q[,1],y=q[,2],node=v,group=group[v])
    })
    tips<-data.frame(feature_id=tr$tip.label,theta=theta[seq_len(nt)],radius=rad[seq_len(nt)],
      leaf_order=pos[seq_len(nt)],step=step,group=group[seq_len(nt)])
    list(edges=edges,arcs=do.call(rbind,arcs),tips=tips,scale=scale)
  }else{
    y<-nt+1-pos
    edges<-data.frame(x=rad[tr$edge[,1]],y=y[tr$edge[,2]],xend=rad[tr$edge[,2]],yend=y[tr$edge[,2]],group=group[tr$edge[,2]])
    arcs<-do.call(rbind,lapply(as.integer(names(children)),function(v){
      yy<-range(y[children[[as.character(v)]]]);data.frame(x=rad[v],y=yy[1],xend=rad[v],yend=yy[2],group=group[v])
    }))
    tips<-data.frame(feature_id=tr$tip.label,x=rad[seq_len(nt)],y=y[seq_len(nt)],group=group[seq_len(nt)])
    list(edges=edges,arcs=arcs,tips=tips,scale=scale)
  }
}

gp_state_colours <- function(cfg)setNames(c(cfg$binary_colors['Comparison'],cfg$binary_colors['Reference'],
  '#BDBDBD','#F3F3F3'),gp_states)
gp_tint <- function(value,colour) {
  stopifnot(all(is.finite(value)),all(value>=0 & value<=1))
  grDevices::rgb(grDevices::colorRamp(c('#FFFFFF',colour))(value)/255)
}
gp_ra_fraction <- function(v,limits) {
  out<-rep(0,length(v));j<-v>0
  out[j]<-(log10(v[j])-limits[1])/diff(limits)
  stopifnot(all(is.finite(out)),all(out>=-1e-10 & out<=1+1e-10));pmin(1,pmax(0,out))
}
gp_tracks <- function(r,tips,cfg) {
  tracks<-do.call(rbind,lapply(seq_len(9),function(i){
    z<-tips;z$track<-i
    if(i<=6){
      g<-c('G0','G1','G2')[(i-1L)%%3+1L];p<-r$profiles[r$profiles$generation==g,]
      at<-match(z$feature_id,p$feature_id);stopifnot(!anyNA(at))
      z$value<-if(i<=3)p$prevalence_percent[at]else p$mean_RA_percent[at]
      v<-if(i<=3)z$value/100 else gp_ra_fraction(z$value,r$ra_limits)
      z$fill<-gp_tint(v,cfg$generation_colors[g]);z$label<-paste(if(i<=3)'P'else'A',g)
    }else{
      p<-r$concordance[r$concordance$comparison==gp_pairs[i-6],];at<-match(z$feature_id,p$feature_id)
      stopifnot(!anyNA(at));z$value<-p$P_BH[at];z$fill<-unname(gp_state_colours(cfg)[p$state[at]])
      z$label<-gp_pairs[i-6]
    };z
  }));tracks
}

gp_annuli <- function(tracks) {
  radii<-c(1.045,1.113,1.181,1.265,1.333,1.401,1.485,1.553,1.621)
  out<-lapply(seq_len(nrow(tracks)),function(i){
    z<-tracks[i,];a<-z$theta+c(-.5,.5)*z$step
    rr<-c(radii[z$track],radii[z$track]+.060)
    aa<-c(a[1],mean(a),a[2],a[2],mean(a),a[1])
    rad<-c(rep(rr[2],3),rep(rr[1],3))
    data.frame(x=rad*cos(aa),y=rad*sin(aa),cell=i,fill=z$fill)
  });do.call(rbind,out)
}

gp_text <- function(p,x,y,label,size=6,hjust=.5,...)p+ggplot2::annotate('text',x=x,y=y,label=label,
  size=size/(72.27/25.4),family='Arial',hjust=hjust,...)
gp_rect <- function(p,x0,x1,y0,y1,fill)p+ggplot2::annotate('rect',xmin=x0,xmax=x1,ymin=y0,ymax=y1,fill=fill,colour=NA)

gp_legends <- function(p,r,cfg,circle=TRUE) {
  # Explicit vector keys: never raster colourbars under/behind tick labels.
  if(circle){sx<-1;ox<-0;sy<-1;oy<-0}else{sx<-.94;ox<-1.66;sy<-11;oy<-19}
  X<-function(x)ox+sx*x;Y<-function(y)oy+sy*y
  txt<-function(x,y,s,size=5.7,h=.5){p<<-gp_text(p,X(x),Y(y),s,size,h)}
  rect<-function(x0,x1,y0,y1,col){p<<-gp_rect(p,X(x0),X(x1),Y(y0),Y(y1),col)}
  txt(0,-1.88,sprintf('%s SGBs; %s participants; stored reference tree',
    format(if(circle)nrow(r$taxonomy)else 30,big.mark=','),format(r$baseline_N,big.mark=',')),6.2)
  centers<-c(-.55,.35,1.25)
  for(j in 1:3){g<-c('G0','G1','G2')[j]
    txt(centers[j],-2.03,sprintf('%s (n = %s)',g,format(r$counts$participants[r$counts$generation==g],big.mark=',')),5.8)
    rect(centers[j]-.32,centers[j]+.32,-2.155,-2.105,unname(cfg$generation_colors[g]))
    for(a in seq(0,1,length.out=31)){
      rect(centers[j]-.32+a*.64,centers[j]-.32+min(1,a+1/30)*.64,-2.155,-2.105,gp_tint(a,cfg$generation_colors[g]))
      rect(centers[j]-.32+a*.64,centers[j]-.32+min(1,a+1/30)*.64,-2.375,-2.325,gp_tint(a,cfg$generation_colors[g]))
    }
    for(a in c(0,.5,1))txt(centers[j]-.32+a*.64,-2.205,as.character(100*a),5.2)
    breaks<-unique(c(r$ra_limits[1],round(mean(r$ra_limits)),r$ra_limits[2]))
    for(a in breaks)txt(centers[j]-.32+(a-r$ra_limits[1])/diff(r$ra_limits)*.64,-2.425,
      if(a>=-2 & a<=1)format(10^a,scientific=FALSE,trim=TRUE)else paste0('1e',a),5.2)
  }
  txt(-1.76,-2.13,'P: detection (%)',5.7,0)
  txt(-1.76,-2.35,'A: mean RA (%)',5.7,0)
  txt(-1.76,-2.45,'Log colour; zero = white',5.3,0)
  # Four clearly separated states; light grey is explicitly unavailable.
  entries<-c('Positive q < 0.05','Negative q < 0.05','Tested, q >= 0.05','Not tested')
  xx<-c(-1.70,-.80,.12,1.10)
  for(i in 1:4){rect(xx[i],xx[i]+.065,-2.62,-2.555,unname(gp_state_colours(cfg)[i]));txt(xx[i]+.10,-2.587,entries[i],5.3,0)}
  cols<-r$phylum_colours
  for(i in seq_along(cols)){
    xx<-if(i<=3)c(-1.70,-.47,.77)[i]else c(-1.15,.40)[i-3]
    yy<-if(i<=3)-2.77 else -2.92
    rect(xx,xx+.065,yy-.025,yy+.025,cols[i]);txt(xx+.10,yy,names(cols)[i],5.5,0)
  }
  p
}

gp_circular <- function(r,cfg) {
  g<-gp_coordinates(r$tree,r$taxonomy,TRUE);tracks<-gp_tracks(r,g$tips,cfg);rings<-gp_annuli(tracks)
  # Dotted terminal guides are not inferred branches; they align the tracks.
  tips<-g$tips;tips$x<-tips$radius*cos(tips$theta);tips$y<-tips$radius*sin(tips$theta)
  tips$xend<-.992*cos(tips$theta);tips$yend<-.992*sin(tips$theta)
  p<-ggplot2::ggplot()+
    ggplot2::geom_segment(data=tips,ggplot2::aes(x,y,xend=xend,yend=yend),colour='#E5E5E5',linewidth=.06,linetype=3)+
    ggplot2::geom_segment(data=g$edges,ggplot2::aes(x,y,xend=xend,yend=yend,colour=group),linewidth=.11)+
    ggplot2::geom_path(data=g$arcs,ggplot2::aes(x,y,group=node,colour=group),linewidth=.11)+
    ggplot2::scale_colour_manual(values=r$phylum_colours,guide='none')+
    ggplot2::geom_polygon(data=rings,ggplot2::aes(x,y,group=cell,fill=fill),colour=NA)+
    ggplot2::scale_fill_identity()+ggplot2::coord_equal(xlim=c(-1.86,1.86),ylim=c(-3.05,1.78),clip='off')+
    ggplot2::labs(title='Reference phylogeny and generation-specific microbiome profiles',x=NULL,y=NULL)+
    ggplot2::theme_void(base_family='Arial')+ggplot2::theme(plot.margin=ggplot2::margin(7,5,5,5,'pt'))
  yy<-c(1.075,1.143,1.211,1.295,1.363,1.431,1.515,1.583,1.651)
  for(i in 1:9)p<-gp_text(p,0,yy[i],unique(tracks$label[tracks$track==i]),5.5)
  # The bar length is in the stored tree's unverified branch-length units.
  p<-p+ggplot2::annotate('segment',x=-.30,xend=-.30+.2*g$scale,y=-1.78,yend=-1.78,linewidth=.3)
  p<-gp_text(p,.13,-1.78,'0.2 stored branch-length units',5.5,0)
  gp_legends(p,r,cfg,TRUE)
}

gp_detail <- function(r,cfg) {
  g<-gp_coordinates(r$detail_tree,r$taxonomy,FALSE);tracks<-gp_tracks(r,g$tips,cfg)
  d<-g$tips;d$label<-sub('^t__','',d$feature_id)
  tracks$x<-1.92+(tracks$track-1)*.16
  p<-ggplot2::ggplot()+
    ggplot2::geom_segment(data=g$edges,ggplot2::aes(x,y,xend=xend,yend=yend,colour=group),linewidth=.3)+
    ggplot2::geom_segment(data=g$arcs,ggplot2::aes(x,y,xend=xend,yend=yend,colour=group),linewidth=.3)+
    ggplot2::geom_segment(data=d,ggplot2::aes(x,y,xend=1.08,yend=y),colour='#D0D0D0',linetype=3,linewidth=.2)+
    ggplot2::geom_text(data=d,ggplot2::aes(x=1.12,y=y,label=label),hjust=0,size=6/(72.27/25.4),family='Arial')+
    ggplot2::geom_tile(data=tracks,ggplot2::aes(x,y,fill=fill),width=.148,height=.85)+
    ggplot2::scale_fill_identity()+ggplot2::scale_colour_manual(values=r$phylum_colours,guide='none')+
    ggplot2::coord_cartesian(xlim=c(-.06,3.43),ylim=c(-15,34),clip='off')+
    ggplot2::theme_void(base_family='Arial')+
    ggplot2::labs(title='Thirty most abundant SGBs: reference phylogeny and generation profiles')+
    ggplot2::theme(plot.margin=ggplot2::margin(7,5,5,5,'pt'))
  for(i in 1:9)p<-gp_text(p,1.92+(i-1)*.16,31.4,
    if(i<=6)c('G0','G1','G2')[(i-1)%%3+1]else gp_pairs[i-6],5.5,angle=90)
  for(i in 1:3)p<-gp_text(p,c(2.08,2.56,3.04)[i],33.0,c('Detection (%)','Mean RA (%)','P2 association')[i],5.7)
  # Separate compact vector legends below the data. Same meanings and scales.
  p<-gp_legends(p,r,cfg,FALSE)
  p
}

gp_reports <- function(r,out) {
  writeLines(c('# Figure 35. Reference phylogeny and generation-specific microbiome profiles.',
    sprintf('All %s SGBs in the adopted microbiome baseline (%s participants: G0, %d; G1, %d; G2, %d) are shown.',
      format(nrow(r$taxonomy),big.mark=','),format(r$baseline_N,big.mark=','),r$counts$participants[1],r$counts$participants[2],r$counts$participants[3]),
    'The backbone was pruned from the reference phylogeny stored in the original MetaPhlAn TreeSummarizedExperiment; branch lengths were retained. Its exact database version and inference algorithm are not verified. Dotted terminal guides only align tips with the annotation rings.',
    'From inside out, rings show detection prevalence (P; G0/G1/G2), mean relative abundance (A; G0/G1/G2), and the saved adjusted older-to-younger abundance associations (P2; G0-G1/G1-G2/G0-G2). Detection is relative abundance >0. Mean relative abundance includes zeros. P uses a common 0-100% linear intensity scale; A uses the same log10-positive-mean colour scale for all generations, with zero shown in white. There is no new abundance transformation or statistical model fit.',
    'Orange and blue denote positive and negative associations with BH q<0.05, respectively; medium grey indicates an eligible non-significant test, and light grey an unavailable test. BH values are copied without modification from module 26, pooling 1,512 eligible SGB-by-comparison tests for P2.',
    'P2 Gaussian GEE models account for primary-family clustering and both endpoints\' age splines, sex, library batch, sequencing pool, sampling season and recent antibiotic use. A primary-family label does not verify a direct pedigree or shared household. The associations are abundance concordance, not evidence of strain transmission or a potential transmission rate.',
    '', '# Supplementary Figure 35. Detail view of the 30 most abundant SGBs.',
    'SGBs are selected only by their overall mean relative abundance across all baseline participants (ties by ID), independently of association P values. The reference subtree preserves pairwise patristic distances; table rows follow that tree. Track definitions and colour scales match Figure 35. Full annotations for all SGBs are provided in the source tables.'),
    file.path(out,'reports/Figure_legends_EN.md'),useBytes=TRUE)
  writeLines(c('# Methods and provenance',
    'This is a descriptive figure adaptation, not a replication of the original strain-level experiment. No sequences were inferred, no maximum-likelihood tree was newly estimated, and no genotype counts, nGDs or pTRs were calculated.',
    'Source: Valles-Colomer M, Bacigalupe R, Vieira-Silva S, et al. Variation and transmission of the human gut microbiota across multiple familial generations. Nature Microbiology. 2022;7:87-96. doi:10.1038/s41564-021-01021-8. Figure 2 supplies the reference-tree-plus-outer-tracks design.',
    'https://www.nature.com/articles/s41564-021-01021-8',
    'An earlier export script imports mpa_vJun23_CHOCOPhlAnSGB_202403.nwk. This is provenance evidence for an earlier export only; it does not establish the exact version of the tree in the current 2025 RDS. The current RDS provides no tree-building metadata or reference sequences.',
    'The adopted module-24 baseline is used without further covariate exclusions or filtering. Families need at least two represented generations. Minors are retained. Means include zeros, and no re-normalization is performed. The full reference plot excludes none of the baseline SGBs.',
    'Official SGB identifiers are retained; a species-level genome bin is not a sample-specific strain. Multiple species annotations are not resolved to invented single names. Source checksums and path-distance preservation checks are stored in the audit directory.',
    'There are no new hypothesis tests in this module, and consequently no additional FDR procedure. The pre-existing module-26 P2 q values retain their original test family.'),
    file.path(out,'reports/Methods_and_provenance_EN.md'),useBytes=TRUE)
  writeLines(c('# 35号：参考系统树与代际菌群概览',
    sprintf('采用现行两代可用菌群基线：%d人、%d个主家系、%d个SGB。G0/G1/G2人数分别为%d/%d/%d。',
      r$baseline_N,r$families,nrow(r$taxonomy),r$counts$participants[1],r$counts$participants[2],r$counts$participants[3]),
    '总览图保留全部SGB。检出率为相对丰度大于0的人数/该代总人数；平均相对丰度包括未检出的0值。丰度最高30个SGB另画细节图，选择与P值无关。',
    '树直接来自原始2025版RDS的rowTree；全体基线SGB均可映射。保留分支长度，但当前RDS未记录确切参考库版本和建树算法，因此不命名为本研究重新构建的最大似然树。',
    '三条最外环复用26号P2模型：橙色为正关联且BH q<0.05；蓝色为负关联且BH q<0.05；中灰色为已检验但未显著；浅灰色为不符合检出率门槛或无可用检验。',
    '这些是跨代丰度关联，不是菌株传播、遗传距离或传播率。没有新增回归、置换或协变量筛选。原始RDS、现行基线和26号结果均保持不变。',
    '重画：运行本脚本时添加 --figures-only；发生绘图失败时用 --resume-figures 读取已保存结果。'),
    file.path(out,'reports/README_zh.md'),useBytes=TRUE)
}

gp_preflight <- function(cfg,out) {
  fun<-c('gp_figures','gp_circular_theme','gp_circular','gp_detail','gp_coordinates','gp_annuli','gp_tracks','gp_legends','gp_text','gp_rect','gp_tint','gp_ra_fraction',
    'lr_export','tf_export','nf_theme','nf_skill_root')
  path<-file.path(out,'audit/actual_plot_and_export_source.R')
  writeLines(unlist(lapply(fun,function(n)c(paste0(n,' <- '),deparse(get(n,mode='function'),width.cutoff=120),''))),path)
  nf_qa_command(cfg,'validate_figure.py',c(shQuote(path),'--json'),file.path(out,'audit/plot_source_preflight.json'))
}

gp_figures <- function(r,cfg,out) {
  v3_log('Drawing all ',nrow(r$taxonomy),' reference tips and nine annotation tracks')
  p<-gp_circular(r,cfg)+lr_theme()+gp_circular_theme()+ggplot2::labs(tag='a')
  a<-tf_export(list(plot=patchwork::wrap_plots(list(p)),width=183,height=240,panels=1L,
    exemptions=list(),fixed_aspect_single_panel=TRUE),'Figure35_Generation_reference_phylogeny',out,cfg)
  v3_log('Drawing the outcome-blind top-30 detail view')
  b<-lr_export(list(gp_detail(r,cfg)),'Supplement35_Thirty_most_abundant_SGBs',out,cfg,width=183,height=210)
  v3_table(rbind(a,b),'Figure_export_manifest_EN',out)
}

gp_circular_theme <- function()ggplot2::theme_void(base_family='Arial')+ggplot2::theme(
  plot.title=ggplot2::element_text(size=6.8,margin=ggplot2::margin(l=12,b=6,unit='pt')),
  plot.tag=ggplot2::element_text(size=8,face='bold'),plot.tag.position=c(0,1),
  plot.margin=ggplot2::margin(7,5,5,5,'pt'))

gp_main <- function(args=character(),root=Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')) {
  if(getRversion()<'4.5'){
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',
      shQuote(file.path(root,'scripts/35_generation_reference_phylogeny.R')),vapply(args,shQuote,character(1))))
    if(status!=0L)stop('Project-R reference-phylogeny entry failed');return(invisible(NULL))
  }
  if(any(!args%in%c('--run','--check','--help','--figures-only','--resume-figures')))stop('Unknown mode')
  if(sum(c('--run','--check','--figures-only','--resume-figures')%in%args)>1L)stop('Choose one execution mode')
  if('--help'%in%args){cat('35_generation_reference_phylogeny.R [--run | --check | --figures-only | --resume-figures]\nNo statistical refits.\n');return(invisible(NULL))}
  cfg<-gp_setup(root);key<-'generation_phylogeny'
  cfg$output_module<-unname(cfg$output_modules[key]);cfg$run_id<-unname(cfg$run_ids[key]);cfg$run_id_explicit<-FALSE
  if('--check'%in%args){parse(file.path(root,'R/generation_reference_phylogeny.R'));cat('PASS: syntax and dependencies; no data loaded\n');return(invisible(NULL))}
  module<-file.path(root,'outputs',cfg$output_module)
  if('--figures-only'%in%args){
    out<-lr_adopted(root,cfg$output_module);r<-readRDS(file.path(out,'private/figure_data.rds'))
    lr_check_hashes(r$protected);gp_preflight(cfg,out);gp_figures(r,cfg,out);lr_check_hashes(r$protected)
    sys.source(file.path(root,'R/result_gallery.R'),envir=.GlobalEnv);fl_publish_organized_results(cfg);return(invisible(out))
  }
  resume<-NULL
  if('--resume-figures'%in%args){
    paths<-list.dirs(file.path(module,'history'),recursive=FALSE,full.names=TRUE)
    paths<-paths[file.exists(file.path(paths,'private/figure_data.rds'))]
    if(!length(paths))stop('No saved figure data found')
    resume<-paths[which.max(file.info(file.path(paths,'private/figure_data.rds'))$mtime)]
  }
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(p in c('tables','figures','reports','audit','private'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
  writeLines(c('RUNNING',gp_protocol),file.path(out,'run_status.txt'))
  tryCatch({
    if(is.null(resume))r<-gp_build(cfg,out)else{
      r<-readRDS(file.path(resume,'private/figure_data.rds'));lr_check_hashes(r$protected);gp_tables(r,out)
      ape::write.tree(r$tree,file=file.path(out,'tables/Reference_tree_all_SGBs.nwk'))
      writeLines(paste('Recovered saved descriptive data; no model fits:',resume),file.path(out,'audit/recovery_mode.txt'))
    }
    saveRDS(r,file.path(out,'private/figure_data.rds'));gp_reports(r,out);gp_preflight(cfg,out)
    gp_figures(r,cfg,out);lr_check_hashes(r$protected)
    fl_write_csv(r$protected,file.path(out,'audit/protected_source_files.csv'))
    writeLines(capture.output(sessionInfo()),file.path(out,'audit/sessionInfo.txt'))
    writeLines('PASS: source checksums unchanged; no new models; strict rendered figure gates passed.',file.path(out,'audit/delivery_validation.txt'))
    writeLines(c('COMPLETED',gp_protocol),file.path(out,'run_status.txt'))
    target<-fl_output_complete_run(ctx);writeLines(basename(target),file.path(module,'LATEST_RUN.txt'))
    v3_log('COMPLETED: ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
