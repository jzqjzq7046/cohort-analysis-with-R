# R-only circular reference-tree adaptation of Valles-Colomer et al. Figure 2.
# The saved species/SGB backbone is not a newly inferred study ML phylogeny.

sgr_selected_scope <- function(cfg,selected,reference) {
  p36 <- lr_adopted(cfg$project_root,'36_Genus_generation_taxonomy')
  files <- c(file.path(p36,'private/analysis_result.rds'),
    file.path(cfg$project_root,'outputs/36_Genus_generation_taxonomy/LATEST_RUN.txt'))
  protected <- lr_hashes(files);genus <- readRDS(files[1])
  stopifnot(identical(genus$protocol,gg_protocol),!anyDuplicated(selected$feature_id),
    all(tf_bool(selected$analysis_valid)),all(selected$P_BH<.05))
  m <- genus$mapping
  m <- m[tf_bool(m$included)&m$genus%in%selected$feature_id,,drop=FALSE]
  stopifnot(!anyDuplicated(m$feature_id),setequal(m$genus,selected$feature_id),
    all(m$feature_id%in%reference$tree$tip.label))
  tr <- reference$tree;nt <- length(tr$tip.label);nn <- nt+tr$Nnode
  g <- m$genus[match(tr$tip.label,m$feature_id)]
  homogeneous <- rep(NA_character_,nn);homogeneous[seq_len(nt)] <- g
  children <- split(tr$edge[,2],tr$edge[,1]);parent <- integer(nn)
  parent[tr$edge[,2]] <- tr$edge[,1];root <- setdiff(tr$edge[,1],tr$edge[,2])[1]
  visit <- function(v) {
    if(v<=nt)return(invisible(NULL))
    ch <- children[[as.character(v)]];for(w in ch)visit(w)
    a <- homogeneous[ch]
    if(!anyNA(a)&&length(unique(a))==1L)homogeneous[v]<<-a[1]
    invisible(NULL)
  }
  visit(root)
  nodes <- which(!is.na(homogeneous));maximal <- nodes[parent[nodes]==0L]
  other <- nodes[parent[nodes]!=0L]
  maximal <- c(maximal,other[is.na(homogeneous[parent[other]])])
  clades <- table(factor(homogeneous[maximal],levels=selected$feature_id))
  counts <- table(factor(m$genus,levels=selected$feature_id))
  a <- selected[,c('feature_id','genus','phylum','family','P_raw','P_BH','beta_log2',
    'abundance_ratio','ratio_CI_low','ratio_CI_high','overall_prevalence_percent','overall_mean_RA_percent')]
  a$mapped_SGBs <- as.integer(counts)
  a$maximal_homogeneous_reference_clades <- as.integer(clades)
  a$monophyletic_in_stored_reference <- a$maximal_homogeneous_reference_clades==1L
  stopifnot(all(a$mapped_SGBs>0L),all(a$maximal_homogeneous_reference_clades>0L),
    sum(a$mapped_SGBs)==nrow(m))
  list(selected=a,mapping=m,protected=protected,
    N_selected_genera=nrow(a),N_mapped_SGBs=nrow(m),
    N_nonmonophyletic_genera=sum(!a$monophyletic_in_stored_reference))
}

sgr_build <- function(cfg) {
  parent <- lr_adopted(cfg$project_root,'35_Generation_reference_phylogeny')
  files <- c(file.path(parent,'private/figure_data.rds'),
    file.path(cfg$project_root,'outputs/35_Generation_reference_phylogeny/LATEST_RUN.txt'))
  protected <- lr_hashes(files);r <- readRDS(files[1]);lr_check_hashes(r$protected)
  stopifnot(identical(r$protocol,gp_protocol),!anyDuplicated(r$tree$tip.label),
    setequal(r$tree$tip.label,r$taxonomy$feature_id),
    all(is.finite(r$tree$edge.length)),all(r$tree$edge.length>=0))
  # Harmonize exact phylum synonyms; do not collapse the genus analysis to a
  # representative tip or mistake the taxonomy hierarchy for branch lengths.
  tax <- r$taxonomy
  tax$phylum_original <- tax$phylum
  tax$phylum[tax$phylum=='Firmicutes'] <- 'Bacillota'
  tax$phylum[tax$phylum=='Actinobacteria'] <- 'Actinomycetota'
  tax$phylum[tax$phylum=='Bacteroidetes'] <- 'Bacteroidota'
  tax$phylum[tax$phylum=='Proteobacteria'] <- 'Pseudomonadota'
  tax$phylum[tax$phylum=='Verrucomicrobia'] <- 'Verrucomicrobiota'
  groups <- c('Bacillota','Bacteroidota','Actinomycetota','Pseudomonadota',
    'Verrucomicrobiota','Archaea','Other phyla')
  tax$branch_colour_group <- ifelse(tax$phylum%in%groups,tax$phylum,'Other phyla')
  tax$branch_colour_group[tax$kingdom%in%c('k__Archaea','Archaea')] <- 'Archaea'
  colours <- setNames(c('#040676','#397FC7','#F1B656','#79CAFB','#F7A24F',
    '#6A95A3','#B4B4B4'),groups)
  # Keep the same observed baseline denominator and include zero abundance.
  stopifnot(all(tax$overall_prevalence_percent>0 & tax$overall_prevalence_percent<=100),
    all(tax$overall_mean_RA_percent>0),abs(sum(tax$overall_mean_RA_percent)-100)<1e-8)
  logs <- log10(tax$overall_mean_RA_percent)
  lim <- c(floor(min(logs)),ceiling(max(logs)))
  if(diff(lim)==0)lim <- lim+c(-.5,.5)
  g <- gp_coordinates(r$tree,tax,TRUE)
  d <- tax[match(g$tips$feature_id,tax$feature_id),,drop=FALSE]
  d$theta <- g$tips$theta;d$radius <- g$tips$radius;d$step <- g$tips$step
  d$leaf_order <- g$tips$leaf_order
  d$prevalence_fill_fraction <- d$overall_prevalence_percent/100
  d$RA_fill_fraction <- (log10(d$overall_mean_RA_percent)-lim[1])/diff(lim)
  stopifnot(all(d$prevalence_fill_fraction>=0 & d$prevalence_fill_fraction<=1),
    all(d$RA_fill_fraction>=0 & d$RA_fill_fraction<=1))
  availability <- data.frame(component=c('Reference tree','Prevalence','Relative abundance',
    'StrainPhlAn-profiled species marker','Potential transmission rate (pTR)','Recovered genotype count'),
    availability=c(rep('Available',3),rep('Unavailable',3)),
    action=c('Show all baseline tips; preserve original branch lengths',
      'Linear 0-100% colour intensity; all baseline participants',
      'Mean percentage including zeros; common log10 colour intensity',
      'Not displayed; species/SGB detection does not establish genotyping',
      'Not displayed; no strain/genotype sharing measurements',
      'Not displayed; detected participants are not recovered genotypes'))
  counts <- as.data.frame(table(factor(d$branch_colour_group,levels=groups)),stringsAsFactors=FALSE)
  names(counts) <- c('branch_colour_group','SGBs')
  list(parent=parent,protected=protected,inherited_protected=r$protected,
    taxonomy=d,geometry=g,tree=r$tree,colours=colours,RA_limits=lim,
    availability=availability,counts=counts,N=r$baseline_N,families=r$families,
    provenance=r$provenance)
}

sgr_rings <- function(r,cfg) {
  x <- r$taxonomy;out <- vector('list',2L)
  for(k in 1:2) {
    lo <- c(1.035,1.165)[k];hi <- lo+.105
    a <- cbind(x$theta-x$step/2,x$theta,x$theta+x$step/2,
      x$theta+x$step/2,x$theta,x$theta-x$step/2)
    radius <- matrix(rep(c(hi,hi,hi,lo,lo,lo),each=nrow(x)),nrow(x),6)
    value <- if(k==1)x$prevalence_fill_fraction else x$RA_fill_fraction
    fill <- gp_tint(value,unname(cfg$binary_colors[c('Reference','Comparison')[k]]))
    out[[k]] <- data.frame(x=as.vector(t(radius*cos(a))),y=as.vector(t(radius*sin(a))),
      cell=rep(seq_len(nrow(x))+(k-1L)*nrow(x),each=6),fill=rep(fill,each=6))
  }
  ov_bind(out)
}

sgr_plot <- function(r,cfg) {
  g <- r$geometry;d <- r$taxonomy;rings <- sgr_rings(r,cfg)
  guide <- data.frame(x=d$radius*cos(d$theta),y=d$radius*sin(d$theta),
    xend=1.02*cos(d$theta),yend=1.02*sin(d$theta))
  p <- ggplot2::ggplot()+
    ggplot2::geom_segment(data=guide,ggplot2::aes(x,y,xend=xend,yend=yend),
      colour='#E6E6E6',linewidth=.045)+
    ggplot2::geom_segment(data=g$edges,ggplot2::aes(x,y,xend=xend,yend=yend,colour=group),linewidth=.095)+
    ggplot2::geom_path(data=g$arcs,ggplot2::aes(x,y,group=node,colour=group),linewidth=.095)+
    ggplot2::scale_colour_manual(values=r$colours,guide='none')+
    ggplot2::geom_polygon(data=rings,ggplot2::aes(x,y,group=cell,fill=fill),colour=NA)+
    ggplot2::scale_fill_identity()+
    ggplot2::coord_equal(xlim=c(-1.40,1.40),ylim=c(-2.65,1.35),expand=FALSE,clip='off')+
    ggplot2::labs(title='Cohort reference phylogeny',x=NULL,y=NULL)+gp_circular_theme()
  # The open angular sector labels observed tracks; no absent pTR/genotype track.
  p <- gp_text(p,0,1.087,'P',6)
  p <- gp_text(p,0,1.217,'A',6)
  p <- p+ggplot2::annotate('segment',x=-.32,xend=-.32+.2*g$scale,
    y=-1.37,yend=-1.37,linewidth=.3)
  p <- gp_text(p,.15,-1.37,'0.2 reference-tree units',5.5,0)
  p <- gp_text(p,0,-1.50,sprintf('%s SGBs; %s participants',
    format(nrow(d),big.mark=','),format(r$N,big.mark=',')),6.2)
  # Vector legends with spatially separated titles and tick labels.
  rects <- texts <- list()
  txt <- function(x,y,label,size=5.7,hjust=.5) {
    texts[[length(texts)+1L]] <<- data.frame(x=x,y=y,label=label,
      size=size/(72.27/25.4),hjust=hjust)
  }
  for(k in 1:2) {
    y <- c(-1.70,-1.98)[k];left <- -.44;right <- 1.24
    v <- seq(0,1,length.out=121)
    rects[[k]] <- data.frame(xmin=left+v*(right-left),
      xmax=left+pmin(1,v+1/120)*(right-left),ymin=y,ymax=y+.055,
      fill=gp_tint(v,unname(cfg$binary_colors[c('Reference','Comparison')[k]])))
    txt(-1.30,y+.027,if(k==1)'P: Prevalence (%)'else'A: Mean RA (%)',5.7,0)
    val <- if(k==1)c(0,50,100)else 10^c(r$RA_limits[1],mean(r$RA_limits),r$RA_limits[2])
    labs <- if(k==1)as.character(val)else format(signif(val,2),scientific=TRUE,trim=TRUE)
    for(j in 1:3)txt(left+c(0,.5,1)[j]*(right-left),y-.070,labs[j],5.3)
  }
  txt(-1.30,-2.18,'Phylum',5.8,0)
  for(i in seq_along(r$colours)) {
    x <- c(-1.30,-.36,.58)[(i-1)%%3+1];y <- -2.31-floor((i-1)/3)*.13
    rects[[length(rects)+1L]] <- data.frame(xmin=x,xmax=x+.07,
      ymin=y-.022,ymax=y+.022,fill=unname(r$colours[i]))
    txt(x+.10,y,names(r$colours)[i],5.5,0)
  }
  p+ggplot2::geom_rect(data=ov_bind(rects),ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill),
    colour=NA,inherit.aes=FALSE,show.legend=FALSE)+
    ggplot2::geom_text(data=ov_bind(texts),ggplot2::aes(x,y,label=label,size=size,hjust=hjust),
      family='Arial',inherit.aes=FALSE,show.legend=FALSE)+ggplot2::scale_size_identity()
}

sgr_tables <- function(r,out) {
  v3_table(r$taxonomy,'Figure39_Reference_phylogeny_all_SGB_source_EN',out)
  v3_table(r$counts,'Figure39_Reference_phylogeny_phylum_counts_EN',out)
  v3_table(r$availability,'Figure39_Reference_phylogeny_data_availability_EN',out)
  v3_table(r$provenance,'Figure39_Reference_phylogeny_provenance_EN',out)
  v3_table(data.frame(track=c('Prevalence','Mean relative abundance'),
    minimum=c(0,10^r$RA_limits[1]),maximum=c(100,10^r$RA_limits[2]),
    colour_intensity=c('Linear percentage','Log10 positive mean percentage')),
    'Figure39_Reference_phylogeny_colour_limits_EN',out)
  ape::write.tree(r$tree,file=file.path(out,'tables/Figure39_Reference_tree_all_SGBs.nwk'))
}

sgr_legend <- function(r,out) {
  writeLines(c('# Figure 39 | Cohort reference phylogeny and abundance profiles.',
    sprintf('All %s detected SGBs in the current microbiome baseline (%s participants from %s primary families) are shown. Branches follow the reference tree saved in the original microbiome object; original branch lengths are retained. The exact database version and tree-inference algorithm are unverified. This is not a newly inferred study-specific maximum-likelihood tree, nor an asserted match to the PhyloPhlAn default database.',
      format(nrow(r$taxonomy),big.mark=','),format(r$N,big.mark=','),format(r$families,big.mark=',')),
    'Branches are coloured by the five named major bacterial phyla, Archaea and other phyla. Firmicutes is displayed as Bacillota. Archaea is a domain-level display group, matching the reference paper\'s broad grouping; its constituent phyla remain explicit in the source table. Mixed ancestral branches use the neutral other-phyla colour.',
    'From inside out, the observed rings show prevalence (P, blue) and mean relative abundance (A, orange). Prevalence is the percentage of all baseline participants with relative abundance >0 and uses a linear 0-100% intensity scale. Mean relative abundance includes non-detection zeros in the denominator and uses a common log10-positive-mean intensity scale; colourbar ticks are actual percentages. Thin pale terminal guides align the non-ultrametric tree with the rings. Individual SGB labels are omitted at this density; every identifier, species annotation and value is provided in the source table and Newick file. No tips are selected by abundance or significance.',
    'The layout is adapted from Valles-Colomer M, Bacigalupe R, Vieira-Silva S, et al. Variation and transmission of the human gut microbiota across multiple familial generations. Nature Microbiology. 2022;7:87-96. Figure 2. doi:10.1038/s41564-021-01021-8. https://www.nature.com/articles/s41564-021-01021-8',
    'The source paper\'s StrainPhlAn-profiled species labels, potential transmission-rate (pTR) ring and recovered-genotype bar plots cannot be reproduced from the available abundance/reference-tree data. They are omitted rather than encoded as zero or replaced by participant detection or abundance concordance. No strain reconstruction, pTR estimation or statistical model fitting is performed. The tree is at species/SGB resolution; the separate significant-genus maps remain genus-level analyses.'),
    file.path(out,'reports/Figure39_Reference_phylogeny_legend_EN.md'),useBytes=TRUE)
}
