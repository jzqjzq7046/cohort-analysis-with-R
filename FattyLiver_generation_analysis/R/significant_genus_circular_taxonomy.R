# Taxonomic hierarchy only: one terminal node per primary BH-supported genus.
# Circular layout is a descriptive adaptation; no sequence phylogeny is inferred.

sgt_fan_geometry <- function(tr,taxonomy) {
  tr <- ape::reorder.phylo(tr,'postorder');nt<-length(tr$tip.label);nn<-nt+tr$Nnode
  children<-split(tr$edge[,2],tr$edge[,1]);root<-setdiff(tr$edge[,1],tr$edge[,2])[1]
  walk<-function(v)if(v<=nt)v else unlist(lapply(children[[as.character(v)]],walk),use.names=FALSE)
  ordertips<-walk(root);pos<-count<-numeric(nn);pos[ordertips]<-seq_len(nt);count[seq_len(nt)]<-1
  group<-rep(NA_character_,nn);group[seq_len(nt)]<-taxonomy$branch_colour_group[match(tr$tip.label,taxonomy$feature_id)]
  visit<-function(v){if(v<=nt)return(invisible(NULL));ch<-children[[as.character(v)]]
    for(w in ch)visit(w)
    pos[v]<<-mean(pos[ch]);count[v]<<-sum(count[ch]);z<-unique(group[ch])
    group[v]<<-if(length(z)==1L)z else 'Other phyla'
  };visit(root)
  # A topological fan layout, NOT evolutionary branch lengths. Large ancestral
  # groups are placed near the centre instead of wasting rings on single ranks.
  radius<-.78*(1-log(count)/log(nt));gap<-12*pi/180;step<-(2*pi-gap)/nt
  theta<-pi/2+gap/2+(pos-.5)*step+60*pi/180
  xy<-function(r,a)cbind(x=r*cos(a),y=r*sin(a))
  a<-xy(radius[tr$edge[,1]],theta[tr$edge[,2]]);b<-xy(radius[tr$edge[,2]],theta[tr$edge[,2]])
  edges<-data.frame(x=a[,1],y=a[,2],xend=b[,1],yend=b[,2],group=group[tr$edge[,2]])
  arcs<-ov_bind(lapply(as.integer(names(children)),function(v){
    aa<-range(theta[children[[as.character(v)]]]);n<-max(3,ceiling(diff(aa)/(.6*pi/180)))
    q<-xy(radius[v],seq(aa[1],aa[2],length.out=n));data.frame(x=q[,1],y=q[,2],node=v,group=group[v])
  }))
  tips<-data.frame(feature_id=tr$tip.label,theta=theta[seq_len(nt)],radius=radius[seq_len(nt)],
    leaf_order=pos[seq_len(nt)],step=step,group=group[seq_len(nt)])
  list(edges=edges,arcs=arcs,tips=tips,
    node_layout=data.frame(node=seq_len(nn),descendant_genera=count,layout_radius=radius))
}

sgt_build <- function(selected,baseline_N,families,limits,cfg) {
  p24 <- lr_adopted(cfg$project_root,'24_Two_generation_reanalysis')
  p36 <- lr_adopted(cfg$project_root,'36_Genus_generation_taxonomy')
  source_files <- c(file.path(p24,'private/stage_prepare.rds'),file.path(p36,'private/analysis_result.rds'))
  protected <- lr_hashes(source_files)
  ft <- readRDS(source_files[1])$features;mapping <- readRDS(source_files[2])$mapping
  s <- selected
  resolve <- function(id,rank,title) {
    ids <- mapping$feature_id[tf_bool(mapping$included)&mapping$genus==id]
    a <- gg_extract(ft$taxonomy[match(ids,ft$feature_id)],rank)
    a <- unique(unlist(a,use.names=FALSE));a <- a[!grepl('unclassified|unknown|unidentified',a,ignore.case=TRUE)]
    if(length(a)==1L)sub(paste0('^',rank,'__'),'',a)else paste('Unresolved',title)
  }
  for(rank in c('k','c','o'))s[[switch(rank,k='kingdom',c='class',o='order')]] <-
    vapply(s$feature_id,resolve,character(1),rank=rank,title=switch(rank,k='kingdom',c='class',o='order'))
  s <- s[order(s$phylum,s$class,s$order,s$family,s$feature_id),,drop=FALSE]
  rownames(s) <- NULL
  stopifnot(!anyDuplicated(s$feature_id),all(s$P_BH<.05),all(tf_bool(s$analysis_valid)),
    all(is.finite(s$overall_prevalence_percent)),all(s$overall_mean_RA_percent>0))
  nt <- nrow(s);root <- nt+1L;cursor <- root
  ranks <- c('kingdom','phylum','class','order','family')
  keys <- lapply(seq_along(ranks),function(k)apply(s[,ranks[seq_len(k)],drop=FALSE],1,paste,collapse='|'))
  nodes <- lapply(keys,function(k){u<-unique(k);z<-setNames(cursor+seq_along(u),u);cursor<<-cursor+length(u);z})
  ed <- list(cbind(root,unname(nodes[[1]])))
  node_table <- list(data.frame(node=root,rank='Root',label='Selected genera'))
  for(k in seq_along(nodes)) {
    u <- names(nodes[[k]]);at <- match(u,keys[[k]])
    if(k>1)ed[[length(ed)+1L]] <- cbind(unname(nodes[[k-1]][keys[[k-1]][at]]),unname(nodes[[k]]))
    node_table[[length(node_table)+1L]] <- data.frame(node=unname(nodes[[k]]),rank=ranks[k],label=s[[ranks[k]]][at])
  }
  ed[[length(ed)+1L]] <- cbind(unname(nodes[[length(nodes)]][keys[[length(nodes)]]]),seq_len(nt))
  edges <- do.call(rbind,ed)
  stopifnot(!anyNA(edges),nrow(edges)==cursor-1L)
  tr <- list(edge=edges,tip.label=s$feature_id,Nnode=cursor-nt,
    edge.length=rep(1,nrow(edges)))
  class(tr) <- 'phylo';storage.mode(tr$edge) <- 'integer';tr <- ape::reorder.phylo(tr)
  s$branch_colour_group <- s$phylum
  s$branch_colour_group[s$phylum%in%c('Methanobacteriota','Thermoplasmatota')] <- 'Archaea'
  groups <- c('Bacillota','Bacteroidota','Actinomycetota','Pseudomonadota',
    'Verrucomicrobiota','Archaea','Other phyla')
  s$branch_colour_group[!s$branch_colour_group%in%groups] <- 'Other phyla'
  colours <- setNames(c('#040676','#397FC7','#F1B656','#79CAFB','#F7A24F',
    '#6A95A3','#B4B4B4'),groups)
  g <- sgt_fan_geometry(tr,s)
  s <- s[match(g$tips$feature_id,s$feature_id),,drop=FALSE]
  s$leaf_order <- as.integer(g$tips$leaf_order)
  s$tip_index <- sprintf('%03d',s$leaf_order)
  s$theta <- g$tips$theta;s$step <- g$tips$step;s$radius <- g$tips$radius
  s$label_shown <- TRUE
  s$evidence_bar <- -log10(s$P_raw)
  s$prevalence_fill_fraction <- s$overall_prevalence_percent/100
  s$RA_fill_fraction <- (log10(s$overall_mean_RA_percent)-limits$RA[1])/diff(limits$RA)
  stopifnot(all(s$prevalence_fill_fraction>=0&s$prevalence_fill_fraction<=1),
    all(s$RA_fill_fraction>=-1e-10&s$RA_fill_fraction<=1+1e-10),
    all(s$beta_log2>=limits$disease[1]&s$beta_log2<=limits$disease[2]))
  s$RA_fill_fraction <- pmin(1,pmax(0,s$RA_fill_fraction))
  node_table[[length(node_table)+1L]] <- data.frame(node=seq_len(nt),rank='Genus',label=s$genus)
  node_table <- ov_bind(node_table)
  edge_table <- data.frame(parent=tr$edge[,1],child=tr$edge[,2],
    meaning='Taxonomic rank connection; length has no evolutionary meaning')
  count <- as.data.frame(table(factor(s$branch_colour_group,levels=groups)),stringsAsFactors=FALSE)
  names(count) <- c('branch_colour_group','genera')
  list(taxonomy=s,tree=tr,geometry=g,colours=colours,limits=limits,
    node_table=node_table,edge_table=edge_table,counts=count,N=baseline_N,families=families,
    model_N=unique(s$N),cases=unique(s$cases),protected=protected,
    evidence_bar_max=ceiling(max(s$evidence_bar)))
}

sgt_cells <- function(r,cfg) {
  s <- r$taxonomy;out <- list()
  association_cfg <- cfg
  if(!is.null(r$association_colors))association_cfg$binary_colors <- r$association_colors
  for(k in 1:3) {
    tip_name_column <- if(is.null(r$tip_name_column))'genus'else r$tip_name_column
    z <- s[,c('feature_id',tip_name_column,'tip_index','leaf_order','theta','step')]
    z$track <- k
    z$value <- switch(k,s$overall_prevalence_percent,s$overall_mean_RA_percent,s$abundance_ratio)
    z$track_name <- c('Prevalence (%)','Mean RA (%)','Model 3 abundance ratio')[k]
    z$fill <- switch(k,gp_tint(s$prevalence_fill_fraction,unname(cfg$binary_colors['Reference'])),
      gp_tint(s$RA_fill_fraction,unname(cfg$binary_colors['Comparison'])),
      sg_ramp(s$beta_log2,r$limits$disease,association_cfg))
    z$q <- if(k==3)s$P_BH else NA_real_
    out[[k]] <- z
  }
  ov_bind(out)
}

sgt_ring_polygons <- function(cells) {
  lows <- c(1.27,1.335,1.40);hi <- lows[cells$track]+.065
  lo <- lows[cells$track]
  a <- cbind(cells$theta-cells$step*.5,cells$theta,cells$theta+cells$step*.5,
    cells$theta+cells$step*.5,cells$theta,cells$theta-cells$step*.5)
  radius <- cbind(hi,hi,hi,lo,lo,lo)
  data.frame(x=as.vector(t(radius*cos(a))),y=as.vector(t(radius*sin(a))),
    cell=rep(seq_len(nrow(cells)),each=6),fill=rep(cells$fill,each=6))
}

sgt_plot_compact <- function(r,cfg) {
  g <- r$geometry;s <- r$taxonomy
  association_cfg <- cfg
  if(!is.null(r$association_colors))association_cfg$binary_colors <- r$association_colors
  cells <- sgt_cells(r,cfg);rings <- sgt_ring_polygons(cells)
  dots <- data.frame(x=.96*cos(s$theta),y=.96*sin(s$theta),group=s$branch_colour_group)
  tags <- s[s$label_shown,,drop=FALSE]
  tags$x <- 1.475*cos(tags$theta);tags$y <- 1.475*sin(tags$theta)
  # Sparse angle-index labels do not remove any terminal node or data cell.
  p <- ggplot2::ggplot()+
    ggplot2::geom_segment(data=g$edges,ggplot2::aes(x,y,xend=xend,yend=yend,colour=group),linewidth=.20)+
    ggplot2::geom_path(data=g$arcs,ggplot2::aes(x,y,group=node,colour=group),linewidth=.20)+
    ggplot2::geom_point(data=dots,ggplot2::aes(x,y,colour=group),size=.45,show.legend=FALSE)+
    ggplot2::scale_colour_manual(values=r$colours,guide='none')+
    ggplot2::geom_polygon(data=rings,ggplot2::aes(x,y,group=cell,fill=fill),colour=NA)+
    ggplot2::scale_fill_identity()+
    ggplot2::geom_text(data=tags,ggplot2::aes(x,y,label=tip_index),family='Arial',
      size=5.2/(72.27/25.4),show.legend=FALSE)+
    ggplot2::coord_equal(xlim=c(-1.60,1.60),ylim=c(-3.08,1.60),expand=FALSE,clip='off')+
    ggplot2::labs(title='Steatosis-associated genera: taxonomic tree',x=NULL,y=NULL)+gp_circular_theme()
  # Three observed data rings; no branch-length scale or strain/genotype proxies.
  for(k in 1:3)p <- gp_text(p,0,c(1.0625,1.2025,1.3425)[k],c('S','P','A')[k],5.5)
  p <- gp_text(p,0,-1.68,sprintf('%d significant genera; %s participants',
    nrow(s),format(r$N,big.mark=',')),6.2)
  rects <- texts <- list()
  txt <- function(x,y,label,size=5.7,hjust=.5) {
    texts[[length(texts)+1L]] <<- data.frame(x=x,y=y,label=as.character(label),
      size=size/(72.27/25.4),hjust=hjust)
  }
  left <- -.42;right <- 1.42
  for(k in 1:3) {
    y <- c(-1.89,-2.14,-2.39)[k];v <- seq(0,1,length.out=121)
    fill <- switch(k,sg_ramp(r$limits$disease[1]+v*diff(r$limits$disease),r$limits$disease,association_cfg),
      gp_tint(v,unname(cfg$binary_colors['Reference'])),gp_tint(v,unname(cfg$binary_colors['Comparison'])))
    rects[[k]] <- data.frame(xmin=left+v*(right-left),xmax=left+pmin(1,v+1/120)*(right-left),
      ymin=y,ymax=y+.052,fill=fill)
    txt(-1.49,y+.026,c('S: Steatosis ratio','P: Prevalence (%)','A: Mean RA (%)')[k],5.7,0)
    vals <- switch(k,2^c(r$limits$disease[1],0,r$limits$disease[2]),c(0,50,100),
      10^c(r$limits$RA[1],mean(r$limits$RA),r$limits$RA[2]))
    labs <- if(k==3)format(signif(vals,2),trim=TRUE,scientific=TRUE)else format(signif(vals,3),trim=TRUE,scientific=FALSE)
    for(j in 1:3)txt(left+c(0,.5,1)[j]*(right-left),y-.065,labs[j],5.3)
  }
  txt(-1.49,-2.60,'Phylum',5.8,0)
  used <- names(r$colours)[names(r$colours)%in%s$branch_colour_group]
  for(i in seq_along(used)) {
    x <- c(-1.49,-.45,.58)[(i-1)%%3+1];y <- -2.73-floor((i-1)/3)*.12
    rects[[length(rects)+1L]] <- data.frame(xmin=x,xmax=x+.06,ymin=y-.02,ymax=y+.02,fill=unname(r$colours[used[i]]))
    txt(x+.09,y,used[i],5.5,0)
  }
  p+ggplot2::geom_rect(data=ov_bind(rects),ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill),
    colour=NA,inherit.aes=FALSE,show.legend=FALSE)+
    ggplot2::geom_text(data=ov_bind(texts),ggplot2::aes(x,y,label=label,size=size,hjust=hjust),
      family='Arial',inherit.aes=FALSE,show.legend=FALSE)+ggplot2::scale_size_identity()
}

sgt_name_positions <- function(r) {
  s <- r$taxonomy
  if(is.null(r$tip_name_column)) {
    s$display_name <- ifelse(s$unnamed_genome_genus_bin,s$genus,gsub('_',' ',s$genus,fixed=TRUE))
    s$fontface <- ifelse(s$unnamed_genome_genus_bin,'plain','italic')
  } else {
    s$display_name <- s[[r$tip_name_column]]
    s$fontface <- 'plain'
  }
  mixed_names <- !is.null(r$tip_scientific_name_column)
  if(mixed_names) {
    scientific <- s[[r$tip_scientific_name_column]]
    named <- !is.na(scientific)&nzchar(scientific)
    s$label_main_text <- ifelse(named,scientific,s$display_name)
    s$label_suffix_text <- ifelse(named,paste0(' (',s[[r$tip_name_column]],')'),'')
    s$display_name <- paste0(s$label_main_text,s$label_suffix_text)
    s$fontface <- ifelse(named,'italic','plain')
  }
  s$label_size_pt <- 5.1
  # Actual Arial text metrics on the same physical named-figure canvas.
  canvas_width <- if(is.null(r$canvas_width_mm))360 else r$canvas_width_mm
  canvas_height <- if(is.null(r$canvas_height_mm))288 else r$canvas_height_mm
  tmp <- tempfile(fileext='.pdf');grDevices::cairo_pdf(tmp,width=canvas_width/25.4,height=canvas_height/25.4,family='Arial')
  on.exit({grDevices::dev.off();unlink(tmp)},add=TRUE)
  measure <- function(label,face,size)grid::convertWidth(grid::grobWidth(
    grid::textGrob(label,gp=grid::gpar(fontfamily='Arial',fontsize=size,fontface=face))),
    'mm',valueOnly=TRUE)
  if(mixed_names) {
    s$label_main_width_mm <- vapply(seq_len(nrow(s)),function(i)
      measure(s$label_main_text[i],s$fontface[i],s$label_size_pt[i]),numeric(1))
    s$label_suffix_width_mm <- vapply(seq_len(nrow(s)),function(i)
      measure(s$label_suffix_text[i],'plain',s$label_size_pt[i]),numeric(1))
    s$label_width_mm <- s$label_main_width_mm+s$label_suffix_width_mm
  } else {
    s$label_width_mm <- vapply(seq_len(nrow(s)),function(i)
      measure(s$display_name[i],s$fontface[i],s$label_size_pt[i]),numeric(1))
  }
  x_span <- if(is.null(r$plot_x_limits))4.37 else diff(r$plot_x_limits)
  unit_mm <- canvas_width/x_span
  s$name_outer_radius <- 1.245
  s$name_inner_radius <- s$name_outer_radius-s$label_width_mm/unit_mm
  stopifnot(all(s$name_inner_radius>.81))
  a <- (s$theta*180/pi)%%360;right <- cos(s$theta)>=0
  s$label_angle <- ifelse(right,ifelse(a>180,a-360,a),a-180)
  s$label_hjust <- ifelse(right,1,0)
  s$label_x <- s$name_outer_radius*cos(s$theta);s$label_y <- s$name_outer_radius*sin(s$theta)
  # Cairo/Arial 5.1-pt metrology: centre the recovered word quad on its radial
  # ray. This prevents a baseline-centering offset at the left/right flip.
  s$normal_centre_correction_pt <- .556
  correction <- s$normal_centre_correction_pt*25.4/72/unit_mm
  aa <- s$label_angle*pi/180
  s$label_x <- s$label_x-sin(aa)*correction
  s$label_y <- s$label_y+cos(aa)*correction
  guide_radius <- if(isTRUE(r$use_actual_tip_radius))s$radius+.01 else rep(.80,nrow(s))
  s$guide_x <- guide_radius*cos(s$theta);s$guide_y <- guide_radius*sin(s$theta)
  s$guide_xend <- (s$name_inner_radius-.025)*cos(s$theta)
  s$guide_yend <- (s$name_inner_radius-.025)*sin(s$theta)
  s$label_colour <- ifelse(s$beta_log2>0,'#B96C20','#262626')
  s
}

sgt_mixed_text_parts <- function(s,canvas_width,x_span=4.37) {
  # Place the italic name and roman identifier on the same radial baseline.
  # Two native R text runs retain editable text without plotmath shrinking.
  unit_mm <- canvas_width/x_span;right <- s$label_hjust==1
  main <- suffix <- s
  main$display_name <- s$label_main_text
  suffix$display_name <- s$label_suffix_text;suffix$fontface <- 'plain'
  main_shift <- ifelse(right,-s$label_suffix_width_mm/unit_mm,0)
  suffix_shift <- ifelse(right,0,-s$label_main_width_mm/unit_mm)
  main$label_x <- s$label_x+main_shift*cos(s$theta)
  main$label_y <- s$label_y+main_shift*sin(s$theta)
  suffix$label_x <- s$label_x+suffix_shift*cos(s$theta)
  suffix$label_y <- s$label_y+suffix_shift*sin(s$theta)
  main$label_piece <- 'Scientific name or unresolved SGB'
  suffix$label_piece <- 'Roman SGB identifier'
  rbind(main,suffix[nzchar(suffix$display_name),,drop=FALSE])
}

sgt_plot <- function(r,cfg) {
  g <- r$geometry;s <- sgt_name_positions(r)
  association_cfg <- cfg
  if(!is.null(r$association_colors))association_cfg$binary_colors <- r$association_colors
  x_limits <- if(is.null(r$plot_x_limits))c(-1.75,2.62)else r$plot_x_limits
  y_limits <- if(is.null(r$plot_y_limits))c(-1.75,1.75)else r$plot_y_limits
  labels <- if(is.null(r$tip_scientific_name_column))s else sgt_mixed_text_parts(s,r$canvas_width_mm,diff(x_limits))
  cells <- sgt_cells(r,cfg);rings <- sgt_ring_polygons(cells)
  # Radial bars map the existing primary-model evidence, never genotype counts.
  a <- cbind(s$theta-s$step*.34,s$theta+s$step*.34,s$theta+s$step*.34,s$theta-s$step*.34)
  hi <- 1.48+.17*s$evidence_bar/r$evidence_bar_max
  radius <- cbind(hi,hi,rep(1.48,nrow(s)),rep(1.48,nrow(s)))
  bars <- data.frame(x=as.vector(t(radius*cos(a))),y=as.vector(t(radius*sin(a))),cell=rep(seq_len(nrow(s)),each=4))
  # Exact circular band outlines in the open fan; genus text stays inside rings.
  bounds <- ov_bind(lapply(c(1.27,1.335,1.40,1.465),function(rad){
    th <- seq(min(s$theta-s$step/2),max(s$theta+s$step/2),length.out=1500)
    data.frame(x=rad*cos(th),y=rad*sin(th),ring=rad)
  }))
  p <- ggplot2::ggplot()+
    ggplot2::geom_segment(data=s,ggplot2::aes(x=guide_x,y=guide_y,xend=guide_xend,yend=guide_yend),
      colour='#B9B9B9',linetype=3,linewidth=.12)+
    ggplot2::geom_segment(data=g$edges,ggplot2::aes(x,y,xend=xend,yend=yend,colour=group),linewidth=.18)+
    ggplot2::geom_path(data=g$arcs,ggplot2::aes(x,y,group=node,colour=group),linewidth=.18)+
    ggplot2::scale_colour_manual(values=r$colours,guide='none')+
    ggplot2::geom_polygon(data=rings,ggplot2::aes(x,y,group=cell,fill=fill),colour=NA)+
    ggplot2::scale_fill_identity()+
    ggplot2::geom_path(data=bounds,ggplot2::aes(x,y,group=ring),colour='#292929',linewidth=.20)+
    ggplot2::geom_polygon(data=bars,ggplot2::aes(x,y,group=cell),fill='#6A95A3',colour='#466871',linewidth=.08)+
    ggplot2::geom_text(data=labels,ggplot2::aes(x=label_x,y=label_y,label=display_name,
      angle=label_angle,hjust=label_hjust,fontface=fontface),colour=labels$label_colour,
      family='Arial',size=5.1/(72.27/25.4),show.legend=FALSE)+
    ggplot2::coord_equal(xlim=x_limits,ylim=y_limits,expand=FALSE,clip='off')+
    ggplot2::labs(title=NULL,x=NULL,y=NULL,tag=NULL)+ggplot2::theme_void(base_family='Arial')+
    ggplot2::theme(plot.margin=ggplot2::margin(1,1,1,1,'pt'))
  # Compact right-hand keys follow the provided reference's visual composition.
  rects <- polys <- texts <- list()
  legend_multiplier <- if(is.null(r$legend_font_multiplier))1 else r$legend_font_multiplier
  stopifnot(is.finite(legend_multiplier),legend_multiplier>0)
  txt <- function(x,y,label,size=5.8,hjust=.5,colour='black') {
    texts[[length(texts)+1L]] <<- data.frame(x=x,y=y,label=as.character(label),
      size=size*legend_multiplier/(72.27/25.4),hjust=hjust,colour=colour)
  }
  x0 <- if(is.null(r$legend_key_left))1.91 else r$legend_key_left
  x1 <- if(is.null(r$legend_key_right))2.49 else r$legend_key_right
  for(k in 1:3) {
    y <- c(.69,.26,-.17)[k];v <- seq(0,1,length.out=121)
    txt((x0+x1)/2,y+.18,c('Prevalence (%)','Mean relative\nabundance (%)','Steatosis\nabundance ratio')[k],6)
    if(k<=2) {
      colour <- unname(cfg$binary_colors[c('Reference','Comparison')[k]])
      for(j in seq_len(length(v)-1L)) {
        t <- c(v[j],v[j+1],v[j+1],v[j]);xx <- x0+t*(x1-x0)
        yy <- y+c(0,0,.075*v[j+1],.075*v[j])
        polys[[length(polys)+1L]] <- data.frame(x=xx,y=yy,cell=paste(k,j),fill=gp_tint(mean(v[j:(j+1)]),colour))
      }
      outline <- data.frame(x=c(x0,x1,x1,x0),y=y+c(0,0,.075,0))
      p <- p+ggplot2::geom_path(data=outline,ggplot2::aes(x,y),inherit.aes=FALSE,colour='#454545',linewidth=.18)
      vals <- if(k==1)c(0,100)else 10^r$limits$RA
      labs <- if(k==1)as.character(vals)else format(signif(vals,2),scientific=TRUE,trim=TRUE)
      txt(x0,y-.055,labs[1],5.3);txt(x1,y-.055,labs[2],5.3)
    }else{
      rects[[length(rects)+1L]] <- data.frame(xmin=x0+v*(x1-x0),xmax=x0+pmin(1,v+1/120)*(x1-x0),
        ymin=y,ymax=y+.065,fill=sg_ramp(r$limits$disease[1]+v*diff(r$limits$disease),r$limits$disease,association_cfg))
      vals <- 2^c(r$limits$disease[1],0,r$limits$disease[2])
      for(j in 1:3)txt(x0+c(0,.5,1)[j]*(x1-x0),y-.055,format(signif(vals[j],3),trim=TRUE),5.3)
    }
  }
  txt((x0+x1)/2,-.45,'-log10(raw P)',5.9)
  for(j in 1:2) {
    v <- c(5,r$evidence_bar_max)[j];y <- -.54-(j-1)*.10
    rects[[length(rects)+1L]] <- data.frame(xmin=x0,xmax=x0+.5*v/r$evidence_bar_max,
      ymin=y-.015,ymax=y+.015,fill='#6A95A3')
    txt(x0-.04,y,as.character(v),5.3,1)
  }
  txt((x0+x1)/2,-.82,'Higher with steatosis',5.5,colour='#B96C20')
  txt(x0,-1.00,'Phylum',5.8,0)
  used <- names(r$colours)[names(r$colours)%in%s$branch_colour_group]
  swatch_half_height <- if(is.null(r$legend_swatch_half_height)).012 else r$legend_swatch_half_height
  for(i in seq_along(used)) {
    y <- -1.10-(i-1)*.075
    rects[[length(rects)+1L]] <- data.frame(xmin=x0,xmax=x0+.045,ymin=y-swatch_half_height,ymax=y+swatch_half_height,fill=unname(r$colours[used[i]]))
    txt(x0+.075,y,used[i],5.3,0)
  }
  p+ggplot2::geom_polygon(data=ov_bind(polys),ggplot2::aes(x,y,group=cell,fill=fill),
    colour=NA,inherit.aes=FALSE,show.legend=FALSE)+
    ggplot2::geom_rect(data=ov_bind(rects),ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill),
      colour=NA,inherit.aes=FALSE,show.legend=FALSE)+
    ggplot2::geom_text(data=ov_bind(texts),ggplot2::aes(x,y,label=label,size=size,hjust=hjust),
      colour=ov_bind(texts)$colour,family='Arial',inherit.aes=FALSE,show.legend=FALSE)+ggplot2::scale_size_identity()
}

sgt_tables <- function(r,out,cfg) {
  s <- r$taxonomy[order(r$taxonomy$leaf_order),,drop=FALSE]
  v3_table(s,'Figure39_Significant_genus_circular_tree_tip_key_EN',out)
  v3_table(r$node_table,'Figure39_Significant_genus_taxonomy_nodes_EN',out)
  v3_table(r$edge_table,'Figure39_Significant_genus_taxonomy_edges_EN',out)
  v3_table(r$geometry$node_layout,'Figure39_Topological_fan_layout_EN',out)
  v3_table(r$counts,'Figure39_Significant_genus_tree_phylum_counts_EN',out)
  v3_table(sgt_cells(r,cfg),
    'Figure39_Significant_genus_circular_tree_cells_EN',out)
  v3_table(sgt_name_positions(r),'Figure39_Full_name_radial_label_positions_EN',out)
  v3_table(s[,c('feature_id','genus','P_raw','P_BH','evidence_bar')],
    'Figure39_Outer_bar_original_P_evidence_EN',out)
  key <- data.frame(Index=s$tip_index,Genus=s$genus,Phylum=s$phylum,
    `Adjusted abundance ratio`=format(signif(s$abundance_ratio,3),trim=TRUE),
    `Prevalence (%)`=format(round(s$overall_prevalence_percent,1),trim=TRUE),
    `Mean RA (%)`=format(signif(s$overall_mean_RA_percent,3),trim=TRUE),
    `BH q`=format(signif(s$P_BH,3),scientific=TRUE,trim=TRUE),check.names=FALSE)
  esc <- function(x)gsub('<','&lt;',gsub('&','&amp;',as.character(x),fixed=TRUE),fixed=TRUE)
  lines <- c('<!doctype html><meta charset="utf-8"><title>Significant genus tree: complete index</title>',
    '<style>body{font:14px Arial;margin:24px}table{border-collapse:collapse}th,td{padding:6px 12px;border-bottom:1px solid #ddd;text-align:left}th{position:sticky;top:0;background:white;border-bottom:2px solid black}</style>',
    '<h2>Significant genera: complete circular-tree index</h2><p>All 265 primary BH discoveries; one terminal node per genus/GGB. The figure prints every full genus/GGB label. Table indices follow counterclockwise order from the upper-left fan gap.</p>',
    paste0('<table><thead><tr>',paste0('<th>',esc(names(key)),'</th>',collapse=''),'</tr></thead><tbody>'))
  for(i in seq_len(nrow(key)))lines <- c(lines,paste0('<tr>',paste0('<td>',esc(unlist(key[i,])),'</td>',collapse=''),'</tr>'))
  writeLines(c(lines,'</tbody></table>'),file.path(out,'tables/Figure39_Significant_genus_circular_tree_tip_key.html'),useBytes=TRUE)
  tr <- r$tree;tr$edge.length <- NULL
  ape::write.tree(tr,file=file.path(out,'tables/Figure39_Significant_genus_taxonomy_no_branch_lengths.nwk'))
}

sgt_legend <- function(r,out) {
  writeLines(c('# Figure 39 | Taxonomic distribution of steatosis-associated genera.',
    sprintf('All %d genus/genome-genus groups meeting the original full-family BH q<0.05 in the primary all-age Model 3 are shown, with exactly one terminal node and one full radial name per group. These are %d named genera and %d unnamed GGB groups. No other taxa or constituent SGB tips are added. Connections encode the recorded kingdom-phylum-class-order-family-genus hierarchy; class/order/kingdom annotations are conservatively resolved from the original constituent-SGB taxonomy, with ambiguous ranks explicitly unresolved. Radial rank spacing has no evolutionary-distance meaning. This is a taxonomic cladogram, not an inferred sequence phylogeny.',
      nrow(r$taxonomy),sum(!r$taxonomy$unnamed_genome_genus_bin),sum(r$taxonomy$unnamed_genome_genus_bin)),
    sprintf('Branches are coloured by phylum, with Archaea and other phyla shown as broad display groups. From inside out, the three rings show prevalence (blue, linear percentage), mean relative abundance (orange, common log10-positive-mean colour intensity), and the adjusted Model 3 abundance ratio for steatosis versus no steatosis (blue-white-orange centred at 1). The first two summaries use all %s baseline microbiome participants, include non-detection zeros, and retain the original full-profile denominator without renormalizing selected genera. Disease estimates use the unchanged model sample of %s participants (%s cases).',
      format(r$N,big.mark=','),format(r$model_N,big.mark=','),format(r$cases,big.mark=',')),
    'Taxa are ordered by their recorded taxonomic hierarchy and identifier, independently of abundance and association strength. Every terminal node, full genus/GGB name and all three data cells are plotted. The topological fan positions internal nodes by descendant-group count for compactness; its radial coordinates have no evolutionary or temporal meaning. Grey dotted lines lead from terminal branches to the full radial names. Orange-brown names identify genera with higher abundance in steatosis; dark names identify lower abundance. Outermost teal bars encode -log10(original raw P), not recovered-genotype counts. All identifiers, label geometry, original values and evidence bars are provided in the accompanying source tables. Shared colour limits match the complete heatmaps. The fully named reference-style canvas is 360 x 288 mm to retain a >=5-pt text floor; this is a large-format reading master, not an asserted 183-mm journal-production figure. No selected-subset black rings are added to the volcano.',
    'Model 3 is the saved official MaAsLin2 log2 relative-abundance mixed model, adjusted for generation, generation-specific age splines, sex, preparation batch, sequencing pool, sampling season and recent antibiotics, with a primary-family random intercept. Ratios refer to zero-replaced relative abundance, not absolute microbial load. No statistical fits or multiple-testing corrections are repeated.',
    'The circular tree-plus-ring composition is visually adapted from Valles-Colomer M, et al. Nature Microbiology. 2022;7:87-96. Figure 2. doi:10.1038/s41564-021-01021-8. The use of a taxonomic hierarchy to organize quantitative community data follows the descriptive framework discussed by Foster ZSL, Sharpton TJ, Gruenwald NJ. Metacoder: An R package for visualization and manipulation of community taxonomic diversity data. PLOS Computational Biology. 2017;13:e1005404. doi:10.1371/journal.pcbi.1005404. The current figure uses custom R geometry, not the metacoder algorithm or a replication of the source paper\'s maximum-likelihood inference, StrainPhlAn genotyping or pTR analysis.'),
    file.path(out,'reports/Figure39_Significant_genus_circular_taxonomy_legend_EN.md'),useBytes=TRUE)
}
