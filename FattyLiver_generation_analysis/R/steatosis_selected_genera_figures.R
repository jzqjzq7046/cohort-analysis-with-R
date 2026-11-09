# R-only structural adaptation of the existing taxonomy map and ratio forest plot.
# Every selected genus is shown. Pagination changes layout, never selection/FDR.
sg_ramp <- function(v, limits, cfg) {
  p <- grDevices::colorRamp(c(cfg$binary_colors['Reference'],'#FFFFFF',cfg$binary_colors['Comparison']))
  ans <- rep('#E9E9E9',length(v));ok <- is.finite(v)
  stopifnot(all(v[ok]>=limits[1]-1e-10),all(v[ok]<=limits[2]+1e-10))
  ans[ok] <- grDevices::rgb(p((v[ok]-limits[1])/diff(limits))/255)
  ans
}

sg_scales <- function(r) {
  p <- r$profiles[r$profiles$feature_id%in%r$selected$feature_id,]
  v <- p$mean_RA_percent[p$mean_RA_percent>0]
  ra <- if(length(v))c(floor(min(log10(v))*2)/2,ceiling(max(log10(v))*2)/2)else c(-2,0)
  if(diff(ra)==0)ra <- ra+c(-.5,.5)
  a <- r$pairs[r$pairs$model=='P2' & r$pairs$analysis_valid,]
  limit <- function(v,step)if(length(v))max(step,ceiling(max(abs(v),na.rm=TRUE)/step)*step)else step
  e <- limit(a$beta_log2_per_log2,.05);d <- limit(r$selected$beta_log2,.5)
  list(RA=ra,disease=c(-d,d),pair=c(-e,e))
}

sg_cells <- function(r,s,cfg,limits) {
  ids <- s$feature_id;tracks <- list()
  tracks[[1]] <- data.frame(feature_id=ids,track=1L,value=s$beta_log2,q=s$P_BH,
    y=(nrow(s)+1-s$display_row)*4)
  for(i in 1:6) {
    gen <- c('G0','G1','G2')[(i-1)%%3+1]
    a <- r$profiles[r$profiles$generation==gen,];at <- match(ids,a$feature_id)
    tracks[[i+1]] <- data.frame(feature_id=ids,track=i+1L,
      value=if(i<=3)a$prevalence_percent[at]else a$mean_RA_percent[at],q=NA_real_,
      y=(nrow(s)+1-s$display_row)*4)
  }
  for(i in 1:3) {
    a <- r$pairs[r$pairs$model=='P2' & r$pairs$comparison==gp_pairs[i],]
    at <- match(ids,a$feature_id);v <- a$beta_log2_per_log2[at]
    valid <- !is.na(at) & tf_bool(a$analysis_valid[at]);v[!valid] <- NA_real_
    tracks[[i+7]] <- data.frame(feature_id=ids,track=i+7L,value=v,q=a$P_BH[at],
      y=(nrow(s)+1-s$display_row)*4)
  }
  z <- ov_bind(tracks);z$x <- c(118,135,142,149,160,167,174,187,194,201)[z$track]
  z$fill <- '#E9E9E9';z$star <- '';z$ratio <- NA_real_
  for(i in 1:10) {
    j <- which(z$track==i)
    if(i%in%c(2:7)) {
      v <- if(i<=4)z$value[j]/100 else gp_ra_fraction(z$value[j],limits$RA)
      z$fill[j] <- gp_tint(v,cfg$generation_colors[c('G0','G1','G2')[(i-2)%%3+1]])
    } else {
      lim <- if(i==1)limits$disease else limits$pair
      z$fill[j] <- sg_ramp(z$value[j],lim,cfg)
      z$ratio[j] <- 2^z$value[j]
      z$star[j[is.finite(z$q[j]) & z$q[j]<.05]] <- '*'
    }
  }
  z$meaning <- c('Model 3 steatosis abundance ratio','G0 detection (%)','G1 detection (%)',
    'G2 detection (%)','G0 mean RA (%)','G1 mean RA (%)','G2 mean RA (%)',
    'G0-G1 Model 6 ratio','G1-G2 Model 6 ratio','G0-G2 Model 6 ratio')[z$track]
  z
}

sg_key <- function(p,cfg,limits) {
  # Vectorize the legend into two layers rather than hundreds of annotations.
  texts <- rectangles <- list()
  txt <- function(x,y,label,size=5.7,hjust=.5) {
    texts[[length(texts)+1L]] <<- data.frame(x=x,y=y,label=as.character(label),
      size=size/(72.27/25.4),hjust=hjust)
  }
  rect <- function(x0,x1,y0,y1,col) {
    rectangles[[length(rectangles)+1L]] <<- data.frame(xmin=x0,xmax=x1,ymin=y0,ymax=y1,fill=col)
  }
  centers <- c(68,119,170);width <- 35
  for(j in 1:3) {
    gen <- c('G0','G1','G2')[j];txt(centers[j],-2,gen)
    for(v in seq(0,1,length.out=81)) {
      for(y in c(-8,-20))rect(centers[j]-width/2+v*width,
        centers[j]-width/2+min(1,v+1/80)*width,y,y+2,gp_tint(v,cfg$generation_colors[gen]))
    }
    for(v in c(0,.5,1))txt(centers[j]-width/2+v*width,-11,as.character(100*v),5.3)
    ticks <- c(10^limits$RA[1],1,10^limits$RA[2]);ticks <- sort(unique(ticks[ticks>0 &
      log10(ticks)>=limits$RA[1] & log10(ticks)<=limits$RA[2]]))
    for(v in ticks)txt(centers[j]-width/2+(log10(v)-limits$RA[1])/diff(limits$RA)*width,
      -23,format(v,trim=TRUE,scientific=FALSE,digits=2),5.3)
  }
  txt(2,-7,'Detection (%)',5.7,0);txt(2,-19,'Mean RA (%)',5.7,0)
  for(k in 1:2) {
    lim <- if(k==1)limits$disease else limits$pair;y <- if(k==1)-32 else -44
    for(v in seq(0,1,length.out=101))rect(68+v*120,68+min(1,v+1/100)*120,
      y,y+2,sg_ramp(lim[1]+v*diff(lim),lim,cfg))
    for(v in c(lim[1],0,lim[2]))txt(68+(v-lim[1])/diff(lim)*120,y-3.5,
      format(signif(2^v,3),trim=TRUE,scientific=FALSE),5.3)
    txt(2,y+1,if(k==1)'Steatosis ratio'else'Cross-generation ratio',5.7,0)
  }
  txt(68,-54,'* BH q < 0.05',5.3,0)
  rect(164,167,-55,-53,'#E9E9E9');txt(169,-54,'Not tested / NE',5.3,0)
  p+ggplot2::geom_rect(data=ov_bind(rectangles),ggplot2::aes(xmin=xmin,xmax=xmax,
    ymin=ymin,ymax=ymax,fill=fill),colour=NA,inherit.aes=FALSE,show.legend=FALSE)+
    ggplot2::geom_text(data=ov_bind(texts),ggplot2::aes(x,y,label=label,size=size,hjust=hjust),
      family='Arial',inherit.aes=FALSE,show.legend=FALSE)+ggplot2::scale_size_identity()
}

sg_map <- function(r,s,cfg,limits,page,out,source_prefix='Figure38') {
  s$display_row <- seq_len(nrow(s));g <- gg_geometry(list(selected=s))
  # Stop connector strokes at the widened phylum boxes, including hidden paths.
  g$edges$xend[g$edges$x==1 & g$edges$xend==5] <- 2
  g$edges$x[g$edges$x==27 & g$edges$xend==30] <- 28
  # Abbreviate only this verbose provisional family label on the drawing.
  g$families$family_label <- ifelse(grepl('_Family_[IVX]+_Incertae_Sedis$',g$families$family),
    sub('^([A-Za-z]+)_Family_([IVX]+)_Incertae_Sedis$','\\1 \\2',g$families$family),
    g$families$family)
  v3_table(unique(g$families[,c('family','family_label')]),
    sprintf('%s_family_display_labels_page_%02d_EN',source_prefix,page),out)
  d <- sg_cells(r,s,cfg,limits);d$page <- page
  v3_table(d,sprintf('%s_map_page_%02d_source_EN',source_prefix,page),out)
  top <- nrow(s)*4
  p <- ggplot2::ggplot()+
    ggplot2::geom_segment(data=g$edges,ggplot2::aes(x,y,xend=xend,yend=yend),colour='#808080',linewidth=.25)+
    ggplot2::geom_rect(data=g$phyla,ggplot2::aes(xmin=2,xmax=28,ymin=y-1.3,ymax=y+1.3),
      fill='white',colour='#B5B5B5',linewidth=.2)+
    ggplot2::geom_rect(data=g$families,ggplot2::aes(xmin=33,xmax=62,ymin=y-1.3,ymax=y+1.3),
      fill='white',colour='#B5B5B5',linewidth=.2)+
    ggplot2::geom_text(data=g$phyla,ggplot2::aes(x=15,y=y,label=phylum),size=5.2/(72.27/25.4),family='Arial')+
    ggplot2::geom_text(data=g$families,ggplot2::aes(x=47.5,y=y,label=family_label),size=5.2/(72.27/25.4),family='Arial')+
    ggplot2::geom_text(data=g$selected,ggplot2::aes(x=70,y=y,label=genus,
      fontface=ifelse(unnamed_genome_genus_bin,'plain','italic')),hjust=0,size=6/(72.27/25.4),family='Arial')+
    ggplot2::geom_tile(data=d,ggplot2::aes(x,y,fill=fill),width=6.3,height=3.5)+
    ggplot2::geom_text(data=d[d$star!='',],ggplot2::aes(x,y,label=star),size=6/(72.27/25.4),family='Arial')+
    ggplot2::scale_fill_identity()+
    ggplot2::coord_cartesian(xlim=c(-1,206),ylim=c(-58,top+24),expand=FALSE,clip='off')+
    ggplot2::theme_void(base_family='Arial')+
    ggplot2::theme(plot.margin=ggplot2::margin(7,5,5,5,'pt'))
  for(i in 1:3)p <- gp_text(p,c(15,47.5,88)[i],top+13,c('Phylum','Family','Genus')[i],6)
  for(i in 1:4)p <- gp_text(p,c(118,142,167,194)[i],top+19,
    c('Steatosis','Detection (%)','Mean RA (%)','Cross-generation')[i],5.7)
  for(i in 1:10)p <- gp_text(p,c(118,135,142,149,160,167,174,187,194,201)[i],top+9,
    if(i==1)'Model 3'else if(i<=7)c('G0','G1','G2')[(i-2)%%3+1]else gp_pairs[i-7],5.5,angle=90,
    colour=if(i%in%2:7)cfg$generation_colors[c('G0','G1','G2')[(i-2)%%3+1]]else'black')
  sg_key(p,cfg,limits)
}

sg_forest <- function(s,cfg,title) {
  stopifnot(all(is.finite(s$abundance_ratio)),all(s$abundance_ratio>0),
    all(is.finite(s$ratio_CI_low)),all(s$ratio_CI_low>0),
    all(is.finite(s$ratio_CI_high)),all(s$ratio_CI_high>=s$ratio_CI_low))
  s <- s[order(s$association_evidence_rank),,drop=FALSE]
  s$genus <- factor(s$genus,levels=rev(s$genus))
  s$direction_label <- factor(ifelse(s$beta_log2<0,'Lower with steatosis','Higher with steatosis'),
    levels=c('Lower with steatosis','Higher with steatosis'))
  ggplot2::ggplot(s,ggplot2::aes(abundance_ratio,genus))+
    ggplot2::geom_vline(xintercept=1,linetype=2,colour='#888888',linewidth=.25)+
    ggplot2::geom_segment(ggplot2::aes(x=ratio_CI_low,xend=ratio_CI_high,yend=genus),colour='black',linewidth=.4)+
    ggplot2::geom_point(ggplot2::aes(fill=direction_label),shape=21,colour='black',stroke=.25,size=2)+
    ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),levels(s$direction_label)),drop=FALSE)+
    ggplot2::scale_x_log10(breaks=scales::breaks_log(n=6),
      labels=function(v)format(v,trim=TRUE,scientific=FALSE,digits=3))+
    ggplot2::labs(x='Adjusted abundance ratio (95% CI)',y=NULL,title=title,fill=NULL)+tf_plot_theme()
}

sg_preflight <- function(cfg,out) {
  path <- file.path(out,'audit/actual_plot_and_export_source.R')
  fns <- c('sg_ramp','sg_scales','sg_cells','sg_key','sg_map','sg_forest','sg_figures',
    'gg_geometry','gp_tint','gp_ra_fraction','gp_text','gp_rect','lr_export','tf_export','nf_theme')
  writeLines(unlist(lapply(fns,function(n)c(paste0(n,' <- '),deparse(get(n),width.cutoff=120),''))),path)
  nf_qa_command(cfg,'validate_figure.py',c(shQuote(path),'--json'),file.path(out,'audit/plot_source_preflight.json'))
}

sg_figures <- function(r,cfg,out) {
  s <- r$selected;manifest <- list();limits <- sg_scales(r)
  v3_table(data.frame(scale=c('Mean RA log10 percent lower','Mean RA log10 percent upper',
    'Steatosis log2 ratio lower','Steatosis log2 ratio upper','Pair log2 ratio lower','Pair log2 ratio upper'),
    value=c(limits$RA,limits$disease,limits$pair)),'Figure38_shared_colour_limits_EN',out)
  if(!nrow(s)) {
    p <- lr_empty('Steatosis-associated genera','No genus met BH q < 0.05 in the primary model')
    manifest[[1]] <- lr_export(list(p),'Figure38_No_BH_supported_genera',out,cfg,width=183,height=70)
  }else{
    # 40 is a per-page physical layout limit, not a discovery/display-count cap.
    groups <- split(seq_len(nrow(s)),ceiling(seq_len(nrow(s))/40L))
    seen <- character()
    for(k in seq_along(groups)) {
      z <- s[groups[[k]],,drop=FALSE];seen <- c(seen,z$feature_id)
      p <- sg_map(r,z,cfg,limits,k,out)
      stem <- if(length(groups)==1)'Figure38_Steatosis_selected_genus_map'else
        sprintf('Figure38_Steatosis_selected_genus_map_page_%02d',k)
      height <- max(120,(4*nrow(z)+82)*.94)
      manifest[[length(manifest)+1L]] <- lr_export(list(p),stem,out,cfg,width=183,height=height)
    }
    stopifnot(identical(seen,s$feature_id),!anyDuplicated(seen))
    evidence <- s[order(s$association_evidence_rank),,drop=FALSE]
    groups <- split(seq_len(nrow(evidence)),ceiling(seq_len(nrow(evidence))/40L))
    for(k in seq_along(groups)) {
      z <- evidence[groups[[k]],,drop=FALSE]
      v3_table(z,sprintf('Figure38_forest_page_%02d_source_EN',k),out)
      p <- sg_forest(z,cfg,'Steatosis-associated genera: Model 3')
      stem <- if(length(groups)==1)'Figure38_Steatosis_genus_abundance_ratios'else
        sprintf('Figure38_Steatosis_genus_abundance_ratios_page_%02d',k)
      manifest[[length(manifest)+1L]] <- lr_export(list(p),stem,out,cfg,width=120,
        height=max(85,35+4.1*nrow(z)))
    }
  }
  v3_table(ov_bind(manifest),'Figure_export_manifest_EN',out)
  writeLines(c('# Figure 38 | Steatosis-associated genera and their cross-generation profiles.',
    sprintf('All %d genera/genome genus groups meeting BH q<0.05 in the all-age Model 3 are shown; no fixed count or abundance ranking is used. The source discovery model includes %d participants and %d ultrasound-defined steatosis cases.',
      nrow(s),r$samples$N[r$samples$model=='B1'],r$samples$cases[r$samples$model=='B1']),
    'The taxonomy map shows phylum-family-genus classification relationships, not an inferred phylogeny. Long provisional family labels are abbreviated on the drawing; full names are retained in source tables. The steatosis column encodes the Model 3 abundance ratio for steatosis versus no steatosis. Generation tracks show detection (%) and mean relative abundance (%) in all baseline participants, including zeros. Cross-generation tracks show Model 6 younger-member abundance ratios per doubling of older-member abundance. The disease and cross-generation ratios have separate symmetric log2 colour scales, both centred at ratio 1. Asterisks indicate BH q<0.05; grey indicates not tested or not estimable. All pages use the same colour limits.',
    'The companion forest plot displays all selected genera, ordered by discovery q, raw P and identifier. Points and bars denote adjusted abundance ratios and 95% confidence intervals. Model 3 adjusts for generation, generation-specific age splines, sex, preparation batch, sequencing pool, sampling season and recent antibiotics, with a primary-family random intercept. Confidence intervals use Satterthwaite degrees of freedom. Disease ratios refer to zero-replaced relative abundance, not absolute microbial load.',
    'Model 6 adjusts both pair endpoints for age splines, sex, batch, pool, season and recent antibiotics using primary-family-clustered Gaussian GEE. Its original all-genus/all-comparison BH values are retained unchanged. Same-cohort discovery and follow-up are exploratory, not independent validation or evidence of strain transmission. All selected taxa remain in the tables even if a pair comparison was unavailable.'),
    file.path(out,'reports/Figure_legend_EN.md'),useBytes=TRUE)
}
