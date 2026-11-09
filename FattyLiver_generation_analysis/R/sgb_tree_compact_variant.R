# R-only display variant. The saved tree, all 379 SGBs and inference stay fixed.
scv_load <- function(root) {
  sys.source(file.path(root,'R/significant_sgb_circular_tree.R'),envir=.GlobalEnv)
  cfg<-sst_setup(root)
  parent<-lr_adopted(root,'40_Steatosis_SGB_circular_tree')
  r<-readRDS(file.path(parent,'private/display_result.rds'))
  stopifnot(identical(r$protocol,sst_protocol),nrow(r$taxonomy)==379L)
  list(cfg=cfg,parent=parent,result=r)
}
scv_probe <- function(x) {
  r<-x$result;s<-sgt_name_positions(r)
  factor<-min((s$name_inner_radius-.05)/s$radius)
  print(data.frame(maximum_collision_free_uniform_tree_scale=factor))
  print(stats::quantile(s$radius))
  print(stats::quantile(s$name_inner_radius-s$radius))
  run<-cumsum(c(TRUE,head(s$branch_colour_group,-1)!=tail(s$branch_colour_group,-1)))
  print(ov_bind(lapply(split(seq_len(nrow(s)),run),function(i)data.frame(
    group=s$branch_colour_group[i[1]],n=length(i),first=min(i),last=max(i),
    theta_lo=min(s$theta[i])*180/pi,theta_hi=max(s$theta[i])*180/pi,
    min_radius=min(s$radius[i]),max_radius=max(s$radius[i])))))
  invisible(x)
}
scv_prepare <- function(original) {
  r<-original;s<-sgt_name_positions(r)
  # One uniform scale preserves every relative reference branch length.
  multiplier<-min(1.30,.98*min((s$name_inner_radius-.05)/s$radius))
  for(k in c('edges','arcs'))for(v in intersect(c('x','y','xend','yend'),names(r$geometry[[k]])))
    r$geometry[[k]][[v]]<-r$geometry[[k]][[v]]*multiplier
  r$geometry$tips$radius<-r$geometry$tips$radius*multiplier
  r$geometry$scale<-r$geometry$scale*multiplier
  r$taxonomy$radius<-r$taxonomy$radius*multiplier
  r$compact_tree_scale<-multiplier;r
}
scv_tree_segments <- function(r) {
  g<-r$geometry;a<-g$arcs
  same<-head(a$node,-1)==tail(a$node,-1)
  arcs<-data.frame(x=head(a$x,-1)[same],y=head(a$y,-1)[same],
    xend=tail(a$x,-1)[same],yend=tail(a$y,-1)[same],group=head(a$group,-1)[same])
  ov_bind(list(g$edges,arcs))
}
scv_box_interval <- function(seg,box) {
  dx<-seg$xend-seg$x;dy<-seg$yend-seg$y
  tx1<-(box$xmin-seg$x)/dx;tx2<-(box$xmax-seg$x)/dx
  ty1<-(box$ymin-seg$y)/dy;ty2<-(box$ymax-seg$y)/dy
  zx<-abs(dx)<1e-12;zy<-abs(dy)<1e-12
  tx1[zx]<-ifelse(seg$x[zx]>=box$xmin&seg$x[zx]<=box$xmax,-Inf,Inf);tx2[zx]<-Inf
  ty1[zy]<-ifelse(seg$y[zy]>=box$ymin&seg$y[zy]<=box$ymax,-Inf,Inf);ty2[zy]<-Inf
  lo<-pmax(0,pmin(tx1,tx2),pmin(ty1,ty2));hi<-pmin(1,pmax(tx1,tx2),pmax(ty1,ty2))
  data.frame(lo=lo,hi=hi,hit=is.finite(lo)&is.finite(hi)&lo<=hi)
}
scv_label_positions <- function(r) {
  s<-sgt_name_positions(r);s<-s[order(s$leaf_order),,drop=FALSE]
  tree<-scv_tree_segments(r);unit_mm<-r$canvas_width_mm/diff(r$plot_x_limits)
  # Conservative boxes for EVERY radial SGB/species label, rather than a
  # nearest-angle corner test alone (long names can reach adjacent sectors).
  ux<-cos(s$theta);uy<-sin(s$theta);pad<-.014
  label_x<-cbind((s$name_inner_radius-.02)*ux-pad*uy,
    (s$name_inner_radius-.02)*ux+pad*uy,(s$name_outer_radius+.015)*ux-pad*uy,
    (s$name_outer_radius+.015)*ux+pad*uy)
  label_y<-cbind((s$name_inner_radius-.02)*uy+pad*ux,
    (s$name_inner_radius-.02)*uy-pad*ux,(s$name_outer_radius+.015)*uy+pad*ux,
    (s$name_outer_radius+.015)*uy-pad*ux)
  species_boxes<-data.frame(xmin=apply(label_x,1,min),xmax=apply(label_x,1,max),
    ymin=apply(label_y,1,min),ymax=apply(label_y,1,max))
  tmp<-tempfile(fileext='.pdf')
  grDevices::cairo_pdf(tmp,width=r$canvas_width_mm/25.4,height=r$canvas_height_mm/25.4,family='Arial')
  on.exit({grDevices::dev.off();unlink(tmp)},add=TRUE)
  groups<-unique(s$branch_colour_group);runs<-cumsum(c(TRUE,head(s$branch_colour_group,-1)!=tail(s$branch_colour_group,-1)))
  anchors<-ov_bind(lapply(groups,function(g){
    blocks<-split(which(s$branch_colour_group==g),runs[s$branch_colour_group==g])
    ix<-blocks[[which.max(lengths(blocks))]]
    data.frame(label=g,n=length(ix),theta=mean(range(s$theta[ix])),
      tip_radius=stats::median(s$radius[ix]),theta_lo=min(s$theta[ix]),theta_hi=max(s$theta[ix]))
  }))
  # Long names first. Search only vacant reference-tree space, keeping labels
  # near their true largest contiguous colour block; leaders identify anchors.
  anchors<-anchors[order(-nchar(anchors$label),anchors$label),,drop=FALSE]
  answer<-list();placed<-data.frame()
  for(j in seq_len(nrow(anchors))){
    a<-anchors[j,];size<-12
    grob<-grid::textGrob(a$label,gp=grid::gpar(fontfamily='Arial',fontsize=size,fontface='bold'))
    halfw<-(grid::convertWidth(grid::grobWidth(grob),'mm',valueOnly=TRUE)/2+.9)/unit_mm
    halfh<-(grid::convertHeight(grid::grobHeight(grob),'mm',valueOnly=TRUE)/2+1.0)/unit_mm
    span<-if(a$n>100)55 else 32
    cand<-expand.grid(theta=a$theta+seq(-span,span,length.out=65)*pi/180,
      radius=seq(.36,1.10,by=.02))
    cand$x<-cand$radius*cos(cand$theta);cand$y<-cand$radius*sin(cand$theta)
    cand$score<-((cand$theta-a$theta)/(.35))^2+((cand$radius-.88)/.4)^2
    cand<-cand[order(cand$score),,drop=FALSE];found<-NULL
    for(i in seq_len(nrow(cand))){
      c<-cand[i,];box<-data.frame(xmin=c$x-halfw,xmax=c$x+halfw,ymin=c$y-halfh,ymax=c$y+halfh)
      if(nrow(placed)&&any(box$xmin<placed$xmax&box$xmax>placed$xmin&
        box$ymin<placed$ymax&box$ymax>placed$ymin))next
      if(any(box$xmin<species_boxes$xmax&box$xmax>species_boxes$xmin&
        box$ymin<species_boxes$ymax&box$ymax>species_boxes$ymin))next
      corners<-expand.grid(x=c(box$xmin,box$xmax),y=c(box$ymin,box$ymax))
      angles<-atan2(corners$y,corners$x)
      near<-vapply(angles,function(t)which.min(abs(atan2(sin(s$theta-t),cos(s$theta-t)))),integer(1))
      if(any(sqrt(corners$x^2+corners$y^2)>=s$name_inner_radius[near]-.03))next
      if(any(scv_box_interval(tree,box)$hit))next
      found<-cbind(a,c[,c('x','y')],box,size_pt=size)
      break
    }
    if(is.null(found))stop('No clear direct-label position for ',a$label)
    placed<-rbind(placed,found[,c('xmin','xmax','ymin','ymax')]);answer[[j]]<-found
  }
  z<-ov_bind(answer)
  text_colours<-c(Bacillota='#040676',Bacteroidota='#397FC7',Actinomycetota='#AC7925',
    Pseudomonadota='#397F9E',Verrucomicrobiota='#C57D36',Archaea='#4E7885','Other phyla'='#707070')
  z$colour<-unname(text_colours[z$label])
  z$anchor_radius<-pmin(1.01,z$tip_radius+.10)
  z$anchor_x<-z$anchor_radius*cos(z$theta);z$anchor_y<-z$anchor_radius*sin(z$theta)
  z
}
scv_clip_guides <- function(lines,boxes) {
  # Clear only non-phylogenetic terminal guides under group names; actual
  # reference edges/arcs are never deleted, clipped or rerouted.
  for(j in seq_len(nrow(boxes))){
    z<-scv_box_interval(lines,boxes[j,]);hit<-which(z$hit)
    if(!length(hit))next
    pieces<-list(lines[!z$hit,,drop=FALSE])
    for(i in hit){
      a<-lines[i,,drop=FALSE];dx<-a$xend-a$x;dy<-a$yend-a$y
      if(z$lo[i]>1e-8){b<-a;b$xend<-a$x+dx*z$lo[i];b$yend<-a$y+dy*z$lo[i];pieces[[length(pieces)+1L]]<-b}
      if(z$hi[i]<1-1e-8){b<-a;b$x<-a$x+dx*z$hi[i];b$y<-a$y+dy*z$hi[i];pieces[[length(pieces)+1L]]<-b}
    }
    lines<-ov_bind(pieces)
  }
  lines
}
scv_align_group_colours <- function(r,labels) {
  # Shared ancestral branches are neutral connectors, not members of the
  # taxonomic catch-all "Other phyla". Terminal categories stay unchanged.
  tr<-ape::reorder.phylo(r$tree,'postorder');nt<-length(tr$tip.label)
  group<-rep(NA_character_,nt+tr$Nnode)
  group[seq_len(nt)]<-r$taxonomy$branch_colour_group[match(tr$tip.label,r$taxonomy$feature_id)]
  children<-split(tr$edge[,2],tr$edge[,1]);root<-setdiff(tr$edge[,1],tr$edge[,2])[1]
  walk<-function(node){
    if(node<=nt)return(group[node])
    values<-vapply(children[[as.character(node)]],walk,character(1))
    group[node]<<-if(length(unique(values))==1L)values[1]else 'Shared ancestor'
    group[node]
  }
  walk(root)
  stopifnot(nrow(r$geometry$edges)==nrow(tr$edge))
  r$geometry$edges$group<-group[tr$edge[,2]]
  r$geometry$arcs$group<-group[r$geometry$arcs$node]
  # Exact same named palette for branches, terminal guides, callouts and text.
  r$colours['Shared ancestor']<-'#D0D0D0'
  labels$colour<-unname(r$colours[labels$label])
  s<-r$taxonomy[order(r$taxonomy$leaf_order),,drop=FALSE]
  radial<-sgt_name_positions(r);radial<-radial[match(s$feature_id,radial$feature_id),]
  finish<-pmax(s$radius+.005,pmin(1.06,radial$name_inner_radius-.045))
  fraction<-if(is.null(r$guide_extension_fraction))1 else r$guide_extension_fraction
  finish<-s$radius+.004+fraction*(finish-s$radius-.004)
  runs<-cumsum(c(TRUE,head(s$branch_colour_group,-1)!=tail(s$branch_colour_group,-1)))
  for(j in seq_len(nrow(labels))){
    ids<-which(s$branch_colour_group==labels$label[j])
    blocks<-split(ids,runs[ids]);ids<-blocks[[which.max(lengths(blocks))]]
    cand<-expand.grid(index=ids,t=seq(.1,.95,length.out=18))
    cand$radius<-s$radius[cand$index]+.004+cand$t*(finish[cand$index]-s$radius[cand$index]-.004)
    cand$x<-cand$radius*cos(s$theta[cand$index]);cand$y<-cand$radius*sin(s$theta[cand$index])
    clear<-vapply(seq_len(nrow(cand)),function(i)!any(cand$x[i]>labels$xmin-.01&cand$x[i]<labels$xmax+.01&
      cand$y[i]>labels$ymin-.01&cand$y[i]<labels$ymax+.01),logical(1))
    cand<-cand[clear,,drop=FALSE];stopifnot(nrow(cand)>0L)
    score<-(cand$x-labels$x[j])^2+(cand$y-labels$y[j])^2+
      .12*atan2(sin(s$theta[cand$index]-labels$theta[j]),cos(s$theta[cand$index]-labels$theta[j]))^2
    best<-cand[which.min(score),]
    labels$anchor_x[j]<-best$x;labels$anchor_y[j]<-best$y;labels$anchor_radius[j]<-best$radius
    labels$anchor_feature_id[j]<-s$feature_id[best$index]
    labels$anchor_branch_group[j]<-s$branch_colour_group[best$index]
  }
  stopifnot(identical(labels$label,labels$anchor_branch_group),
    identical(labels$colour,unname(r$colours[labels$label])))
  r$align_group_colours<-TRUE
  list(result=r,labels=labels)
}
scv_plot <- function(r,cfg,labels) {
  s<-sgt_name_positions(r)
  extended<-pmax(s$radius+.005,pmin(1.06,s$name_inner_radius-.045))
  # Optional display-only follow-up: shorten only the added coloured guide;
  # keep the actual tree, frozen group names and outer tracks identical.
  if(!is.null(r$guide_extension_fraction)){
    stopifnot(r$guide_extension_fraction>0,r$guide_extension_fraction<=1)
    start<-s$radius+.004
    extended<-start+r$guide_extension_fraction*(extended-start)
  }
  coloured<-data.frame(feature_id=s$feature_id,x=(s$radius+.004)*cos(s$theta),
    y=(s$radius+.004)*sin(s$theta),xend=extended*cos(s$theta),yend=extended*sin(s$theta),
    group=s$branch_colour_group)
  grey<-data.frame(feature_id=s$feature_id,x=(extended+.005)*cos(s$theta),
    y=(extended+.005)*sin(s$theta),xend=(s$name_inner_radius-.025)*cos(s$theta),
    yend=(s$name_inner_radius-.025)*sin(s$theta),group=s$branch_colour_group)
  coloured<-scv_clip_guides(coloured,labels);grey<-scv_clip_guides(grey,labels)
  p<-sgt_plot(r,cfg)
  # All source polygons, radial species labels, bars and numeric keys are reused.
  p$layers[[1]]<-ggplot2::geom_segment(data=coloured,
    ggplot2::aes(x,y,xend=xend,yend=yend,colour=group),linewidth=.22,
    alpha=if(isTRUE(r$align_group_colours))1 else .65)
  p$layers[[2]]$aes_params$linewidth<-.34;p$layers[[3]]$aes_params$linewidth<-.34
  p<-p+ggplot2::geom_segment(data=grey,ggplot2::aes(x,y,xend=xend,yend=yend),
    inherit.aes=FALSE,colour='#B9B9B9',linetype=3,linewidth=.12)
  # The direct labels replace the categorical right-side phylum key only.
  n<-length(p$layers)
  for(i in seq_len(n)){
    d<-p$layers[[i]]$data
    if(is.data.frame(d)&&all(c('label','hjust','size','colour')%in%names(d))){
      keep<-!d$label%in%c('Phylum',names(r$colours))
      p$layers[[i]]$data<-d[keep,,drop=FALSE]
      if(length(p$layers[[i]]$aes_params$colour)==nrow(d))
        p$layers[[i]]$aes_params$colour<-p$layers[[i]]$aes_params$colour[keep]
    }
    if(is.data.frame(d)&&all(c('xmin','xmax','ymin','ymax','fill')%in%names(d)))
      p$layers[[i]]$data<-d[d$ymin> -1,,drop=FALSE]
  }
  leaders<-data.frame(x=labels$anchor_x,y=labels$anchor_y,xend=labels$x,yend=labels$y,
    group=labels$label)
  leaders<-scv_clip_guides(leaders,labels)
  p<-p+ggplot2::geom_segment(data=leaders,ggplot2::aes(x,y,xend=xend,yend=yend,colour=group),
    inherit.aes=FALSE,linewidth=.24)+
    ggplot2::geom_text(data=labels,ggplot2::aes(x,y,label=label),inherit.aes=FALSE,
      colour=labels$colour,family='Arial',fontface='bold',size=12/(72.27/25.4))
  if(isTRUE(r$align_group_colours))p<-p+ggplot2::geom_point(data=labels,
    ggplot2::aes(anchor_x,anchor_y,colour=label),inherit.aes=FALSE,shape=21,
    fill='white',size=.85,stroke=.26)
  p
}
scv_main <- function(args=character(),root) {
  if(getRversion()<'4.5'){
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',
      shQuote(file.path(root,'scripts/40b_SGB_tree_compact_variant.R')),vapply(args,shQuote,character(1))))
    if(status!=0L)stop('Project-R compact SGB-tree variant failed')
    return(invisible(NULL))
  }
  stopifnot(length(args)<=1L,all(args%in%c('--inspect','--run','--preview','--publish')))
  x<-scv_load(root);if('--inspect'%in%args)return(scv_probe(x))
  if('--publish'%in%args){
    variants<-list.dirs(file.path(dirname(dirname(x$parent)),'variants'),recursive=FALSE,full.names=TRUE)
    variants<-variants[order(file.info(variants)$mtime,decreasing=TRUE)]
    stem<-'Figure40_SGB_circular_tree_compact_phylum_labels_20261103'
    complete<-vapply(variants,function(d)all(file.exists(file.path(d,'figures',paste0(stem,c('.pdf','.png','.tiff','.svg'))))),logical(1))
    stopifnot(any(complete));out<-variants[which(complete)[1]]
    report<-jsonlite::fromJSON(file.path(out,'audit',paste0(stem,'.collision-audit.json')))
    stopifnot(identical(report$verdict,'PASS'),report$summary$fail==0,report$summary$warn==0)
    hashes<-tf_read(file.path(out,'audit/original_file_fingerprints.csv'))
    return(scv_publish_variant(x,out,stem,hashes,root))
  }
  r<-scv_prepare(x$result);cfg<-x$cfg
  id<-paste0('compact_phylum_labels_',format(Sys.time(),'%Y%m%d_%H%M%S'))
  out<-file.path(dirname(dirname(x$parent)),'variants',id)
  for(d in c('private','audit','figures','tables','reports'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
  stem<-'Figure40_SGB_circular_tree_compact_phylum_labels_20261103'
  public<-file.path(root,'Results_表格与图片/01_当前分析/16_脂肪肝相关SGB_环形树')
  originals<-c(file.path(x$parent,'figures',paste0(sst_stem,c('.pdf','.png','.tiff','.svg'))),
    file.path(public,'figures',paste0(sst_stem,c('.pdf','.png','.tiff','.svg'))),
    file.path(x$parent,'private/display_result.rds'))
  stopifnot(all(file.exists(originals)))
  hashes<-lr_hashes(originals);fl_write_csv(hashes,file.path(out,'audit/original_file_fingerprints.csv'))
  labels<-scv_label_positions(r);p<-scv_plot(r,cfg,labels)
  fl_write_csv(labels,file.path(out,'tables',paste0(stem,'_direct_group_labels_EN.csv')))
  saveRDS(list(source_run=x$parent,result=r,labels=labels),file.path(out,'private/display_variant.rds'))
  grDevices::png(file.path(out,'audit/compact_variant_preview.png'),width=r$canvas_width_mm,
    height=r$canvas_height_mm,units='mm',res=300,type='cairo',bg='white')
  tryCatch(print(p),finally=grDevices::dev.off())
  cat('VARIANT:',out,'\n');print(labels[,c('label','x','y')])
  if('--preview'%in%args)return(invisible(out))
  manifest<-tf_export(list(plot=patchwork::wrap_plots(list(p)),width=r$canvas_width_mm,
    height=r$canvas_height_mm,panels=1L,exemptions=list(),fixed_aspect_single_panel=TRUE,
    no_panel_tags=TRUE,oriented_text=TRUE,max_body_pt=16.51,
    large_format_font_reason='Separate user-requested layout variant of the existing 590x448-mm reading master; preserve original 5.1-pt SGB labels and enlarged numeric keys; new direct phylum names are 12 pt'),
    stem,out,cfg)
  source0<-sgt_cells(x$result,cfg);source1<-sgt_cells(r,cfg)
  stopifnot(identical(source0,source1),identical(r$tree,x$result$tree),
    identical(r$all_results,x$result$all_results),identical(r$taxonomy$feature_id,x$result$taxonomy$feature_id),
    identical(r$taxonomy$figure_display_label,x$result$taxonomy$figure_display_label),
    setequal(labels$label,unique(r$taxonomy$branch_colour_group)))
  lr_check_hashes(hashes);lr_check_hashes(r$protected);lr_check_hashes(r$inherited_protected)
  fl_write_csv(manifest,file.path(out,'tables',paste0(stem,'_export_manifest_EN.csv')))
  fl_write_csv(source1,file.path(out,'tables',paste0(stem,'_unchanged_ring_source_EN.csv')))
  writeLines(c('# Compact reference-tree layout variant',
    'All 379 original significant SGBs, exact saved tree topology and relative branch lengths, original annotations, labels, ring values, colours and P/q are retained.',
    sprintf('The actual reference tree has a single uniform display scale multiplier of %.6f; branch strokes increase from 0.18 to 0.34 mm.',r$compact_tree_scale),
    'Lighter coloured terminal guide extensions connect the phylogenetic tips toward the unchanged name band. These display-only extensions are not additional phylogenetic branch lengths. A short grey dotted tail distinguishes the guide-to-name connection.',
    'Direct phylum/broad-group labels use the largest contiguous tip block of each recorded colour group, with short leaders. Archaea and Other phyla are broad display groups, as in the original. Only non-phylogenetic guides are interrupted beneath these names; actual reference-tree segments are never removed or rerouted.',
    'This separate layout is saved beside the original figure. Original images, source RDS and inference are byte-preserved. No statistical models, BH adjustment or tree inference are rerun.'),
    file.path(out,'reports',paste0(stem,'_layout_note_EN.md')))
  scv_publish_variant(x,out,stem,hashes,root)
}
scv_publish_variant <- function(x,out,stem,hashes,root,replace_existing=FALSE) {
  # Add files only. No run adoption, historical move or overwrite of originals.
  public<-file.path(root,'Results_表格与图片/01_当前分析/16_脂肪肝相关SGB_环形树')
  lr_check_hashes(hashes)
  for(d in c('figures','tables','reports')){
    files<-list.files(file.path(out,d),full.names=TRUE)
    dest<-file.path(x$parent,d,basename(files))
    existing<-file.exists(dest)
    if(any(existing)&&!replace_existing)stopifnot(identical(unname(tools::md5sum(files[existing])),unname(tools::md5sum(dest[existing]))))
    if(any(existing)&&replace_existing){
      # User explicitly authorizes replacing only this latest short-line stem.
      stopifnot(all(startsWith(basename(dest[existing]),stem)))
      backup<-file.path(out,'private/prior');dir.create(backup,showWarnings=FALSE)
      before<-lr_hashes(dest[existing])
      before$backup_file<-paste0(d,'_',seq_len(nrow(before)),'.',tools::file_ext(dest[existing]))
      stopifnot(all(file.copy(dest[existing],file.path(backup,before$backup_file),overwrite=FALSE)))
      fl_write_csv(before,file.path(out,'audit',paste0('replaced_',d,'_fingerprints.csv')))
      stopifnot(all(file.copy(files[existing],dest[existing],overwrite=TRUE)))
    }
    if(any(!existing))stopifnot(all(file.copy(files[!existing],file.path(x$parent,d),overwrite=FALSE)))
  }
  # Keep the full strict audit in the variant directory. Nesting another long
  # filename below the old adopted-run name exceeds Windows MAX_PATH in R.
  writeLines(c('Full strict export audit and preserved original fingerprints:',
    file.path(out,'audit')),file.path(x$parent,'audit/Figure40b_compact_variant_audit_location.txt'))
  lr_check_hashes(hashes)
  fl_refresh_result_gallery(root)
  lr_check_hashes(hashes)
  expected<-file.path(public,'figures',paste0(stem,c('.pdf','.png','.tiff','.svg')))
  stopifnot(all(file.exists(expected)),identical(unname(tools::md5sum(expected)),
    unname(tools::md5sum(file.path(out,'figures',basename(expected))))))
  writeLines(c('PASS: strict font, oriented collision and alignment checks.',
    'PASS: 379 SGBs, original tree and inference, three data rings and names unchanged.',
    'PASS: all protected original/preceding images and the saved RDS byte-preserved.',
    'PASS: new four-format files synchronized to the original viewing folder.'),
    file.path(out,'audit/delivery_validation.txt'))
  cat('COMPLETED:',file.path(public,'figures',stem),'\n');invisible(out)
}

scv_short_main <- function(args=character(),root) {
  if(getRversion()<'4.5'){
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',
      shQuote(file.path(root,'scripts/40c_SGB_tree_shorter_lines.R')),vapply(args,shQuote,character(1))))
    if(status!=0L)stop('Project-R shorter SGB-tree guides failed')
    return(invisible(NULL))
  }
  stopifnot(length(args)<=1L,all(args%in%c('--run','--publish')))
  x<-scv_load(root);cfg<-x$cfg
  module<-dirname(dirname(x$parent));variants<-list.dirs(file.path(module,'variants'),recursive=FALSE,full.names=TRUE)
  stem<-'Figure40_SGB_circular_tree_short_lines_20261103'
  if('--publish'%in%args){
    keep<-vapply(variants,function(p)all(file.exists(file.path(p,'figures',paste0(stem,c('.pdf','.png','.tiff','.svg'))))),logical(1))
    stopifnot(any(keep));available<-variants[keep];out<-available[which.max(file.info(available)$mtime)]
    qa<-jsonlite::fromJSON(file.path(out,'audit',paste0(stem,'.collision-audit.json')))
    stopifnot(qa$verdict=='PASS',qa$summary$fail==0,qa$summary$warn==0)
    return(scv_publish_variant(x,out,stem,tf_read(file.path(out,'audit/original_file_fingerprints.csv')),root,replace_existing=TRUE))
  }
  previous_stem<-'Figure40_SGB_circular_tree_compact_phylum_labels_20261103'
  complete<-vapply(variants,function(p)all(file.exists(file.path(p,'figures',paste0(previous_stem,c('.pdf','.png','.tiff','.svg')))))&&
    file.exists(file.path(p,'private/display_variant.rds')),logical(1))
  stopifnot(any(complete));previous<-variants[complete][which.max(file.info(variants[complete])$mtime)]
  old<-readRDS(file.path(previous,'private/display_variant.rds'))
  r<-old$result;labels<-old$labels;r$guide_extension_fraction<-.50
  aligned<-scv_align_group_colours(r,labels);r<-aligned$result;labels<-aligned$labels
  out<-file.path(module,'variants',paste0('short_lines_',format(Sys.time(),'%Y%m%d_%H%M%S')))
  for(d in c('private','audit','figures','tables','reports'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
  public<-file.path(root,'Results_表格与图片/01_当前分析/16_脂肪肝相关SGB_环形树')
  names<-as.vector(outer(c(sst_stem,previous_stem),c('.pdf','.png','.tiff','.svg'),paste0))
  protected<-c(file.path(x$parent,'figures',names),
    file.path(public,'figures',paste0(sst_stem,c('.pdf','.png','.tiff','.svg'))),
    file.path(x$parent,'private/display_result.rds'))
  stopifnot(all(file.exists(protected)))
  optional_previous_view<-file.path(public,'figures',paste0(previous_stem,c('.pdf','.png','.tiff','.svg')))
  protected<-c(protected,optional_previous_view[file.exists(optional_previous_view)])
  hashes<-lr_hashes(protected)
  fl_write_csv(hashes,file.path(out,'audit/original_file_fingerprints.csv'))
  p<-scv_plot(r,cfg,labels)
  cat('SHORTER-LINE VARIANT:',out,'\n')
  manifest<-tf_export(list(plot=patchwork::wrap_plots(list(p)),width=r$canvas_width_mm,
    height=r$canvas_height_mm,panels=1L,exemptions=list(),fixed_aspect_single_panel=TRUE,
    no_panel_tags=TRUE,oriented_text=TRUE,max_body_pt=16.51,
    large_format_font_reason='Same previously approved 590x448-mm reading-master canvas, 5.1-pt SGB names, 12-pt group names and enlarged keys; only the coloured terminal guide length is halved'),
    stem,out,cfg)
  # Explicitly freeze every previously approved component except guide length.
  a<-sgt_name_positions(r);start<-a$radius+.004
  before<-pmax(a$radius+.005,pmin(1.06,a$name_inner_radius-.045))
  after<-start+.5*(before-start)
  fl_write_csv(data.frame(feature_id=a$feature_id,start_radius=start,
    previous_guide_end=before,new_guide_end=after,coloured_extension_fraction=.5),
    file.path(out,'tables',paste0(stem,'_guide_lengths_EN.csv')))
  geom_same<-all(vapply(c('edges','arcs'),function(n)identical(
    r$geometry[[n]][,setdiff(names(r$geometry[[n]]),'group'),drop=FALSE],
    old$result$geometry[[n]][,setdiff(names(old$result$geometry[[n]]),'group'),drop=FALSE]),logical(1)))
  frozen_label_columns<-c('label','x','y','xmin','xmax','ymin','ymax','size_pt')
  stopifnot(geom_same,identical(r$taxonomy,old$result$taxonomy),
    identical(r$tree,old$result$tree),identical(sgt_cells(r,cfg),sgt_cells(old$result,cfg)),
    identical(r$all_results,old$result$all_results),identical(labels[,frozen_label_columns],old$labels[,frozen_label_columns]),
    identical(labels$colour,unname(r$colours[labels$label])),
    identical(labels$label,labels$anchor_branch_group),
    all(abs((after-start)/pmax(before-start,1e-12)-.5)<1e-10))
  lr_check_hashes(hashes);lr_check_hashes(r$protected);lr_check_hashes(r$inherited_protected)
  saveRDS(list(source_run=x$parent,previous_variant=previous,result=r,labels=labels),file.path(out,'private/display_variant.rds'))
  fl_write_csv(labels,file.path(out,'tables',paste0(stem,'_unchanged_group_labels_EN.csv')))
  fl_write_csv(data.frame(group=names(r$colours),hex_colour=unname(r$colours),
    mapping='One exact palette for reference branches, terminal extensions, group names, callout leaders and anchor marks; shared ancestors neutral'),
    file.path(out,'tables',paste0(stem,'_colour_mapping_EN.csv')))
  fl_write_csv(manifest,file.path(out,'tables',paste0(stem,'_export_manifest_EN.csv')))
  writeLines(c('# Shorter terminal-guide layout',
    'Relative to the preceding compact layout, each added coloured terminal guide is shortened to 50% of its former length; its remaining connection is a faint grey dotted guide.',
    'The actual reference-tree coordinates, line widths, all 379 SGBs, group-name positions, outer rings/bars, numerical values, canvas and keys are unchanged.',
    'Names, terminal branches, extensions, leader lines and anchor marks use the exact same per-group HEX colour, without a darker approximate text shade or alpha fading. Each callout ends on an actual annotated terminal branch of that group, with its SGB ID recorded in the source table. Heterogeneous shared ancestors have a separate neutral grey and are not labelled as Other phyla.',
    'The user explicitly requested overwriting the preceding short-line four-format results. Earlier original and long-line compact figures remain byte-preserved, and the replaced short-line payload is retained privately for audit. No statistical models or BH corrections are rerun.'),
    file.path(out,'reports',paste0(stem,'_layout_note_EN.md')))
  scv_publish_variant(x,out,stem,hashes,root,replace_existing=TRUE)
}
