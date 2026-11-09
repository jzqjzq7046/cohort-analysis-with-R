# Presentation identifiers only. Never recode analysis data, formulas or saved fits.
fl_model_name_map<-c(C1='Model 1',C2='Model 2',B1='Model 3',B2='Model 4',P1='Model 5',P2='Model 6')
fl_model_labels<-function(x){
 if(!is.character(x))return(x)
 for(k in names(fl_model_name_map))x<-gsub(paste0('(?<![A-Za-z0-9])',k,'(?![A-Za-z0-9])'),fl_model_name_map[[k]],x,perl=TRUE)
 x
}
fl_model_label_grob<-function(g){
 if(!is.null(g$label)&&is.character(g$label))g$label<-fl_model_labels(g$label)
 if(length(g$grobs))g$grobs<-lapply(g$grobs,fl_model_label_grob)
 if(length(g$children))for(i in seq_along(g$children))g$children[[i]]<-fl_model_label_grob(g$children[[i]])
 g
}

# Transform labels before ggplot builds its guide/axis layout. The post-build
# fallback remains for custom grobs, but never replaces statistical identifiers.
fl_model_plot_labels<-function(p){
 if(!inherits(p,'ggplot'))return(p)
 if(length(p$patches$plots))p$patches$plots<-lapply(p$patches$plots,fl_model_plot_labels)
 if(length(p$labels))for(k in names(p$labels))if(is.character(p$labels[[k]]))p$labels[[k]]<-fl_model_labels(p$labels[[k]])
 wrap_mapping<-function(mapping){
  if(!is.null(mapping$label)){
   q<-mapping$label;mapping$label<-rlang::new_quosure(rlang::call2('fl_model_labels',rlang::get_expr(q)),rlang::get_env(q))
  };mapping
 }
 p$mapping<-wrap_mapping(p$mapping)
 if(length(p$layers))for(i in seq_along(p$layers))p$layers[[i]]$mapping<-wrap_mapping(p$layers[[i]]$mapping)
 if(length(p$scales$scales))for(i in seq_along(p$scales$scales)){
  s<-p$scales$scales[[i]]$clone();if(!s$is_discrete())next;previous<-s$labels
  if(inherits(previous,'waiver'))s$labels<-fl_model_labels else if(is.character(previous))s$labels<-fl_model_labels(previous)else if(is.function(previous)){
   s$labels<-local({f<-previous;function(x)fl_model_labels(f(x))})
  };p$scales$scales[[i]]<-s
 }
 # Missing discrete x/y scales otherwise format axis text after the existing
 # plot has chosen too little label space for the longer model names.
 b<-ggplot2::ggplot_build(p)
 # Implicit shape/linetype scales must receive the same display labels as an
 # explicit colour scale, otherwise ggplot splits their previously shared guide.
 for(ss in b$plot$scales$scales){
  axis_name<-ss$aesthetics[1]
  if(!axis_name%in%c('x','y')&&is.null(p$scales$get_scales(axis_name))&&ss$is_discrete()){
   sc<-ss$clone();sc$labels<-fl_model_labels;p$scales$add(sc)
  }
 }
 for(axis in c('x','y'))if(is.null(p$scales$get_scales(axis))){
  ss<-if(axis=='x')b$layout$panel_scales_x else b$layout$panel_scales_y
  if(length(ss)&&ss[[1]]$is_discrete()){
   if(axis=='x')p<-p+ggplot2::scale_x_discrete(labels=fl_model_labels)else p<-p+ggplot2::scale_y_discrete(labels=fl_model_labels)
  }
 }
 p
}
