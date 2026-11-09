# Inert definitions. User-approved genus-level replacement for the dense overview.
# Classificatory hierarchy, NOT a newly estimated genus phylogeny.
gg_protocol<-'genus_aggregation_all_tests_P1_P2_taxonomy_map_20261028_v1'
gg_setup<-function(root){
  sys.source(file.path(root,'R/generation_reference_phylogeny.R'),envir=.GlobalEnv)
  cfg<-gp_setup(root)
  for(f in c('R/reanalysis_numerical_helpers.R','R/taxa_followup.R'))sys.source(file.path(root,f),envir=.GlobalEnv)
  stopifnot(requireNamespace('geepack',quietly=TRUE))
  cfg
}
gg_extract<-function(s,rank){
  z<-regmatches(s,gregexpr(paste0(rank,'__[^|,; ]+'),s))
  # Exact nomenclatural synonym only: do not merge GTDB suffix groups or genera.
  if(rank=='p')z<-lapply(z,function(v)sub('^p__Firmicutes$','p__Bacillota',v))
  lapply(z,unique)
}
gg_unique<-function(z,prefix){
  z<-unique(unlist(z,use.names=FALSE));z<-z[!grepl('unclassified|unknown|unidentified',z,ignore.case=TRUE)]
  if(length(z)==1L)sub(paste0('^',prefix,'__'),'',z)else paste('Unresolved',switch(prefix,p='phylum',f='family'))
}
gg_aggregate<-function(x,out){
  ft<-x$features[match(rownames(x$relative),x$features$feature_id),];stopifnot(!anyNA(ft$feature_id))
  gs<-gg_extract(ft$taxonomy,'g');fs<-gg_extract(ft$taxonomy,'f');ps<-gg_extract(ft$taxonomy,'p')
  ok<-vapply(gs,function(v)length(v)==1L&&!grepl('unclassified|unknown|unidentified',v,ignore.case=TRUE),logical(1))
  genus<-vapply(gs,function(v)if(length(v)==1L)v else NA_character_,character(1))
  mapping<-data.frame(feature_id=ft$feature_id,genus=genus,included=ok,
    genus_annotations=vapply(gs,paste,collapse=';',character(1)),
    reason=ifelse(ok,'Unique genus annotation',ifelse(lengths(gs)>1,'Ambiguous genus','Unclassified / absent genus')))
  G<-rowsum(x$relative[ok,,drop=FALSE],genus[ok],reorder=TRUE)
  fraction<-colSums(G)
  stopifnot(all(is.finite(G)),all(G>=0),all(rowSums(G)>0),all(fraction>0),all(fraction<=1+1e-10),
    max(abs(fraction-colSums(x$relative[ok,,drop=FALSE])))<1e-10)
  groups<-split(which(ok),genus[ok])
  taxonomy<-do.call(rbind,lapply(rownames(G),function(k){ii<-groups[[k]]
    data.frame(feature_id=k,genus=sub('^g__','',k),phylum=gg_unique(ps[ii],'p'),family=gg_unique(fs[ii],'f'),
      SGBs_aggregated=length(ii),unnamed_genome_genus_bin=grepl('^g__GGB',k),
      overall_mean_RA_percent=100*mean(G[k,]),overall_prevalence_percent=100*mean(G[k,]>0))
  }))
  profile<-do.call(rbind,lapply(c('G0','G1','G2'),function(g){a<-G[,x$metadata$generation==g,drop=FALSE]
    data.frame(feature_id=rownames(G),generation=g,N=ncol(a),detected_n=rowSums(a>0),
      prevalence_percent=100*rowMeans(a>0),mean_RA_percent=100*rowMeans(a))
  }))
  coverage<-do.call(rbind,lapply(c('Overall','G0','G1','G2'),function(g){
    j<-if(g=='Overall')rep(TRUE,ncol(G))else x$metadata$generation==g
    data.frame(generation=g,N=sum(j),mean_classified_percent=100*mean(fraction[j]),
      median_classified_percent=100*median(fraction[j]),Q1_classified_percent=100*unname(quantile(fraction[j],.25)),
      Q3_classified_percent=100*unname(quantile(fraction[j],.75)),min_classified_percent=100*min(fraction[j]),
      mean_unassigned_or_ambiguous_percent=100*mean(1-fraction[j]))
  }))
  # Freeze the display set before fitting; association P values are not consulted.
  selected<-head(taxonomy[order(-taxonomy$overall_mean_RA_percent,taxonomy$feature_id),],30L)
  selected$abundance_rank<-seq_len(nrow(selected))
  selected$selection_rule<-'Top 30 genera/genome genus bins by overall mean full-profile RA; ties by ID; no P-value selection'
  # Taxonomy-based row order: descending summed display abundance at each node.
  selected$family_key<-paste(selected$phylum,selected$family,sep='|')
  pw<-tapply(selected$overall_mean_RA_percent,selected$phylum,sum)
  fw<-tapply(selected$overall_mean_RA_percent,selected$family_key,sum)
  selected<-selected[order(-pw[selected$phylum],selected$phylum,-fw[selected$family_key],selected$family,
    -selected$overall_mean_RA_percent,selected$feature_id),]
  selected$display_row<-seq_len(nrow(selected))
  taxonomy$main_display<-taxonomy$feature_id%in%selected$feature_id
  count<-data.frame(SGBs_baseline=nrow(x$relative),SGBs_assigned=sum(ok),SGBs_ambiguous_or_unclassified=sum(!ok),
    genus_groups=nrow(G),named_genus_groups=sum(!taxonomy$unnamed_genome_genus_bin),
    genome_genus_bins=sum(taxonomy$unnamed_genome_genus_bin),main_display_genera=nrow(selected),
    mean_total_RA_percent_in_display=sum(selected$overall_mean_RA_percent),renormalized=FALSE)
  saveRDS(list(relative=G,metadata=x$metadata,unassigned_fraction=1-fraction),file.path(out,'private/genus_relative_abundance.rds'))
  list(relative=G,metadata=x$metadata,mapping=mapping,taxonomy=taxonomy,profiles=profile,coverage=coverage,selected=selected,counts=count)
}

gg_harmonize_taxonomy<-function(r,out){
  # Saved-fit display repair only. Abundances, selected IDs and inference freeze.
  inference<-r$results;profiles<-r$profiles;display_ids<-sort(r$selected$feature_id)
  source<-r$protected$path[basename(r$protected$path)=='stage_prepare.rds'];stopifnot(length(source)==1L)
  ft<-readRDS(source)$features;k<-match(r$mapping$feature_id,ft$feature_id);stopifnot(!anyNA(k))
  ps<-gg_extract(ft$taxonomy[k],'p');map<-r$mapping
  groups<-split(which(map$included),map$genus[map$included])
  r$taxonomy$phylum<-vapply(r$taxonomy$feature_id,function(id)gg_unique(ps[groups[[id]]],'p'),character(1))
  r$selected$phylum<-r$taxonomy$phylum[match(r$selected$feature_id,r$taxonomy$feature_id)]
  s<-r$selected;s$family_key<-paste(s$phylum,s$family,sep='|')
  pw<-tapply(s$overall_mean_RA_percent,s$phylum,sum);fw<-tapply(s$overall_mean_RA_percent,s$family_key,sum)
  s<-s[order(-pw[s$phylum],s$phylum,-fw[s$family_key],s$family,-s$overall_mean_RA_percent,s$feature_id),]
  s$display_row<-seq_len(nrow(s));r$selected<-s
  stopifnot(identical(inference,r$results),identical(profiles,r$profiles),identical(display_ids,sort(r$selected$feature_id)))
  v3_table(data.frame(original_phylum='Firmicutes',display_phylum='Bacillota',
    relationship='Homotypic nomenclatural synonyms; exact strings only',
    source='https://lpsn.dsmz.de/phylum/bacillota'),'Taxonomy_name_harmonization_EN',out)
  r
}

gg_archive_table_aliases<-function(out){
  aliases<-c(Genus_count_EN='Genus_counts_EN',Genus_profile_EN='Genus_profiles_EN')
  for(old in names(aliases)){
    a<-file.path(out,'tables',paste0(old,'.csv'));b<-file.path(out,'tables',paste0(aliases[[old]],'.csv'))
    if(file.exists(a)&&file.exists(b)){
      stopifnot(identical(tf_read(a),tf_read(b)))
      dest<-file.path(out,'audit','legacy_table_name_aliases');dir.create(dest,showWarnings=FALSE)
      for(ext in c('.csv','.html')){f<-file.path(out,'tables',paste0(old,ext))
        if(file.exists(f)&&!file.rename(f,file.path(dest,basename(f))))stop('Cannot archive duplicated table alias')
      }
    }
  }
}

gg_fit<-function(a,p26,cfg,out){
  ids<-rownames(a$relative);A<-a$relative
  replacement<-apply(A,1,function(v)min(v[v>0])/2);stopifnot(all(is.finite(replacement)&replacement>0))
  Z<-log2(A);zero<-which(A==0,arr.ind=TRUE);Z[zero]<-log2(replacement[zero[,1]])
  fl_write_csv(data.frame(feature_id=ids,zero_replacement=replacement),file.path(out,'audit/fixed_genus_zero_replacements.csv'))
  old_counts<-tf_read(file.path(p26,'tables/Family_pair_sample_sizes_EN.csv'))
  results<-samples<-list()
  for(g in v3_pairs){
    key<-paste(g,collapse='-');m<-v3_scope(a$metadata,g)
    m<-v3_scope(m[v3_cc(m,v3_P2),,drop=FALSE],g);m<-m[order(m$family_id,m$sample_id),,drop=FALSE]
    ia<-which(m$generation==g[1]);ib<-which(m$generation==g[2])
    edge<-which(outer(as.character(m$family_id[ia]),as.character(m$family_id[ib]),'=='),arr.ind=TRUE)
    edge<-edge[order(m$family_id[ia[edge[,1]]],m$sample_id[ib[edge[,2]]],m$sample_id[ia[edge[,1]]]),,drop=FALSE]
    older<-ia[edge[,1]];younger<-ib[edge[,2]];fam<-as.character(m$family_id[older])
    stopifnot(all(fam==as.character(m$family_id[younger])))
    pairs<-data.frame(older_id=m$sample_id[older],younger_id=m$sample_id[younger],family_id=fam)
    old_pairs<-readRDS(file.path(p26,'private',paste0('Observed_pairs_',key,'.rds')))
    stopifnot(identical(pairs,old_pairs))
    saveRDS(pairs,file.path(out,'private',paste0('Observed_pairs_',key,'.rds')))
    j<-match(m$sample_id,colnames(A));stopifnot(!anyNA(j))
    pvA<-rowMeans(A[,j[ia],drop=FALSE]>0);pvB<-rowMeans(A[,j[ib],drop=FALSE]>0)
    eligible<-pvA>=.10 & pvB>=.10;kept<-ids[eligible];stopifnot(length(kept)>0)
    X<-t(Z[kept,j[older],drop=FALSE]);Y<-t(Z[kept,j[younger],drop=FALSE])
    oa<-v3_nuisance(m[ia,,drop=FALSE],v3_P2);ob<-v3_nuisance(m[ib,,drop=FALSE],v3_P2)
    WA<-oa$Z[edge[,1],,drop=FALSE];WB<-ob$Z[edge[,2],,drop=FALSE]
    colnames(WA)<-paste0('older_',colnames(WA));colnames(WB)<-paste0('younger_',colnames(WB));W2<-cbind(WA,WB)
    fl_write_csv(data.frame(endpoint=c(rep('older',ncol(WA)),rep('younger',ncol(WB))),column=c(colnames(WA),colnames(WB))),
      file.path(out,'audit',paste0('P2_endpoint_design_',key,'.csv')))
    fl_write_csv(data.frame(endpoint=c('older','younger'),dropped=c(paste(oa$dropped,collapse=';'),paste(ob$dropped,collapse=';'))),
      file.path(out,'audit',paste0('P2_endpoint_dropped_terms_',key,'.csv')))
    samples[[key]]<-data.frame(comparison=key,participants=nrow(m),older_N=length(ia),younger_N=length(ib),
      families=length(unique(fam)),observed_within_family_pairs=nrow(edge),eligible_genera=sum(eligible))
    for(model in c('P1','P2')){
      W<-if(model=='P1')matrix(1,nrow(edge),1,dimnames=list(NULL,'Intercept'))else W2
      v3_log('Genus ',key,' ',model,': ',length(kept),' eligible genera; ',nrow(edge),' observed pairs')
      fit<-tf_cluster_regression(X,Y,W,fam);tf_verify_gee(X,Y,fit,fam,key,model,out)
      at<-match(ids,kept);b<-fit$beta[at];se<-fit$SE[at];valid<-eligible & !is.na(at) & tf_bool(fit$valid[at])
      results[[paste(key,model)]]<-data.frame(comparison=key,model=model,feature_id=ids,
        genus=sub('^g__','',ids),eligible=eligible,older_prevalence=pvA,younger_prevalence=pvB,
        participants=nrow(m),families=length(unique(fam)),observed_pairs=nrow(edge),
        beta_log2_per_log2=b,SE=se,CI_low=b-qnorm(.975)*se,CI_high=b+qnorm(.975)*se,
        P_raw=fit$P[at],analysis_valid=valid,status=ifelse(!eligible,'FILTERED_PREVALENCE',ifelse(valid,'ESTIMATED','NOT_ESTIMABLE')),
        nuisance_rank=ncol(fit$W),predictor_residual_SS=fit$predictor_residual_SS[at],
        younger_ratio_per_older_doubling=2^b,ratio_CI_low=2^(b-qnorm(.975)*se),ratio_CI_high=2^(b+qnorm(.975)*se))
    }
  }
  t<-ov_bind(results);t$P_BH<-NA_real_;t$FDR_tests<-0L
  for(model in c('P1','P2')){
    k<-which(t$model==model & t$eligible);good<-t$analysis_valid[k]&is.finite(t$P_raw[k])
    p<-rep(1,length(k));p[good]<-t$P_raw[k][good];q<-p.adjust(p,'BH')
    t$P_BH[k[good]]<-q[good];t$FDR_tests[t$model==model]<-length(k)
  }
  t$FDR_family<-paste0(t$model,'_all_eligible_genus_by_three_comparison_tests')
  t$FDR_lt_005<-is.finite(t$P_BH)&t$P_BH<.05
  summary<-ov_bind(lapply(split(t,list(t$model,t$comparison),drop=TRUE),function(z)data.frame(
    model=z$model[1],comparison=z$comparison[1],eligible_genera=sum(z$eligible),estimable=sum(z$analysis_valid),
    BH_positive=sum(z$FDR_lt_005 & z$beta_log2_per_log2>0,na.rm=TRUE),
    BH_negative=sum(z$FDR_lt_005 & z$beta_log2_per_log2<0,na.rm=TRUE),
    raw_P_lt_005=sum(z$analysis_valid & z$P_raw<.05,na.rm=TRUE),BH_family_tests=z$FDR_tests[1])))
  samples<-ov_bind(samples)
  for(v in c('participants','older_N','younger_N','families','observed_within_family_pairs'))
    stopifnot(identical(samples[[v]],old_counts[[v]][match(samples$comparison,old_counts$comparison)]))
  v3_table(t,'All_genus_cross_generation_GEE_results_EN',out)
  v3_table(samples,'Genus_pair_sample_sizes_EN',out);v3_table(summary,'Genus_concordance_summary_EN',out)
  v3_table(t[t$model=='P2' & t$FDR_lt_005,],'P2_BH_supported_genus_associations_EN',out)
  list(results=t,samples=samples,summary=summary)
}

gg_build<-function(cfg,out){
  p24<-lr_adopted(cfg$project_root,'24_Two_generation_reanalysis')
  p26<-lr_adopted(cfg$project_root,'26_Taxon_family_sharing')
  source<-file.path(p24,'private/stage_prepare.rds')
  paths<-c(source,file.path(p26,'tables/Family_pair_sample_sizes_EN.csv'),
    file.path(p26,'private',paste0('Observed_pairs_',gp_pairs,'.rds')),
    file.path(cfg$project_root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),
    file.path(cfg$project_root,'outputs/26_Taxon_family_sharing/LATEST_RUN.txt'))
  protected<-lr_hashes(paths);x<-readRDS(source)
  stopifnot(identical(x$protocol,v3_protocol),!anyDuplicated(x$metadata$sample_id),
    identical(as.character(x$metadata$sample_id),colnames(x$relative)),all(is.finite(x$relative)),
    all(x$relative>=0),max(abs(colSums(x$relative)-1))<1e-10,
    all(rowSums(table(x$metadata$family_id,x$metadata$generation)>0)>=2))
  a<-gg_aggregate(x,out);fit<-gg_fit(a,p26,cfg,out)
  # Save only aggregated source values in the portable figure checkpoint.
  r<-c(a[c('mapping','taxonomy','profiles','coverage','selected','counts')],fit,
    list(protocol=gg_protocol,protected=protected,baseline_N=nrow(x$metadata),
      families=length(unique(x$metadata$family_id))))
  gg_source_tables(r,out);lr_check_hashes(protected);r
}

gg_source_tables<-function(r,out){
  for(n in c('mapping','taxonomy','profiles','coverage','selected','counts'))v3_table(r[[n]],paste0('Genus_',n,'_EN'),out)
  v3_table(r$results,'All_genus_cross_generation_GEE_results_EN',out)
  v3_table(r$samples,'Genus_pair_sample_sizes_EN',out);v3_table(r$summary,'Genus_concordance_summary_EN',out)
  v3_table(r$results[r$results$model=='P2' & r$results$FDR_lt_005,],'P2_BH_supported_genus_associations_EN',out)
  d<-merge(r$profiles[r$profiles$feature_id%in%r$selected$feature_id,],
    r$selected[,c('feature_id','genus','phylum','family','abundance_rank','display_row')],by='feature_id')
  v3_table(d,'Figure36_generation_profile_source_EN',out)
  t<-r$results[r$results$model=='P2' & r$results$feature_id%in%r$selected$feature_id,]
  v3_table(t,'Figure36_adjusted_association_source_EN',out)
}

gg_geometry<-function(r){
  s<-r$selected;s$y<-(nrow(s)+1-s$display_row)*4
  p<-aggregate(y~phylum,s,mean);f<-aggregate(y~phylum+family+family_key,s,mean)
  root<-mean(p$y);edges<-list(data.frame(x=1,y=min(p$y),xend=1,yend=max(p$y)))
  edges[[2]]<-data.frame(x=1,y=p$y,xend=5,yend=p$y)
  for(i in seq_len(nrow(p))){j<-which(f$phylum==p$phylum[i]);yy<-f$y[j]
    edges[[length(edges)+1L]]<-data.frame(x=27,y=p$y[i],xend=30,yend=p$y[i])
    edges[[length(edges)+1L]]<-data.frame(x=30,y=min(c(yy,p$y[i])),xend=30,yend=max(c(yy,p$y[i])))
    edges[[length(edges)+1L]]<-data.frame(x=30,y=yy,xend=33,yend=yy)
  }
  for(i in seq_len(nrow(f))){yy<-s$y[s$family_key==f$family_key[i]]
    edges[[length(edges)+1L]]<-data.frame(x=62,y=f$y[i],xend=65,yend=f$y[i])
    edges[[length(edges)+1L]]<-data.frame(x=65,y=min(yy),xend=65,yend=max(yy))
    edges[[length(edges)+1L]]<-data.frame(x=65,y=yy,xend=68,yend=yy)
  }
  list(selected=s,phyla=p,families=f,edges=ov_bind(edges))
}

gg_heat_data<-function(r,cfg,g){
  ids<-g$selected$feature_id;tracks<-list()
  for(i in seq_len(6)){
    gen<-c('G0','G1','G2')[(i-1)%%3+1];p<-r$profiles[r$profiles$generation==gen,];at<-match(ids,p$feature_id)
    value<-if(i<=3)p$prevalence_percent[at]else p$mean_RA_percent[at]
    tracks[[i]]<-data.frame(feature_id=ids,y=g$selected$y,track=i,value=value,q=NA_real_)
  }
  for(i in 7:9){t<-r$results[r$results$model=='P2' & r$results$comparison==gp_pairs[i-6],]
    at<-match(ids,t$feature_id);v<-t$beta_log2_per_log2[at]
    v[!t$analysis_valid[at]]<-NA_real_
    tracks[[i]]<-data.frame(feature_id=ids,y=g$selected$y,track=i,value=v,q=t$P_BH[at])
  }
  z<-ov_bind(tracks);z$x<-c(113,120,127,137,144,151,162,169,176)[z$track]
  ra<-z$value[z$track%in%4:6 & z$value>0]
  lim<-c(floor(min(log10(ra))*2)/2,ceiling(max(log10(ra))*2)/2);if(diff(lim)==0)lim<-lim+c(-.5,.5)
  effect<-z$value[z$track>=7 & is.finite(z$value)]
  el<-if(length(effect))max(.1,ceiling(max(abs(effect))/.05)*.05)else .1
  z$fill<-'#E9E9E9';z$star<-''
  pal<-grDevices::colorRamp(c(cfg$binary_colors['Reference'],'#FFFFFF',cfg$binary_colors['Comparison']))
  for(i in 1:9){j<-which(z$track==i)
    if(i<=6){v<-if(i<=3)z$value[j]/100 else gp_ra_fraction(z$value[j],lim)
      z$fill[j]<-gp_tint(v,cfg$generation_colors[c('G0','G1','G2')[(i-1)%%3+1]])
    }else{ok<-j[is.finite(z$value[j])];z$fill[ok]<-grDevices::rgb(pal((z$value[ok]+el)/(2*el))/255)
      z$star[j[is.finite(z$q[j]) & z$q[j]<.05]]<-'*'
    }
  }
  list(data=z,RA_log_limits=lim,effect_limits=c(-el,el))
}

gg_legend<-function(p,r,cfg,h){
  text<-function(x,y,label,size=5.7,hjust=.5){p<<-gp_text(p,x,y,label,size,hjust)}
  rect<-function(x0,x1,y0,y1,col){p<<-gp_rect(p,x0,x1,y0,y1,col)}
  centers<-c(64,111,158);width<-33
  for(j in 1:3){gen<-c('G0','G1','G2')[j]
    text(centers[j],-3,gen,5.7)
    for(a in seq(0,1,length.out=41)){
      rect(centers[j]-width/2+a*width,centers[j]-width/2+min(1,a+1/40)*width,-9,-7,gp_tint(a,cfg$generation_colors[gen]))
      rect(centers[j]-width/2+a*width,centers[j]-width/2+min(1,a+1/40)*width,-21,-19,gp_tint(a,cfg$generation_colors[gen]))
    }
    for(a in c(0,.5,1))text(centers[j]-width/2+a*width,-12,as.character(a*100),5.3)
    breaks<-10^seq(ceiling(h$RA_log_limits[1]),floor(h$RA_log_limits[2]))
    for(v in breaks)text(centers[j]-width/2+(log10(v)-h$RA_log_limits[1])/diff(h$RA_log_limits)*width,-24,
      format(v,scientific=FALSE,trim=TRUE),5.3)
  }
  text(2,-8,'Detection (%)',5.7,0);text(2,-20,'Mean RA (%)',5.7,0)
  pal<-grDevices::colorRamp(c(cfg$binary_colors['Reference'],'#FFFFFF',cfg$binary_colors['Comparison']))
  for(a in seq(0,1,length.out=61))rect(64+a*94,64+min(1,a+1/60)*94,-33,-30,
    grDevices::rgb(pal(a)/255))
  for(v in c(h$effect_limits[1],0,h$effect_limits[2]))text(64+(v-h$effect_limits[1])/diff(h$effect_limits)*94,
    -36,format(round(2^v,2),nsmall=2),5.3)
  text(2,-31.5,'Adjusted ratio',5.7,0)
  text(64,-42,'* BH q < 0.05',5.3,0)
  rect(155,158,-43,-41,'#E9E9E9');text(160,-42,'Not tested / NE',5.3,0)
  p
}

gg_plot<-function(r,cfg,out){
  g<-gg_geometry(r);h<-gg_heat_data(r,cfg,g);d<-h$data
  v3_table(d,'Figure36_cell_values_and_colours_EN',out)
  v3_table(data.frame(scale=c('Mean RA log10 percent lower','Mean RA log10 percent upper',
    'P2 log2 ratio lower','P2 log2 ratio upper'),value=c(h$RA_log_limits,h$effect_limits)),
    'Figure36_shared_colour_limits_EN',out)
  p<-ggplot2::ggplot()+
    ggplot2::geom_segment(data=g$edges,ggplot2::aes(x,y,xend=xend,yend=yend),colour='#808080',linewidth=.25)+
    ggplot2::geom_rect(data=g$phyla,ggplot2::aes(xmin=5,xmax=27,ymin=y-1.3,ymax=y+1.3),
      fill='white',colour='#B5B5B5',linewidth=.2)+
    ggplot2::geom_rect(data=g$families,ggplot2::aes(xmin=33,xmax=62,ymin=y-1.3,ymax=y+1.3),
      fill='white',colour='#B5B5B5',linewidth=.2)+
    ggplot2::geom_text(data=g$phyla,ggplot2::aes(x=16,y=y,label=phylum),size=5.2/(72.27/25.4),family='Arial')+
    ggplot2::geom_text(data=g$families,ggplot2::aes(x=47.5,y=y,label=family),size=5.2/(72.27/25.4),family='Arial')+
    ggplot2::geom_text(data=g$selected,ggplot2::aes(x=70,y=y,label=genus,
      fontface=ifelse(unnamed_genome_genus_bin,'plain','italic')),hjust=0,size=6/(72.27/25.4),family='Arial')+
    ggplot2::geom_tile(data=d,ggplot2::aes(x,y,fill=fill),width=6.4,height=3.5)+
    ggplot2::geom_text(data=d[d$star!='',],ggplot2::aes(x,y,label=star),size=6/(72.27/25.4),family='Arial')+
    ggplot2::scale_fill_identity()+
    ggplot2::coord_cartesian(xlim=c(-1,181),ylim=c(-46,139),expand=FALSE,clip='off')+
    ggplot2::labs(title=NULL)+
    ggplot2::theme_void(base_family='Arial')+ggplot2::theme(plot.margin=ggplot2::margin(7,5,5,5,'pt'))
  for(k in 1:3)p<-gp_text(p,c(16,47.5,87)[k],132,c('Phylum','Family','Genus')[k],6)
  for(k in 1:3)p<-gp_text(p,c(120,144,169)[k],136,c('Detection (%)','Mean RA (%)','Adjusted association')[k],5.7)
  for(i in 1:9)p<-gp_text(p,c(113,120,127,137,144,151,162,169,176)[i],129,
    if(i<=6)c('G0','G1','G2')[(i-1)%%3+1]else gp_pairs[i-6],5.5,angle=90,
    colour=if(i<=6)cfg$generation_colors[c('G0','G1','G2')[(i-1)%%3+1]]else 'black')
  gg_legend(p,r,cfg,h)
}

gg_reports<-function(r,out){
  writeLines(c('# Figure 36. Genus-level profiles and cross-generation abundance associations.',
    sprintf('The 30 most abundant genera or genome genus bins are displayed from %s baseline participants (G0, 1,130; G1, 1,028; G2, 1,059). Selection uses overall mean relative abundance, including zeros, independently of association results.',format(r$baseline_N,big.mark=',')),
    'Left, phylum-family-genus classification relationships; branch lengths have no evolutionary or temporal interpretation. Formal genus names are italicized; GGB/FGB identifiers remain explicitly unnamed genome bins. Firmicutes and Bacillota are nomenclatural synonyms and are displayed consistently as Bacillota. Ambiguous higher taxonomy is labelled unresolved rather than inferred.',
    'Middle, detection prevalence (>0) and mean relative abundance (%) in each generation. Within-genus SGB abundances were summed for each participant before calculating either quantity; genus detection is not a sum of SGB detection rates. No re-normalization to the assigned genus fraction was performed. Each quantity uses one common colour intensity scale across generations; the mean-RA colour scale is log10 for positive values and zero is white.',
    'Right, newly fitted P2 Gaussian GEE associations for G0-G1, G1-G2 and G0-G2. Colours encode the adjusted younger-member abundance ratio per doubling of older-member abundance, on a symmetric log2-ratio colour scale centred at ratio 1; * denotes BH q<0.05. Grey denotes not tested (prevalence filter) or non-estimable (NE). All eligible genera, not only the displayed 30, enter the multiple-testing family. Complete results, including confidence intervals, are provided in the source tables.',
    'P2 adjusts both endpoints\' age splines, sex, preparation batch, sequencing pool, sampling season sine/cosine and recent antibiotic use. Primary-family cluster-robust standard errors account for repeated participants and pairs within a primary label. Labels do not verify direct parent-child relationships or cohabitation. These are abundance associations, not strain transmission, heritability or absolute microbial load.'),
    file.path(out,'reports/Figure_legend_EN.md'),useBytes=TRUE)
  writeLines(c('# Methods, references and study adaptations',
    'Aggregation: use every baseline SGB with exactly one non-unknown genus annotation; sum full-profile relative abundances per person. Ambiguous/unclassified mass is retained as an explicit coverage deficit rather than redistributed. GGB bins are not renamed to species or known genera. All aggregate genera are screened in each comparison; no significance or original-SGB filter is inherited.',
    'Eligibility: the genus must be detected in >=10% of distinct older participants and >=10% of distinct younger participants in the P2-complete comparison sample. The adopted baseline, endpoint covariates and exact observed pair lists match module 26. P1 uses the same complete sample with only an intercept; P2 adds both endpoint adjustment sets. Each observed pair is equally weighted; people and families may recur.',
    'For genus g, zero abundance is replaced by half its smallest positive full-baseline genus abundance. Fit log2(younger genus abundance) ~ log2(older genus abundance) + covariates. Gaussian independence GEE is estimated by the existing audited cluster-score kernel and checked against actual geepack::geeglm for the first three estimable genera in each comparison/model, with maximum absolute coefficient/SE/P discrepancy <1e-6.',
    'BH pools all eligible genus-by-three-comparison tests separately for P1 and P2. Non-estimable eligible tests conservatively count as P=1 in the BH denominator; their reported q remains unavailable. Original SGB coefficients or P values are never averaged or reused as genus inference.',
    'Genus aggregation and family-related microbiome questions are motivated by Valles-Colomer et al. Figure 1g. The specific older-to-younger exposure, YFS covariates, primary labels and regression formula are study adaptations of the module-26 framework; the paper did not use this exact formula. The new tree is a taxonomic hierarchy, not a reconstruction or collapse of a sequence-based maximum-likelihood phylogeny.',
    'Nomenclatural harmonization: exact phylum labels Firmicutes and Bacillota denote homotypic synonyms and are displayed as Bacillota, according to LPSN. No suffix groups, genus identities, abundance values or statistical results are merged by this display harmonization. https://lpsn.dsmz.de/phylum/bacillota',
    '',
    'Valles-Colomer M, Bacigalupe R, Vieira-Silva S, et al. Variation and transmission of the human gut microbiota across multiple familial generations. Nature Microbiology. 2022;7:87-96. doi:10.1038/s41564-021-01021-8. https://www.nature.com/articles/s41564-021-01021-8',
    'Liang KY, Zeger SL. Longitudinal data analysis using generalized linear models. Biometrika. 1986;73(1):13-22. doi:10.1093/biomet/73.1.13.',
    'Halekoh U, Hojsgaard S, Yan J. The R Package geepack for Generalized Estimating Equations. Journal of Statistical Software. 2006;15(2):1-11. doi:10.18637/jss.v015.i02. https://www.jstatsoft.org/article/view/v015i02'),
    file.path(out,'reports/Methods_and_references_EN.md'),useBytes=TRUE)
  z<-r$summary[r$summary$model=='P2',]
  lines<-c('# 36号：属水平代际分析',
    sprintf('基线%d人、%d个主家系。%d个SGB中，%d个能明确归入%d个属/基因组属分组；%d个SGB存在歧义或未分类。',
      r$baseline_N,r$families,r$counts$SGBs_baseline,r$counts$SGBs_assigned,r$counts$genus_groups,r$counts$SGBs_ambiguous_or_unclassified),
    sprintf('每人可明确归属到属的丰度比例中位数为%.1f%%；没有重新归一化。主图30属合计占原完整丰度的平均%.1f%%。',
      r$coverage$median_classified_percent[r$coverage$generation=='Overall'],r$counts$mean_total_RA_percent_in_display),
    '主图按总体平均丰度选择前30个属，不先筛显著。GGB是未正式命名的基因组属分组，不是具体菌株；无法可靠确定的科/门标注Unresolved。分类图连线长度不表示进化距离。',
    '原注释中的Firmicutes与Bacillota为同义门名，依据LPSN统一显示为Bacillota；不改变属的聚合、丰度或任何统计估计。',
    '属丰度逐人求和后重新拟合P1/P2，并重新对全部合格属×三组比较做BH；未沿用或平均SGB效应和P值。原始数据及24—35号结果保留。',
    'P2调整配对双方年龄样条、性别、建库批次、测序池、季节正余弦和近期抗生素，并采用主家系聚类稳健标准误。具体配对逐行核对与26号完全一致。',
    paste0(z$comparison,'：检验',z$eligible_genera,'个属，BH q<0.05正关联',z$BH_positive,'个，负关联',z$BH_negative,'个。'),
    '图中星号表示全属检验范围的BH q<0.05。颜色为年长成员丰度每翻倍时年轻成员丰度比；不是OR或菌株传播率。',
    '描述性检出率/丰度差异尚未调整年龄，不能直接解释为代际或年龄的独立效应。')
  writeLines(lines,file.path(out,'reports/Results_summary_zh.md'),useBytes=TRUE)
}

gg_preflight<-function(cfg,out){
  fun<-c('gg_plot','gg_legend','gg_geometry','gg_heat_data','gg_figures','gp_text','gp_rect','gp_tint',
    'gp_ra_fraction','lr_export','tf_export','nf_theme','nf_skill_root')
  path<-file.path(out,'audit/actual_plot_and_export_source.R')
  writeLines(unlist(lapply(fun,function(n)c(paste0(n,' <- '),deparse(get(n,mode='function'),width.cutoff=120L),''))),path)
  invisible(parse(path,encoding='UTF-8'))
  nf_qa_command(cfg,'validate_figure.py',c(shQuote(path),'--json'),file.path(out,'audit/plot_source_preflight.json'))
}
gg_figures<-function(r,cfg,out){
  p<-gg_plot(r,cfg,out)
  a<-lr_export(list(p),'Figure36_Genus_taxonomy_and_generation_profiles',out,cfg,width=183,height=185)
  v3_table(a,'Figure_export_manifest_EN',out)
}

gg_main<-function(args=character(),root=Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')){
  if(getRversion()<'4.5'){
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',
      shQuote(file.path(root,'scripts/36_genus_generation_taxonomy.R')),vapply(args,shQuote,character(1))))
    if(status!=0L)stop('Project-R genus-generation entry failed');return(invisible(NULL))
  }
  if(any(!args%in%c('--run','--check','--help','--figures-only','--resume-figures')))stop('Unknown mode')
  if(sum(c('--run','--check','--figures-only','--resume-figures')%in%args)>1L)stop('Choose one execution mode')
  if('--help'%in%args){cat('36_genus_generation_taxonomy.R [--run | --check | --figures-only | --resume-figures]\n');return(invisible(NULL))}
  cfg<-gg_setup(root);key<-'generation_genus';cfg$output_module<-unname(cfg$output_modules[key])
  cfg$run_id<-unname(cfg$run_ids[key]);cfg$run_id_explicit<-FALSE
  if('--check'%in%args){invisible(parse(file.path(root,'R/genus_generation_taxonomy.R')));cat('PASS: syntax/dependencies; no data or fits loaded\n');return(invisible(NULL))}
  module<-file.path(root,'outputs',cfg$output_module)
  if('--figures-only'%in%args){
    out<-lr_adopted(root,cfg$output_module);r<-readRDS(file.path(out,'private/analysis_result.rds'));lr_check_hashes(r$protected)
    stamp<-format(Sys.time(),'%Y%m%d_%H%M%S');dir.create(file.path(out,'audit',stamp))
    file.copy(list.files(file.path(out,'figures'),full.names=TRUE),file.path(out,'audit',stamp))
    file.copy(file.path(out,'private/analysis_result.rds'),file.path(out,'audit',stamp,'analysis_result_before_display_update.rds'))
    fitted_table<-file.path(out,'tables/All_genus_cross_generation_GEE_results_EN.csv')
    fit_md5<-unname(tools::md5sum(fitted_table))
    r<-gg_harmonize_taxonomy(r,out);gg_source_tables(r,out);gg_archive_table_aliases(out)
    stopifnot(identical(fit_md5,unname(tools::md5sum(fitted_table))))
    saveRDS(r,file.path(out,'private/analysis_result.rds'))
    gg_preflight(cfg,out);gg_figures(r,cfg,out);gg_reports(r,out);lr_check_hashes(r$protected)
    sys.source(file.path(root,'R/result_gallery.R'),envir=.GlobalEnv);fl_publish_organized_results(cfg);return(invisible(out))
  }
  resume<-NULL
  if('--resume-figures'%in%args){
    paths<-list.dirs(file.path(module,'history'),recursive=FALSE,full.names=TRUE)
    paths<-paths[file.exists(file.path(paths,'private/analysis_result.rds'))]
    if(!length(paths))stop('No saved genus analysis to recover')
    resume<-paths[which.max(file.info(file.path(paths,'private/analysis_result.rds'))$mtime)]
  }
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(p in c('tables','figures','reports','audit','private'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
  writeLines(c('RUNNING',gg_protocol),file.path(out,'run_status.txt'))
  tryCatch({
    if(is.null(resume))r<-gg_build(cfg,out)else{
      r<-readRDS(file.path(resume,'private/analysis_result.rds'));lr_check_hashes(r$protected)
      for(p in c('private','audit')){
        dest<-if(p=='private')file.path(out,p)else file.path(out,p,'saved_fit_audit')
        dir.create(dest,showWarnings=FALSE)
        files<-list.files(file.path(resume,p),full.names=TRUE)
        if(length(files)&&!all(file.copy(files,dest,recursive=TRUE)))stop('Saved-result recovery failed')
      }
      gg_source_tables(r,out);writeLines('Saved genus results recovered; no regression refit.',file.path(out,'audit/recovery_mode.txt'))
    }
    r<-gg_harmonize_taxonomy(r,out);gg_source_tables(r,out);gg_archive_table_aliases(out)
    saveRDS(r,file.path(out,'private/analysis_result.rds'));gg_reports(r,out);gg_preflight(cfg,out)
    gg_figures(r,cfg,out);lr_check_hashes(r$protected)
    fl_write_csv(r$protected,file.path(out,'audit/protected_source_files.csv'))
    writeLines(capture.output(sessionInfo()),file.path(out,'audit/sessionInfo.txt'))
    writeLines('PASS: sources unchanged; exact module-26 pairs; genus GEE verified against geepack; rendered QA passed.',file.path(out,'audit/delivery_validation.txt'))
    writeLines(c('COMPLETED',gg_protocol),file.path(out,'run_status.txt'))
    target<-fl_output_complete_run(ctx);writeLines(basename(target),file.path(module,'LATEST_RUN.txt'))
    v3_log('COMPLETED: ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
