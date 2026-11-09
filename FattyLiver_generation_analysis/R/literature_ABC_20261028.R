# Inert definitions only. The three independent numbered entries call lr_main().
# Written at the user's request; no analysis or plotting runs while sourcing
# this support file. Original module-29 inference is reused, never overwritten.
lr_protocol<-'literature_ABC_redesign_20261028_v1'
lr_modules<-c(A='literature_within',B='literature_severity',C='literature_metabolic')
lr_entries<-c(A='32_A_within_family_literature.R',B='33_B_steatosis_severity.R',C='34_C_metabolic_phenotype_map.R')
lr_models<-c('B1','Adult_B1_same_sample','B2')

lr_setup<-function(root){
  # RStudio can pass an older user library to the project-R subprocess.
  # Prefer this R's own Matrix before any inherited library (no reinstall).
  .libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),file.path(R.home(),'library')))
  sys.source(file.path(root,'R/abc_analysis.R'),envir=.GlobalEnv)
  cfg<-abc_setup(root)
  sys.source(file.path(root,'R/result_gallery.R'),envir=.GlobalEnv)
  for(p in c('ggplot2','patchwork','svglite','jsonlite','digest','Maaslin2','lme4','lmerTest'))
    if(!requireNamespace(p,quietly=TRUE))stop('Missing project dependency: ',p)
  matrix_abi<-as.character(Matrix::Matrix.Version()$abi)
  lme4_abi<-as.character(get('.Matrix.abi.build.version',asNamespace('lme4')))
  if(!identical(matrix_abi,lme4_abi))stop('Matrix/lme4 ABI mismatch: run in the clean project R 4.5.2 session')
  cfg
}
lr_adopted<-function(root,module){
  base<-file.path(root,'outputs',module);id<-readLines(file.path(base,'LATEST_RUN.txt'),n=1,warn=FALSE)
  stopifnot(length(id)==1L,grepl('^[A-Za-z0-9_.-]+$',id))
  path<-file.path(base,'current',id)
  if(!file.exists(file.path(path,'run_status.txt'))||readLines(file.path(path,'run_status.txt'),1)!='COMPLETED')stop('Completed adopted run required: ',module)
  path
}
lr_table<-function(x,name,out)v3_table(x,name,out)
lr_theme<-function()tf_plot_theme()+ggplot2::theme(legend.text=ggplot2::element_text(size=6),legend.title=ggplot2::element_text(size=6.5))
lr_empty<-function(title,note){
  ggplot2::ggplot(data.frame(x=.5,y=.5),ggplot2::aes(x,y))+
    ggplot2::geom_text(label=note,size=6/(72.27/25.4),family='Arial')+
    ggplot2::xlim(0,1)+ggplot2::ylim(0,1)+ggplot2::labs(title=title,x=NULL,y=NULL)+ggplot2::theme_void(base_size=6.5)
}
lr_stars<-function(q){
  ifelse(!is.finite(q),'',ifelse(q<.001,'****',ifelse(q<.01,'***',ifelse(q<.05,'**',ifelse(q<.1,'*','')))))
}
lr_export<-function(plots,stem,out,cfg,width=183,height=130,nrow=1,design=NULL,exemptions=list(),uniform_y_axes=FALSE){
  # Apply the common theme first; retain panel-specific blank axes, legend
  # positions and hidden duplicate legends. Previously these were overwritten.
  plots<-lapply(seq_along(plots),function(i){
    local_theme<-plots[[i]]$theme
    plots[[i]]+lr_theme()+local_theme+ggplot2::theme(
      text=ggplot2::element_text(family='Arial'),
      plot.title=ggplot2::element_text(size=6.8),
      plot.tag=ggplot2::element_text(size=8,face='bold'),
      legend.text=ggplot2::element_text(size=6),legend.title=ggplot2::element_text(size=6.5))+
      ggplot2::labs(tag=letters[i])
  })
  p<-if(is.null(design))patchwork::wrap_plots(plots,nrow=nrow,guides='keep')else patchwork::wrap_plots(plots,design=design,guides='keep')
  tf_export(list(plot=p,width=width,height=height,panels=length(plots),exemptions=exemptions,uniform_y_axes=uniform_y_axes),stem,out,cfg)
}
lr_label<-function(ids)sub('^t__','',ids)
lr_hashes<-function(paths)data.frame(path=normalizePath(paths,winslash='/',mustWork=TRUE),md5=unname(tools::md5sum(paths)))
lr_check_hashes<-function(h){
  if(any(!file.exists(h$path))||!identical(unname(tools::md5sum(h$path)),h$md5))stop('Protected source files changed; publication stopped')
  invisible(TRUE)
}

# ----- A: exact reuse of completed LDM, literature-style presentation -----
lr_A<-function(cfg,out){
  parent<-lr_adopted(cfg$project_root,'29_Within_family_steatosis')
  paths<-file.path(parent,c('private/analysis_result.rds','private/B1_analysis_design.rds',
    'private/Adult_B1_same_sample_analysis_design.rds','private/B2_analysis_design.rds','audit/input_fingerprints.csv'))
  stopifnot(all(file.exists(paths)));protected<-lr_hashes(paths)
  recorded<-tf_read(tail(paths,1));lr_check_hashes(recorded)
  src<-abc_source(cfg$project_root)
  if(!any(normalizePath(recorded$path,winslash='/',mustWork=TRUE)==normalizePath(src$file,winslash='/',mustWork=TRUE)&recorded$md5==src$md5))stop('Module 29 and the current source population differ; no stale inference reuse')
  r<-readRDS(paths[1]);t<-r$results
  stopifnot(all(c('model','feature_id','P_raw','P_BH','beta_percentage_points','scale_direction_agrees')%in%names(t)))
  for(model in unique(t$model)){
    j<-t$model==model
    if(max(abs(p.adjust(t$P_raw[j],'BH')-t$P_BH[j]))>1e-12)stop('Saved BH correction mismatch: ',model)
  }
  designs<-setNames(lapply(paths[2:4],readRDS),lr_models)
  adult<-designs$Adult_B1_same_sample$metadata;adult2<-designs$B2$metadata
  stopifnot(identical(adult$sample_id,adult2$sample_id),all(adult$age>=18))
  a<-t[t$model=='Adult_B1_same_sample',];b<-t[t$model=='B2',]
  stopifnot(setequal(a$feature_id,b$feature_id))
  t$source_style_q010<-is.finite(t$P_BH)&t$P_BH<.1
  t$strict_q005<-is.finite(t$P_BH)&t$P_BH<.05
  lr_table(t,'A_All_saved_within_family_SGB_tests_EN',out)
  s<-r$summary;s$source_style_q010<-vapply(s$model,function(k)if(!any(t$model==k))NA_integer_ else sum(t$source_style_q010[t$model==k]),integer(1))
  lr_table(s,'A_Within_family_model_summary_EN',out)
  ids1<-a$feature_id[is.finite(a$P_BH)&a$P_BH<.1];ids2<-b$feature_id[is.finite(b$P_BH)&b$P_BH<.1]
  overlap<-data.frame(comparison='Same adult participants: B1 versus B2',adult_B1_only=length(setdiff(ids1,ids2)),
    common=length(intersect(ids1,ids2)),adult_B2_only=length(setdiff(ids2,ids1)),threshold='BH q < 0.10; source-style exploratory')
  lr_table(overlap,'A_Same_adult_model_overlap_EN',out)
  z<-t[t$model%in%lr_models,];ids<-unique(z$feature_id[z$source_style_q010])
  rank<-data.frame(feature_id=ids)
  for(k in lr_models)rank[[k]]<-z$P_BH[match(paste(ids,k),paste(z$feature_id,z$model))]
  # One display-only order, frozen before drawing; missing estimates sort last.
  rank<-rank[order(rank$B2,rank$Adult_B1_same_sample,rank$B1,rank$feature_id,na.last=TRUE),,drop=FALSE]
  rank$main_display<-seq_len(nrow(rank))<=30L
  rank$selection_rule<-rep('Union q<0.10; B2 q, adult-B1 q, all-age-B1 q, ID; first 30 main, all rows supplemental',nrow(rank))
  lr_table(rank,'A_Heatmap_display_selection_EN',out)
  common<-intersect(ids1,ids2);rep<-b[b$feature_id%in%common,]
  rep<-head(rep[order(rep$P_BH,rep$P_raw,rep$feature_id),],2L)
  paired<-lapply(rep$feature_id,function(id){
    y<-as.numeric(src$x$relative[id,match(adult$sample_id,colnames(src$x$relative))])
    d<-data.frame(family_id=as.character(adult$family_id),case=adult$fatty_liver,abundance=100*y)
    m<-aggregate(abundance~family_id+case,d,mean)
    stopifnot(all(table(m$family_id)==2L))
    m$feature_id<-id;m
  })
  paired<-ov_bind(paired)
  saveRDS(paired,file.path(out,'private/A_family_mean_plot_data.rds'))
  lr_table(rep,'A_Representative_shared_SGBs_EN',out)
  fl_write_csv(protected,file.path(out,'audit/A_Reused_source_files_EN.csv'))
  lr_check_hashes(protected);lr_check_hashes(recorded)
  list(results=t,summary=s,overlap=overlap,rank=rank,representatives=rep,paired=paired,protected=protected,
    provenance=list(method='Reuse completed official LDM; no refit',parent=parent,source=src$file))
}
lr_venn<-function(r,cfg){
  theta<-seq(0,2*pi,length.out=241);pal<-unname(cfg$binary_colors)
  circles<-rbind(data.frame(x=-.55+cos(theta),y=sin(theta),set='Adults: B1'),data.frame(x=.55+cos(theta),y=sin(theta),set='Adults: B2'))
  counts<-data.frame(x=c(-.88,0,.88),y=0,n=unname(unlist(r$overlap[1,c('adult_B1_only','common','adult_B2_only')])))
  ggplot2::ggplot(circles,ggplot2::aes(x,y,group=set))+
    ggplot2::geom_polygon(ggplot2::aes(fill=set),alpha=.28,colour=NA)+
    ggplot2::geom_path(ggplot2::aes(colour=set),linewidth=.5)+
    ggplot2::geom_text(data=counts,ggplot2::aes(x,y,label=n),inherit.aes=FALSE,size=6.5/(72.27/25.4))+
    ggplot2::scale_fill_manual(values=setNames(pal,c('Adults: B1','Adults: B2')))+
    ggplot2::scale_colour_manual(values=setNames(pal,c('Adults: B1','Adults: B2')))+
    ggplot2::coord_cartesian(xlim=c(-1.7,1.7),ylim=c(-1.2,1.2))+
    ggplot2::labs(title='Same-adult model overlap',subtitle='BH q < 0.10',x=NULL,y=NULL)+
    ggplot2::theme_void(base_size=6.5)+ggplot2::theme(legend.position='bottom',legend.title=ggplot2::element_blank())
}
lr_A_heatmap<-function(r,ids,cfg){
  if(!length(ids))return(lr_empty('Within-family differences','No SGB with BH q < 0.10'))
  d<-expand.grid(feature_id=ids,model=lr_models,stringsAsFactors=FALSE)
  at<-match(paste(d$feature_id,d$model),paste(r$results$feature_id,r$results$model))
  d$effect<-r$results$beta_percentage_points[at];d$q<-r$results$P_BH[at]
  d$label<-paste0(lr_stars(d$q),ifelse(!is.na(at)&!r$results$scale_direction_agrees[at],'\u2020',''))
  d$SGB<-factor(lr_label(d$feature_id),levels=rev(lr_label(ids)))
  d$model<-factor(d$model,lr_models,c('All ages\nB1','Adults\nB1','Adults\nB2'))
  lim<-max(abs(d$effect),na.rm=TRUE);if(!is.finite(lim)||lim==0)lim<-1
  ggplot2::ggplot(d,ggplot2::aes(model,SGB,fill=effect))+
    ggplot2::geom_tile(colour='white',linewidth=.12)+
    ggplot2::geom_text(ggplot2::aes(label=label),size=5.2/(72.27/25.4))+
    ggplot2::scale_fill_gradient2(low=unname(cfg$binary_colors['Reference']),mid='white',high=unname(cfg$binary_colors['Comparison']),
      midpoint=0,limits=c(-lim,lim),na.value='#D5D5D5',name='Adjusted difference (percentage points)',
      breaks=c(-lim,0,lim),labels=scales::label_number(accuracy=10^(floor(log10(lim))-1)),
      guide=ggplot2::guide_colourbar(display='rectangles',direction='horizontal',title.position='top',label.position='bottom',ticks=FALSE,
        barwidth=grid::unit(28,'mm'),barheight=grid::unit(2.6,'mm'),
        theme=ggplot2::theme(legend.text=ggplot2::element_text(margin=ggplot2::margin(t=3,unit='pt')),
          legend.title=ggplot2::element_text(margin=ggplot2::margin(b=3,unit='pt')))))+
    ggplot2::scale_y_discrete(position='right')+ggplot2::labs(title='Within-family abundance differences',x=NULL,y=NULL)+
    ggplot2::theme(axis.ticks=ggplot2::element_blank(),axis.line=ggplot2::element_blank(),legend.position='top',
      legend.key.width=grid::unit(10,'mm'),axis.text.y=ggplot2::element_text(size=5.7))
}
lr_A_pairplot<-function(r,id,cfg){
  if(is.na(id)||!nzchar(id))return(lr_empty('Family-group abundances','No shared adult B1/B2 hit'))
  d<-r$paired[r$paired$feature_id==id,];d$group<-factor(d$case,0:1,c('No steatosis','Steatosis'))
  q<-r$representatives$P_BH[match(id,r$representatives$feature_id)]
  ggplot2::ggplot(d,ggplot2::aes(group,abundance))+
    ggplot2::geom_line(ggplot2::aes(group=family_id),colour='#999999',alpha=.25,linewidth=.18)+
    ggplot2::geom_boxplot(ggplot2::aes(fill=group),colour='black',width=.45,outlier.shape=NA,linewidth=.3)+
    ggplot2::geom_point(colour='#777777',size=.3,alpha=.35)+
    ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),levels(d$group)))+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.03,.14)),breaks=scales::breaks_pretty(3))+
    ggplot2::labs(title=lr_label(id),subtitle=sprintf('Individual-level adult B2: BH q = %.3g',q),x=NULL,y='Family-mean relative abundance (%)')+
    ggplot2::theme(legend.position='none',plot.subtitle=ggplot2::element_text(size=5.7))
}
lr_A_figures<-function(r,cfg,out){
  ids<-head(r$rank$feature_id,30);reps<-c(r$representatives$feature_id,rep('',2))[1:2]
  plots<-list(lr_venn(r,cfg),lr_A_heatmap(r,ids,cfg),lr_A_pairplot(r,reps[1],cfg),lr_A_pairplot(r,reps[2],cfg))
  # Match the second left-column gutter to the first, which also contains the
  # Venn legend. The final-device alignment audit verifies this physical gap.
  plots[[4]]<-plots[[4]]+ggplot2::theme(plot.margin=ggplot2::margin(t=12.6,r=7,b=5,l=5,unit='pt'))
  # The heatmap intentionally spans the three left-hand evidence rows.
  exported<-list(lr_export(plots,'Figure32_Within_family_Chen_style',out,cfg,height=170,design='AB\nCB\nDB'))
  if(length(r$rank$feature_id)>30){
    pages<-split(r$rank$feature_id,ceiling(seq_len(nrow(r$rank))/30L))
    for(i in seq_along(pages))exported[[length(exported)+1]]<-lr_export(list(lr_A_heatmap(r,pages[[i]],cfg)),
      sprintf('Supplement32_All_supported_SGBs_%02d',i),out,cfg,width=120,height=max(75,35+length(pages[[i]])*3.5))
  }
  lr_table(ov_bind(exported),'Figure_export_manifest_EN',out)
  writeLines(c('# Figure 32. Within-primary-family steatosis associations.',
    'Visual structure follows Chen et al. Cell Reports Medicine 2024, Figure 2A-C: set overlap, differential-abundance heatmap and connected abundance boxplots. Disease, SGB resolution, comparison columns and unequal family sizes are study adaptations.',
    '(a) Non-area-proportional Venn counts compare adult B1 and B2 in exactly the same people. This is model overlap, not active/quiescent disease, enrichment or independent validation.',
    '(b) Individual-level LDM frequency coefficients, steatosis minus no steatosis, in percentage points. Stars use calibrated LDM-omnibus P values followed by BH: * q<0.10, ** q<0.05, *** q<0.01, **** q<0.001. A dagger indicates frequency/AST direction disagreement. Missing estimates are grey. The display union and deterministic 30-row rule are tabulated; supplemental pages retain all supported SGBs.',
    '(c-d) Up to two shared adult B1/B2 hits ranked by B2 q, raw P and ID. Grey lines connect the case and non-case means within a primary family; boxes show median/IQR and whiskers 1.5 IQR across family means. Inference is the saved individual-level B2 LDM, not a paired test of these plotted means.',
    'q<0.10 is source-style exploratory; all original q<0.05 flags are also retained. Primary labels do not verify siblings, cohabitation or genetics. Sources: doi:10.1016/j.xcrm.2024.101624; doi:10.1186/s40168-021-01034-9.'),
    file.path(out,'reports/Figure_legends_EN.md'),useBytes=TRUE)
}

# ----- B: official MaAsLin2, observed ultrasound grade -----
lr_B_design<-function(x,m,fields,kind){
  m<-droplevels(m[order(m$family_id,m$sample_id),]);nz<-v3_nuisance(m,fields);m<-nz$metadata
  Z<-nz$Z[,colnames(nz$Z)!='(Intercept)',drop=FALSE]
  centers<-colMeans(Z);ss<-apply(Z,2,sd);stopifnot(all(is.finite(ss)&ss>0))
  map<-data.frame(column=sprintf('nuis_%03d',seq_len(ncol(Z))),original=colnames(Z),center=centers,scale=ss)
  Z<-sweep(sweep(Z,2,centers,'-'),2,ss,'/');colnames(Z)<-map$column
  base<-if(kind=='ordinal')'grade_score'else'grade'
  md<-data.frame(family_id=factor(m$family_id),grade_score=as.numeric(m$fat_code-1),
    grade=factor(m$fat_code,1:3,c('No','Mild','Clear')),Z,row.names=m$sample_id,check.names=FALSE)
  M<-model.matrix(reformulate(c(base,map$column)),md)
  if(qr(M)$rank!=ncol(M)||nrow(M)<=ncol(M)+10)stop('Non-identifiable grade design')
  A<-t(x$relative[,match(m$sample_id,colnames(x$relative)),drop=FALSE]);pv<-colMeans(A>=.0002)
  list(metadata=m,design=md,relative=A,bases=base,terms=if(kind=='ordinal')'grade_score'else c('gradeMild','gradeClear'),
    nuisance=map$column,mapping=map,random=any(table(m$family_id)>1),
    filter=data.frame(feature_id=colnames(A),prevalence_at_RA_00002=pv,retained=pv>=.10),dropped=nz$dropped)
}
lr_B_extract<-function(ans,ids,d,expected,model,kind){
  rows<-diagnostics<-list()
  for(id in ids){
    fit<-ans$fits[[id]];err<-'';messages<-character();singular<-NA;hmin<-NA_real_;opt<-0;cs<-NULL;response_error<-NA_real_
    if(inherits(fit,'merMod')){
      messages<-as.character(fit@optinfo$conv$lme4$messages);opt<-unlist(fit@optinfo$conv$opt)
      H<-fit@optinfo$derivs$Hessian;if(is.matrix(H)&&all(is.finite(H)))hmin<-min(eigen(H,symmetric=TRUE,only.values=TRUE)$values)
      singular<-lme4::isSingular(fit,tol=1e-4)
      data<-model.frame(fit);formula<-formula(fit);na.action<-stats::na.fail
      ft<-tryCatch(withCallingHandlers(lmerTest::as_lmerModLmerTest(fit),warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);NULL})
      if(!is.null(ft))cs<-coef(summary(ft))
    }else if(inherits(fit,'lm'))cs<-coef(summary(fit))else err<-'Official fit unavailable'
    if(!is.null(cs)){
      mf<-model.frame(fit);at<-match(rownames(mf),rownames(expected))
      if(anyNA(at))stop('Official response sample-order mismatch')
      response_error<-max(abs(as.numeric(model.response(mf))-expected[at,id]))
      if(response_error>1e-10)stop('Official LOG abundance response mismatch')
    }
    conv<-messages[!grepl('boundary.*singular|singular fit',messages,ignore.case=TRUE)]
    bad<-nzchar(err)||length(conv)>0||any(opt!=0)||(!is.na(hmin)&&hmin< -1e-6)
    for(term in d$terms){
      beta<-se<-p<-df<-lo<-hi<-NA_real_;finite<-matched<-FALSE
      if(!is.null(cs)&&term%in%rownames(cs)){
        beta<-cs[term,'Estimate'];se<-cs[term,'Std. Error'];df<-if('df'%in%colnames(cs))cs[term,'df']else df.residual(fit)
        pcol<-grep('^Pr',colnames(cs),value=TRUE)[1];p<-cs[term,pcol]
        finite<-all(is.finite(c(beta,se,p,df)))&&se>0&&df>0&&p>=0&&p<=1
        if(finite){lo<-beta-qt(.975,df)*se;hi<-beta+qt(.975,df)*se}
        off<-ans$results[ans$results$feature==id&ans$results$metadata==d$bases,,drop=FALSE]
        if(kind=='categorical')off<-off[as.character(off$value)==sub('^grade','',term),,drop=FALSE]
        if(nrow(off)==1L&&finite)matched<-isTRUE(all.equal(unname(beta),off$coef[1],tolerance=1e-7))&&isTRUE(all.equal(unname(p),off$pval[1],tolerance=1e-7))
      }
      status<-if(!finite)'NOT_ESTIMABLE'else if(!matched)'OFFICIAL_RESULT_MISMATCH'else if(bad)'UNRELIABLE_CONVERGENCE'else if(isTRUE(singular))'ESTIMATED_SINGULAR'else'ESTIMATED'
      rows[[paste(id,term)]]<-data.frame(model=model,kind=kind,feature_id=id,term=term,beta_log2=beta,SE=se,df=df,
        CI_low=lo,CI_high=hi,P_raw=p,analysis_valid=status%in%c('ESTIMATED','ESTIMATED_SINGULAR'),status=status,singular=singular)
    }
    diagnostics[[id]]<-data.frame(model=model,kind=kind,feature_id=id,optimizer=paste(opt,collapse=';'),
      minimum_Hessian_eigenvalue=hmin,singular=singular,response_max_error=response_error,warnings=paste(unique(messages),collapse='; '),error=err)
  }
  list(effects=ov_bind(rows),diagnostics=ov_bind(diagnostics))
}
lr_B_fit<-function(x,m,fields,model,kind,cfg,out,cache){
  prefix<-paste(model,kind,sep='_');counts<-table(factor(m$fat_code,1:3))
  if(nrow(m)<30||any(counts==0))return(list(results=NULL,status=data.frame(model=model,kind=kind,status='NOT_ESTIMABLE: missing grade or N<30')))
  d<-tryCatch(lr_B_design(x,m,fields,kind),error=function(e)e)
  if(inherits(d,'error'))return(list(results=NULL,status=data.frame(model=model,kind=kind,status=conditionMessage(d))))
  ids<-d$filter$feature_id[d$filter$retained]
  lr_table(d$filter,paste0('B_Filter_',prefix,'_EN'),out)
  fl_write_csv(d$mapping,file.path(out,'audit',paste0(prefix,'_design.csv')))
  saveRDS(d$metadata,file.path(out,'private',paste0(prefix,'_metadata.rds')))
  if(!length(ids))return(list(results=NULL,status=data.frame(model=model,kind=kind,status='NOT_ESTIMABLE: no eligible species')))
  signature<-digest::digest(list(lr_protocol,model,kind,d$relative[,ids,drop=FALSE],d$design,
    Maaslin2=as.character(packageVersion('Maaslin2')),extract=body(lr_B_extract)),algo='sha256')
  batches<-split(ids,ceiling(seq_along(ids)/20L));parts<-list()
  for(i in seq_along(batches)){
    keys<-batches[[i]];folder<-file.path(cache,prefix,substr(signature,1,18),sprintf('batch_%03d',i))
    parts[[i]]<-tf_checkpoint(file.path(folder,'checkpoint.rds'),paste(signature,paste(keys,collapse=';')),function(){
      dir.create(folder,recursive=TRUE,showWarnings=FALSE);tf_log('B ',prefix,': batch ',i,'/',length(batches))
      A<-d$relative[,keys,drop=FALSE]
      expected<-vapply(keys,function(k)get('LOG',asNamespace('Maaslin2'))(A[,k]),numeric(nrow(A)));dimnames(expected)<-dimnames(A)
      con<-file(file.path(folder,'console.log'),'wt');sink(con)
      ans<-tryCatch(Maaslin2::Maaslin2(input_data=as.data.frame(A),input_metadata=d$design,output=file.path(folder,'official'),
        fixed_effects=c(d$bases,d$nuisance),random_effects=if(d$random)'family_id'else NULL,
        reference=if(kind=='categorical')'grade,No'else NULL,
        normalization='NONE',transform='LOG',standardize=FALSE,min_prevalence=0,min_abundance=0,min_variance=0,
        analysis_method='LM',cores=1,correction='BH',plot_heatmap=FALSE,plot_scatter=FALSE,save_models=TRUE),finally={sink();close(con)})
      lr_B_extract(ans,keys,d,expected,model,kind)
    })
  }
  estimated<-ov_bind(lapply(parts,`[[`,'effects'));diagnostics<-ov_bind(lapply(parts,`[[`,'diagnostics'))
  fl_write_csv(diagnostics,file.path(out,'audit',paste0(prefix,'_diagnostics.csv')))
  all<-expand.grid(feature_id=d$filter$feature_id,term=d$terms,stringsAsFactors=FALSE)
  idx<-match(paste(all$feature_id,all$term),paste(estimated$feature_id,estimated$term))
  t<-estimated[idx,,drop=FALSE];t$feature_id<-all$feature_id;t$term<-all$term;t$model<-model;t$kind<-kind
  t$eligible<-d$filter$retained[match(t$feature_id,d$filter$feature_id)]
  t$status[!t$eligible]<-'FILTERED_DETECTION_PREVALENCE';t$analysis_valid[is.na(t$analysis_valid)]<-FALSE
  good<-t$eligible&t$analysis_valid&is.finite(t$P_raw);p<-rep(1,sum(t$eligible));p[good[t$eligible]]<-t$P_raw[good]
  q<-p.adjust(p,'BH');t$P_BH<-NA_real_;t$P_BH[which(t$eligible)[good[t$eligible]]]<-q[good[t$eligible]]
  t$source_style_q010<-is.finite(t$P_BH)&t$P_BH<.1;t$strict_q005<-is.finite(t$P_BH)&t$P_BH<.05
  t$FDR_family<-paste0(model,'_',kind,'_all_eligible_SGBs_and_contrasts')
  t$N<-nrow(m);t$families<-length(unique(m$family_id));t$No_N<-unname(counts[1]);t$Mild_N<-unname(counts[2]);t$Clear_N<-unname(counts[3])
  t$SGB_id<-lr_label(t$feature_id);t$species_label<-ov_species_label(x$features$taxonomy[match(t$feature_id,x$features$feature_id)])
  t$abundance_ratio<-2^t$beta_log2;t$ratio_CI_low<-2^t$CI_low;t$ratio_CI_high<-2^t$CI_high
  t$effect_unit<-if(kind=='ordinal')'Zero-replaced relative-abundance ratio per one ordered grade step'else'Zero-replaced relative-abundance ratio versus No steatosis'
  list(results=t,status=data.frame(model=model,kind=kind,status='FITTED; see per-SGB diagnostics'))
}
lr_B<-function(cfg,out,cache){
  src<-tf_source(cfg$project_root);protected<-data.frame(path=src$files,md5=unname(src$md5));x<-tf_load(src);m<-x$metadata
  stopifnot('fat_code'%in%names(m),all(m$fat_code%in%1:3),all(m$fatty_liver==as.integer(m$fat_code>1)))
  samples<-v3_micro_samples(x);rows<-fits<-list()
  for(model in lr_models){
    mm<-droplevels(m[samples[[model]]$mask,,drop=FALSE])
    for(code in 1:3)rows[[paste(model,code)]]<-data.frame(model=model,grade=c('No','Mild','Clear')[code],N=sum(mm$fat_code==code),
      families=length(unique(mm$family_id[mm$fat_code==code])),minors=sum(mm$age[mm$fat_code==code]<18))
    for(kind in c('ordinal','categorical'))fits[[paste(model,kind)]]<-lr_B_fit(x,mm,samples[[model]]$fields,model,kind,cfg,out,cache)
  }
  stopifnot(identical(which(samples$Adult_B1_same_sample$mask),which(samples$B2$mask)))
  t<-ov_bind(lapply(fits,`[[`,'results'));status<-ov_bind(lapply(fits,`[[`,'status'))
  lr_table(ov_bind(rows),'B_Grade_sample_sizes_EN',out);lr_table(status,'B_Model_status_EN',out)
  lr_table(t,'B_All_grade_associations_EN',out)
  if(nrow(t))lr_table(t[t$source_style_q010,],'B_Source_style_supported_SGBs_EN',out)
  saveRDS(list(metadata=m[samples$B1$mask,],relative=x$relative[,samples$B1$mask,drop=FALSE]),file.path(out,'private/B_Grade_distribution_data.rds'))
  lr_check_hashes(protected)
  list(results=t,summary=ov_bind(rows),status=status,protected=protected)
}
lr_B_figures<-function(r,cfg,out){
  t<-r$results;exports<-list()
  z<-if(nrow(t))t[t$model=='B1'&t$kind=='ordinal'&t$source_style_q010&t$analysis_valid,]else data.frame()
  if(nrow(z)){
    z<-head(z[order(z$P_BH,z$P_raw,z$feature_id),],12L);z$SGB<-factor(z$SGB_id,levels=rev(z$SGB_id))
    z$direction<-factor(ifelse(z$beta_log2<0,'Lower','Higher'),levels=c('Lower','Higher'))
    p<-ggplot2::ggplot(z,ggplot2::aes(abundance_ratio,SGB))+
      ggplot2::geom_vline(xintercept=1,linetype=2,colour='#999999',linewidth=.25)+
      ggplot2::geom_segment(ggplot2::aes(x=ratio_CI_low,xend=ratio_CI_high,yend=SGB),colour='black',linewidth=.35)+
      ggplot2::geom_point(ggplot2::aes(fill=direction),shape=21,colour='black',size=1.8,stroke=.25)+
      ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),c('Lower','Higher')))+
      ggplot2::scale_x_log10(breaks=scales::breaks_log(4))+
      ggplot2::labs(title='Ultrasound steatosis-grade trend',x='Adjusted abundance ratio per grade step (95% CI)',y=NULL)
  }else p<-lr_empty('Ultrasound steatosis-grade trend','No eligible B1 SGB with BH q < 0.10')
  exports[[1]]<-lr_export(list(p),'Figure33_Grade_trend_effects',out,cfg,width=120,height=110)
  if(nrow(z)){
    source<-readRDS(file.path(out,'private/B_Grade_distribution_data.rds'));ids<-head(z$feature_id,6L)
    plots<-lapply(ids,function(id){
      m<-source$metadata;d<-data.frame(grade=factor(m$fat_code,1:3,c('No','Mild','Clear')),
        abundance=100*as.numeric(source$relative[id,match(m$sample_id,colnames(source$relative))]))
      qq<-z$P_BH[match(id,z$feature_id)]
      ggplot2::ggplot(d,ggplot2::aes(grade,abundance,fill=grade))+
        ggplot2::geom_boxplot(colour='black',linewidth=.3,width=.5,outlier.shape=NA)+
        ggplot2::geom_point(position=ggplot2::position_jitter(width=.12,height=0,seed=33),colour='#777777',size=.3,alpha=.2)+
        ggplot2::scale_fill_manual(values=c(No=unname(cfg$binary_colors['Reference']),Mild='#FBD1A7',Clear=unname(cfg$binary_colors['Comparison'])))+
        ggplot2::labs(title=lr_label(id),subtitle=sprintf('Ordinal Model 3: BH q = %.3g',qq),x='Observed ultrasound grade',y='Relative abundance (%)')+
        ggplot2::theme(legend.position='none',plot.subtitle=ggplot2::element_text(size=5.7))
    })
    exports[[2]]<-lr_export(plots,'Figure33_Grade_abundance_profiles',out,cfg,height=if(length(plots)>3)145 else 90,nrow=if(length(plots)>3)2 else 1,uniform_y_axes=TRUE)
    lr_table(z,'B_Displayed_ordinal_SGBs_EN',out)
  }
  lr_table(ov_bind(exports),'Figure_export_manifest_EN',out)
  writeLines(c('# Figure 33. Microbial abundance across observed ultrasound steatosis grades.',
    'Official MaAsLin2 LOG models follow the ordinal/categorical disease-state framework of Mei et al. Nature Medicine 2024 (doi:10.1038/s41591-024-03067-7). The single-cohort, disease, family intercept and B1/B2 covariates are explicitly documented applications; no independent-cohort meta-analysis is claimed.',
    'No, Mild and Clear map exactly to original ultrasound codes 1, 2 and 3. Grade steps are ordinal, not measured equal increments in liver-fat concentration. The adjusted ratio is 2^beta on zero-replaced relative abundance. Ordinary 95% confidence intervals are not multiplicity-adjusted.',
    'Source-style exploratory discovery is BH q<0.10; q<0.05 flags are also tabulated. All eligible SGBs are tested before display selection. Up to 12 top ordinal B1 hits are shown in the effect plot, with raw distributions for up to six. Boxes show median/IQR, whiskers 1.5 IQR, grey dots individuals. These raw distributions are unadjusted; their annotations come from the adjusted ordinal model.',
    'B1 includes generation, generation-specific age splines, sex, batch, pool, season and antibiotics with primary-family dependence. Adult common-sample B2 adds smoking, alcohol, diet, physical activity, energy and BMI. Cross-sectional grade differences do not demonstrate progression, fibrosis, causation or transmission.'),file.path(out,'reports/Figure_legends_EN.md'),useBytes=TRUE)
}

# ----- C: source-style partial Spearman, outcome-blind one adult per label -----
lr_partial_spearman<-function(x,y,controls){
  keep<-is.finite(x)&is.finite(y)&apply(as.matrix(controls),1,function(v)all(is.finite(v)))
  x<-x[keep];y<-y[keep];controls<-as.matrix(controls[keep,,drop=FALSE]);n<-length(x)
  if(n<10||length(unique(x))<2||length(unique(y))<2)return(data.frame(N=n,df=NA_real_,rho=NA_real_,P_raw=NA_real_,status='NOT_ESTIMABLE: small N or constant variable'))
  rankcols<-vapply(seq_len(ncol(controls)),function(j)rank(controls[,j],ties.method='average'),numeric(n))
  Z<-cbind(Intercept=1,rankcols);q<-qr(Z,tol=1e-9);df<-n-q$rank-1
  if(q$rank!=ncol(Z)||df<=3)return(data.frame(N=n,df=df,rho=NA_real_,P_raw=NA_real_,status='NOT_ESTIMABLE: nuisance rank/df'))
  rx<-qr.resid(q,rank(x,ties.method='average'));ry<-qr.resid(q,rank(y,ties.method='average'))
  if(sd(rx)<1e-12||sd(ry)<1e-12)return(data.frame(N=n,df=df,rho=NA_real_,P_raw=NA_real_,status='NOT_ESTIMABLE: no residual rank variation'))
  rho<-max(-1,min(1,cor(rx,ry)));stat<-rho*sqrt(df/pmax(1-rho^2,.Machine$double.eps));p<-2*pt(-abs(stat),df)
  data.frame(N=n,df=df,rho=rho,P_raw=p,status='ESTIMATED')
}
lr_C<-function(cfg,out){
  src<-tf_source(cfg$project_root);protected<-data.frame(path=src$files,md5=unname(src$md5));x<-tf_load(src);m<-x$metadata
  fields<-c(BMI18='BMI',WAIST18='Waist',TRIGLY18='Triglycerides',HDLKOL18='HDL-C',ALAT18='ALT',ASAT18='AST')
  stopifnot(all(names(fields)%in%names(m)))
  eligible<-m$age>=18&v3_cc(m,c('age','BMI18','family_id'))&m$BMI18>0
  mm<-droplevels(m[eligible,]);mm<-mm[order(mm$family_id,mm$sample_id),]
  set.seed(20261028L);indices<-sort(vapply(split(seq_len(nrow(mm)),mm$family_id),function(ii)ii[sample.int(length(ii),1)],integer(1)))
  mm<-droplevels(mm[indices,]);stopifnot(!anyDuplicated(mm$family_id),!anyDuplicated(mm$sample_id))
  A<-x$relative[,match(mm$sample_id,colnames(x$relative)),drop=FALSE];pv<-rowMeans(A>0);tested<-pv>.20
  lr_table(data.frame(feature_id=rownames(A),prevalence=pv,eligible=tested,rule='Present in >20% of outcome-blind selected adults'),
    'C_Feature_filter_EN',out)
  saveRDS(list(metadata=mm,selection_seed=20261028L,selection='One observed age/BMI-complete adult per primary label before phenotype filtering'),
    file.path(out,'private/C_Selected_adult_population.rds'))
  flow<-data.frame(stage=c('Adopted microbiome participants','Observed age/BMI-complete adults','One adult per primary family'),
    N=c(nrow(m),sum(eligible),nrow(mm)))
  lr_table(flow,'C_One_per_family_sample_flow_EN',out)
  rows<-list();ids<-rownames(A)[tested]
  for(field in names(fields)){
    controls<-if(field=='BMI18')'age'else c('age','BMI18')
    y<-mm[[field]];valid<-if(field%in%c('ALAT18','ASAT18'))y>=0 else y>0;y[!is.finite(y)|is.na(valid)|!valid]<-NA_real_
    for(id in ids){
      result<-lr_partial_spearman(as.numeric(A[id,]),y,mm[,controls,drop=FALSE])
      rows[[paste(field,id)]]<-cbind(feature_id=id,phenotype=field,phenotype_label=unname(fields[field]),
        adjustment=paste(controls,collapse=' + '),result)
    }
  }
  t<-ov_bind(rows);if(!nrow(t))stop('No eligible SGB for the predeclared metabolic map')
  good<-t$status=='ESTIMATED'&is.finite(t$P_raw);p<-rep(1,nrow(t));p[good]<-t$P_raw[good]
  t$P_BH<-NA_real_;t$P_BH[good]<-p.adjust(p,'BH')[good]
  t$source_style_q020<-is.finite(t$P_BH)&t$P_BH<.2;t$strict_q005<-is.finite(t$P_BH)&t$P_BH<.05
  t$FDR_family<-'All eligible SGB x six predeclared phenotypes; source-style exploratory q<0.20'
  t$SGB_id<-lr_label(t$feature_id);t$species_label<-ov_species_label(x$features$taxonomy[match(t$feature_id,x$features$feature_id)])
  ranking<-aggregate(source_style_q020~feature_id,t,sum);names(ranking)[2]<-'supported_phenotypes'
  ranking<-ranking[order(-ranking$supported_phenotypes,ranking$feature_id),];ranking$display<-ranking$supported_phenotypes>0&seq_len(nrow(ranking))<=30L
  lr_table(t,'C_All_partial_Spearman_associations_EN',out);lr_table(ranking,'C_Top30_display_selection_EN',out)
  summary<-ov_bind(lapply(names(fields),function(field){z<-t[t$phenotype==field,];data.frame(phenotype=field,label=unname(fields[field]),
    tested=nrow(z),estimable=sum(z$status=='ESTIMATED'),positive_q020=sum(z$source_style_q020&z$rho>0,na.rm=TRUE),
    negative_q020=sum(z$source_style_q020&z$rho<0,na.rm=TRUE),strict_q005=sum(z$strict_q005))}))
  lr_table(summary,'C_Phenotype_association_summary_EN',out)
  lr_check_hashes(protected)
  list(results=t,ranking=ranking,summary=summary,flow=flow,protected=protected)
}
lr_C_order<-function(M){
  if(nrow(M)<2)return(list(order=seq_len(nrow(M)),status='Single row; no clustering required'))
  d<-dist(M,method='euclidean')
  if(any(!is.finite(d)))return(list(order=seq_len(nrow(M)),status='Clustering not identifiable from finite pairwise distances; deterministic unclustered order retained'))
  list(order=hclust(d,method='complete')$order,status='Complete-linkage Euclidean; pairwise finite coordinates; no zero imputation')
}
lr_C_figures<-function(r,cfg,out){
  ids<-r$ranking$feature_id[r$ranking$display];fields<-r$summary$phenotype
  if(length(ids)){
    d<-r$results[r$results$feature_id%in%ids,];M<-matrix(NA_real_,length(ids),length(fields),dimnames=list(ids,fields))
    M[cbind(match(d$feature_id,ids),match(d$phenotype,fields))]<-d$rho
    ro<-lr_C_order(M);co<-lr_C_order(t(M));ids<-ids[ro$order];fields<-fields[co$order]
    lr_table(data.frame(axis=c('Rows','Columns'),method=c(ro$status,co$status)),'C_Heatmap_clustering_EN',out)
    d$SGB<-factor(lr_label(d$feature_id),levels=rev(lr_label(ids)))
    d$trait<-factor(d$phenotype,levels=fields,labels=r$summary$label[match(fields,r$summary$phenotype)])
    d$star<-ifelse(d$source_style_q020,'*','')
    # One common, symmetric linear scale, rounded outward in 0.05 rho steps.
    # Preserve every value: do not clip cells or scale individual columns.
    colour_limit<-max(.05,ceiling(max(abs(d$rho),na.rm=TRUE)/.05)*.05)
    stopifnot(is.finite(colour_limit),all(abs(d$rho[is.finite(d$rho)])<=colour_limit))
    lr_table(data.frame(lower=-colour_limit,midpoint=0,upper=colour_limit,
      observed_min=min(d$rho,na.rm=TRUE),observed_max=max(d$rho,na.rm=TRUE),
      rule='Common symmetric linear rho scale; 0.05 outward rounding; no clipping or column standardization',
      low_colour=unname(cfg$binary_colors['Reference']),high_colour=unname(cfg$binary_colors['Comparison'])),
      'C_Heatmap_colour_scale_EN',out)
    p<-ggplot2::ggplot(d,ggplot2::aes(trait,SGB,fill=rho))+
      ggplot2::geom_tile(colour='white',linewidth=.15)+ggplot2::geom_text(ggplot2::aes(label=star),size=6/(72.27/25.4))+
      ggplot2::scale_fill_gradient2(low=unname(cfg$binary_colors['Reference']),mid='white',high=unname(cfg$binary_colors['Comparison']),
        midpoint=0,limits=c(-colour_limit,colour_limit),breaks=c(-colour_limit,0,colour_limit),
        labels=scales::label_number(accuracy=.01),na.value='#D5D5D5',name='Partial Spearman rho',
        guide=ggplot2::guide_colourbar(display='rectangles',direction='horizontal',title.position='top',label.position='bottom',ticks=FALSE,
          barwidth=grid::unit(28,'mm'),barheight=grid::unit(2.6,'mm'),
          theme=ggplot2::theme(legend.text=ggplot2::element_text(margin=ggplot2::margin(t=3,unit='pt')),
            legend.title=ggplot2::element_text(margin=ggplot2::margin(b=3,unit='pt')))))+
      ggplot2::labs(title='Microbiome-metabolic phenotype map',subtitle='* BH q < 0.20; source-style exploratory',x=NULL,y=NULL)+
      ggplot2::theme(axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=45,hjust=1),
        axis.text.y=ggplot2::element_text(size=5.8),legend.position='top',plot.subtitle=ggplot2::element_text(size=6))
    lr_table(d,'C_Displayed_heatmap_cells_EN',out)
  }else p<-lr_empty('Microbiome-metabolic phenotype map','No supported SGB-phenotype association at BH q < 0.20')
  record<-lr_export(list(p),'Figure34_Metabolic_phenotype_map',out,cfg,width=150,height=max(90,45+length(ids)*3.5))
  lr_table(record,'Figure_export_manifest_EN',out)
  writeLines(c('# Figure 34. Partial-Spearman map of microbial and metabolic phenotypes.',
    'Following Asnicar et al. Nature Medicine 2021 (doi:10.1038/s41591-020-01183-8), taxa present in >20% of selected adults are correlated with six predeclared phenotypes using partial Spearman, BH correction and complete-linkage Euclidean clustering.',
    'One age/BMI-complete adult per primary family is selected with seed 20261028 before inspecting phenotype-specific values or associations. This is a family-structure adaptation, not the original PREDICT sampling scheme. BMI correlations adjust age only; waist, triglycerides, HDL cholesterol, ALT and AST adjust age and BMI. Pairwise-complete N and df are tabulated.',
    'Colour encodes rho; one star denotes source-style exploratory BH q<0.20. Stricter q<0.05 flags are retained in the full table. The 30 taxa with the most supported phenotype associations are displayed, with ID tie-breaking. This controls display only; BH includes all eligible SGB-by-phenotype tests. Grey cells are unestimable, not zero. Clustering status is reported, and missing rho is never filled with zero.',
    'The colour scale is linear, symmetric around zero and common to all six phenotypes. Its limit is rounded outward in 0.05-rho steps to contain every displayed correlation; the plotted limits are shown on the colourbar and recorded in C_Heatmap_colour_scale_EN. No clipping, nonlinear colour transformation or column standardization is used.',
    'No HOMA, inflammatory-marker, mediation, causal or prospective-prediction analysis is performed. Primary labels do not verify unrelatedness across labels. This map describes selected observed clinical phenotypes, not all metabolic pathways or microbial functions.'),file.path(out,'reports/Figure_legends_EN.md'),useBytes=TRUE)
}

lr_reports<-function(r,which,cfg,out){
  file.copy(file.path(cfg$project_root,'analysis_plan/adopted_ABC_redesign_20261028.md'),file.path(out,'reports/Methods_and_references_EN.md'),overwrite=TRUE)
  writeLines(c(paste0('# ',which,' 结果摘要'),'',
    switch(which,A='家系内个人脂肪肝比较：复用已完成29号官方LDM，以Chen Figure 2A-C形式展示。',
      B='超声脂肪变三级程度：官方MaAsLin2有序趋势及分类比较。',C='成人代谢表型图谱：一名成人/主家系、偏Spearman、BH校正。'),
    '',capture.output(print(r$summary,row.names=FALSE)),'',
    '文献式探索阈值与既有BH q<0.05标记分别报告；不把阈值放宽后的结果包装成原有主要发现。',
    '结局仅为超声脂肪变，不是MASLD；主家系标签未核实直接亲缘、同住或遗传匹配。',
    '原始RDS与既往结果保留。方法、图注、阈值、展示筛选及本队列适配见英文报告和完整表。'),
    file.path(out,'reports/Results_summary_zh.md'),useBytes=TRUE)
  capture.output(sessionInfo(),file=file.path(out,'audit/sessionInfo.txt'))
  lr_table_paths<-c('Matrix','lme4')
  fl_write_csv(data.frame(package=lr_table_paths,version=vapply(lr_table_paths,function(p)as.character(packageVersion(p)),character(1)),
    path=vapply(lr_table_paths,find.package,character(1)),ABI=c(as.character(Matrix::Matrix.Version()$abi),
      as.character(get('.Matrix.abi.build.version',asNamespace('lme4'))))),file.path(out,'audit/Matrix_lme4_runtime.csv'))
}
lr_failed_saved_run<-function(module){
  candidates<-list.dirs(file.path(module,'history'),recursive=FALSE,full.names=TRUE)
  candidates<-candidates[vapply(candidates,function(p){
    status<-file.path(p,'run_status.txt')
    file.exists(status)&&identical(readLines(status,1,warn=FALSE),'FAILED')&&
      file.exists(file.path(p,'private/analysis_result.rds'))
  },logical(1))]
  if(!length(candidates))stop('No failed run with saved analysis results; use --run')
  candidates[which.max(file.info(file.path(candidates,'private/analysis_result.rds'))$mtime)]
}
lr_preflight<-function(cfg,out){
  # Audit the actual drawing/export functions together. Auditing only this
  # analytical support file would miss the shared Cairo/vector/alignment gates.
  functions<-c('lr_theme','lr_empty','lr_stars','lr_label','lr_export','lr_venn',
    'lr_A_heatmap','lr_A_pairplot','lr_A_figures','lr_B_figures','lr_C_order','lr_C_figures',
    'tf_export','nf_theme','nf_skill_root')
  text<-unlist(lapply(functions,function(n)c(paste0(n,' <- '),deparse(get(n,mode='function'),width.cutoff=120L),'')))
  path<-file.path(out,'audit/actual_plot_and_export_source.R');writeLines(text,path,useBytes=TRUE)
  invisible(parse(path,encoding='UTF-8'))
  nf_qa_command(cfg,'validate_figure.py',c(shQuote(path),'--json'),file.path(out,'audit/plot_source_preflight.json'))
}
lr_main<-function(which,args=character(),root=Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')){
  stopifnot(which%in%c('A','B','C'))
  matrix_loaded_elsewhere<-'Matrix'%in%loadedNamespaces()&&
    !identical(normalizePath(getNamespaceInfo(asNamespace('Matrix'),'path'),winslash='/'),
      normalizePath(file.path(root,'tools/R-4.5.2/library/Matrix'),winslash='/',mustWork=FALSE))
  if(getRversion()<'4.5'||matrix_loaded_elsewhere){
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts',lr_entries[[which]])),vapply(args,shQuote,character(1))))
    if(status!=0L)stop('Project-R literature analysis failed');return(invisible(NULL))
  }
  if(any(!args%in%c('--run','--check','--help','--figures-only','--resume-figures','--fresh-analysis','--resume-fresh-analysis')))stop('Unsupported execution mode')
  if(sum(c('--run','--check','--figures-only','--resume-figures','--fresh-analysis','--resume-fresh-analysis')%in%args)>1L)stop('Choose one execution mode')
  if(any(c('--fresh-analysis','--resume-fresh-analysis')%in%args)&&which!='B')stop('Fresh analysis modes are implemented only for the explicitly requested B/33 reanalysis')
  if('--help'%in%args){cat(lr_entries[[which]],'[--run | --check | --figures-only | --resume-figures | --fresh-analysis | --resume-fresh-analysis (B only)]\nOnly this route executes. No new pipeline.\n');return(invisible(NULL))}
  cfg<-lr_setup(root);key<-lr_modules[[which]];cfg$output_module<-unname(cfg$output_modules[key]);cfg$run_id<-unname(cfg$run_ids[key]);cfg$run_id_explicit<-FALSE
  if(anyNA(c(cfg$output_module,cfg$run_id)))stop('Register the literature module in config.R first')
  if('--check'%in%args){invisible(parse(file.path(root,'R/literature_ABC_20261028.R'),encoding='UTF-8'));cat('Static/runtime check only; no participant data or models loaded\n');return(invisible(NULL))}
  module<-file.path(root,'outputs',cfg$output_module)
  resume<-if('--resume-figures'%in%args)lr_failed_saved_run(module)else NULL
  fresh_resume <- NULL
  if('--resume-fresh-analysis'%in%args) {
    candidates <- list.dirs(file.path(module,'history'),recursive=FALSE,full.names=TRUE)
    candidates <- candidates[vapply(candidates,function(p)
      file.exists(file.path(p,'audit/analysis_execution_mode.txt'))&&
      length(list.files(file.path(p,'private'),pattern='checkpoint[.]rds$',recursive=TRUE))>0L,logical(1))]
    if(!length(candidates))stop('No checkpoints from this fresh reanalysis were found; old analysis caches will not be used')
    fresh_resume <- candidates[which.max(file.info(file.path(candidates,'audit/analysis_execution_mode.txt'))$mtime)]
  }
  draw<-switch(which,A=lr_A_figures,B=lr_B_figures,C=lr_C_figures)
  if('--figures-only'%in%args){
    out<-lr_adopted(root,cfg$output_module);r<-readRDS(file.path(out,'private/analysis_result.rds'));lr_check_hashes(r$protected)
    lr_preflight(cfg,out);draw(r,cfg,out);lr_check_hashes(r$protected)
    fl_publish_organized_results(cfg);return(invisible(out))
  }
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(p in c('tables','figures','reports','audit','private'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
  writeLines(c('RUNNING',lr_protocol,which),file.path(out,'run_status.txt'))
  tryCatch({
    fresh <- any(c('--fresh-analysis','--resume-fresh-analysis')%in%args)
    # Short path avoids Windows MAX_PATH in MaAsLin2's feature export names.
    cache<-if(fresh)file.path(out,'private/cp')else file.path(module,'_checkpoints/private')
    dir.create(cache,recursive=TRUE,showWarnings=FALSE)
    if(fresh) {
      writeLines(c('Fresh Model-3/all-age, Model-3/adult and Model-4/adult ordinal/categorical fits requested 2026-11-01.',
        'Historical model checkpoints are NOT reused. Only new checkpoints inside this run are allowed.',
        'Cohort, input abundance, filters, covariates, statistical method, BH families and display rules remain unchanged.',
        paste('Started:',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'))),file.path(out,'audit/analysis_execution_mode.txt'))
      tf_log('Fresh analysis: all six grade model sets will be fitted; historical checkpoints are not read')
      if(!is.null(fresh_resume)) {
        old_cache <- if(dir.exists(file.path(fresh_resume,'private/cp')))
          file.path(fresh_resume,'private/cp')else file.path(fresh_resume,'private/fresh_analysis_checkpoints')
        copied <- list();source_files <- list.files(old_cache,pattern='checkpoint[.]rds$',recursive=TRUE,full.names=TRUE)
        for(source_file in source_files) {
          relative <- substring(source_file,nchar(old_cache)+2L);target_file <- file.path(cache,relative)
          dir.create(dirname(target_file),recursive=TRUE,showWarnings=FALSE)
          stopifnot(file.copy(source_file,target_file),identical(unname(tools::md5sum(source_file)),unname(tools::md5sum(target_file))))
          copied[[length(copied)+1L]] <- data.frame(source=source_file,destination=target_file,
            original_mtime=format(file.info(source_file)$mtime,'%Y-%m-%d %H:%M:%S %Z'),
            md5=unname(tools::md5sum(source_file)))
        }
        fl_write_csv(ov_bind(copied),file.path(out,'audit/recovered_same_request_fresh_checkpoints.csv'))
        file.copy(file.path(fresh_resume,'audit/analysis_execution_mode.txt'),
          file.path(out,'audit/original_fresh_request_execution_mode.txt'))
        tf_log('Resuming ',length(source_files),' completed NEW checkpoints from the same 2026-11-01 request; 2026-10-28 caches remain unused')
      }
    }
    if(!is.null(resume)){
      r<-readRDS(file.path(resume,'private/analysis_result.rds'));lr_check_hashes(r$protected)
      originals<-unlist(lapply(c('tables','private'),function(p)list.files(file.path(resume,p),recursive=TRUE,full.names=TRUE)))
      fingerprints<-lr_hashes(originals)
      for(p in c('tables','private')){
        files<-list.files(file.path(resume,p),full.names=TRUE)
        if(length(files)&&!all(file.copy(files,file.path(out,p),recursive=TRUE,overwrite=FALSE)))stop('Saved-result recovery copy failed')
      }
      # Retain the original fitting runtime and per-feature convergence/design
      # records separately from the new redraw-only runtime and rendered QA.
      source_audit<-file.path(out,'audit','saved_analysis_audit')
      dir.create(source_audit,showWarnings=FALSE)
      audit_files<-list.files(file.path(resume,'audit'),full.names=TRUE)
      if(length(audit_files)&&!all(file.copy(audit_files,source_audit,recursive=TRUE,overwrite=FALSE)))stop('Original analysis audit recovery failed')
      fl_write_csv(fingerprints,file.path(out,'audit/recovered_result_fingerprints.csv'))
      writeLines(c('Figure recovery only; no data preparation or statistical model fit.',paste('Saved run:',resume)),
        file.path(out,'audit/recovery_mode.txt'))
      lr_check_hashes(fingerprints)
      tf_log('Redrawing saved ',which,' results; no statistical model refit')
    }else r<-switch(which,A=lr_A(cfg,out),B=lr_B(cfg,out,cache),C=lr_C(cfg,out))
    saveRDS(r,file.path(out,'private/analysis_result.rds'));lr_reports(r,which,cfg,out)
    # Static figure gate runs before the first graphics device; rendered gates
    # are inside lr_export and block COMPLETED publication if they fail.
    lr_preflight(cfg,out)
    draw(r,cfg,out);lr_check_hashes(r$protected)
    writeLines('PASS: source hashes unchanged; all enabled rendered QA gates passed.',file.path(out,'audit/delivery_validation.txt'))
    writeLines(c('COMPLETED',lr_protocol,which),file.path(out,'run_status.txt'));target<-fl_output_complete_run(ctx)
    writeLines(basename(target),file.path(module,'LATEST_RUN.txt'));tf_log('COMPLETED: ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
