# R-only figures reuse the adopted visual style and rendered QA/export helper.
abc_wrap <- function(plots,nrow=1) {
  plots<-lapply(seq_along(plots),function(i)plots[[i]]+tf_plot_theme()+ggplot2::labs(tag=letters[i])+
    ggplot2::theme(legend.position='none',axis.title.x=ggplot2::element_text(margin=ggplot2::margin(t=7,unit='pt')),
      axis.text.x=ggplot2::element_text(margin=ggplot2::margin(t=3,unit='pt'))))
  patchwork::wrap_plots(plots,nrow=nrow,guides='keep')
}
abc_sig_palette <- function(cfg) c('Not BH-supported'='#C6C6C6','Lower'=unname(cfg$binary_colors['Reference']),
  'Higher'=unname(cfg$binary_colors['Comparison']),'Scale-dependent'='#777777')
abc_box <- function(d,y,title,ylab,cfg) {
  d$group<-factor(d$exposure,0:1,c('No older case','Older case'))
  ggplot2::ggplot(d,ggplot2::aes(group,.data[[y]],fill=group))+
    ggplot2::geom_boxplot(width=.5,colour='black',linewidth=.3,outlier.shape=NA)+
    ggplot2::geom_point(position=ggplot2::position_jitter(width=.14,height=0,seed=28),colour='#777777',size=.3,alpha=.25)+
    ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),levels(d$group)))+
    ggplot2::labs(title=title,x=NULL,y=ylab)
}
abc_volcano <- function(t,x,title,xlab,cfg,p='P_raw',sig='BH_lt_005',direction=NULL) {
  if(is.null(direction))direction<-ifelse(t[[x]]<0,'Lower','Higher')
  t$significance<-ifelse(!t[[sig]],'Not BH-supported',direction)
  ggplot2::ggplot(t,ggplot2::aes(.data[[x]],-log10(pmax(.data[[p]],1e-300))))+
    ggplot2::geom_vline(xintercept=0,colour='#999999',linetype=2,linewidth=.25)+
    ggplot2::geom_point(ggplot2::aes(fill=significance),shape=21,colour='#555555',size=1.2,stroke=.12,alpha=.8)+
    ggplot2::scale_fill_manual(values=abc_sig_palette(cfg))+
    ggplot2::scale_x_continuous(breaks=scales::breaks_pretty(4))+
    ggplot2::scale_y_continuous(breaks=scales::breaks_pretty(4),expand=ggplot2::expansion(mult=c(0,.08)))+
    ggplot2::labs(title=title,x=xlab,y='-log10(raw P)')
}

abc_figures <- function(r,which,cfg,out) {
  manifest<-list();legends<-c('# Figure legends','')
  if(which=='A'){
    keys<-c('G0-G1','G1-G2','G0-G2')
    plots<-lapply(keys,function(k){m<-readRDS(file.path(out,'private',paste0(k,'_genus_shannon.rds')))
      abc_box(data.frame(exposure=m$older_any_case,shannon=m$shannon),'shannon',k,'Genus Shannon',cfg)})
    manifest[[1]]<-tf_export(list(plot=abc_wrap(plots),width=183,height=88,panels=3,exemptions=list()),'Figure28_Genus_Shannon',out,cfg)
    plots<-lapply(keys,function(k){z<-r$ordinations[[k]]
      ggplot2::ggplot(z,ggplot2::aes(PC1,PC2,fill=exposure))+
        ggplot2::geom_point(shape=21,colour='#555555',stroke=.12,size=1,alpha=.55)+
        ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),levels(z$exposure)))+
        ggplot2::scale_x_continuous(breaks=scales::breaks_pretty(3))+
        ggplot2::scale_y_continuous(breaks=scales::breaks_pretty(3))+
        ggplot2::labs(title=k,x='Genus Bray-Curtis PCoA 1',y='Genus Bray-Curtis PCoA 2')})
    manifest[[2]]<-tf_export(list(plot=abc_wrap(plots),width=183,height=88,panels=3,exemptions=list()),'Figure28_Genus_community',out,cfg)
    plots<-lapply(keys,function(k){t<-r$results[r$results$model=='B1'&r$results$comparison==k&r$results$analysis_valid,]
      t$BH_lt_005<-t$FDR_lt_005;abc_volcano(t,'beta_log2',k,'Adjusted coefficient (log2)',cfg)})
    manifest[[3]]<-tf_export(list(plot=abc_wrap(plots),width=183,height=88,panels=3,exemptions=list()),'Figure28_Genus_associations',out,cfg)
    legends<-c(legends,
      '**Figure28. Genus-level microbiome differences by observed older-generation steatosis among unaffected younger members.**',
      'Shannon: boxes show median and IQR, whiskers 1.5 IQR, grey points individual younger members. Blue: no observed older case; orange: at least one observed older case. Diversity uses the explicitly normalized subcomposition of unambiguously assigned genera/genome genus bins; it is not rarefied 16S diversity. Participants can recur across different generation comparisons.',
      'Community: descriptive Bray-Curtis PCoA (additive correction for non-Euclidean distances) on one seed-selected younger person per family. Colours as above. Original-style adonis2, adjusted dmanova and dispersion tests are reported in tables; the coordinates are unadjusted. No within-family permutation of the family-constant older exposure was used.',
      'Associations: each point is a genus tested using official MaAsLin2 LOG on full-profile relative abundance, B1 adjustment and family intercept where estimable. Blue/orange points are lower/higher with older-case exposure at joint BH q<0.05 across all eligible genera and three comparisons; grey points are not BH-supported. Genus and SGB results are different resolutions of the same data, not independent replications.')
  }
  if(which=='B'){
    t<-r$results[r$results$model=='B1',];direction<-ifelse(!t$scale_direction_agrees,'Scale-dependent',ifelse(t$beta_relative<0,'Lower','Higher'))
    plots<-list(abc_volcano(t,'beta_percentage_points','Frequency scale','Within-family difference (percentage points)',cfg,direction=direction),
      abc_volcano(t,'beta_AST','Arcsine-root scale','Within-family coefficient (AST)',cfg,direction=direction))
    manifest[[1]]<-tf_export(list(plot=abc_wrap(plots),width=183,height=105,panels=2,exemptions=list()),'Figure29_Within_family_LDM',out,cfg)
    s<-r$summary[is.finite(r$summary$BH_lt_005),]
    display<-c(B1='All ages: B1',G0_same_generation='G0 only: B1',Adult_B1_same_sample='Adults: B1',B2='Adults: B2')
    plotdata<-data.frame(model=factor(s$model,levels=rev(s$model),labels=unname(display[rev(s$model)])),count=s$BH_lt_005,N=s$N)
    p<-ggplot2::ggplot(plotdata,ggplot2::aes(count,model))+
      ggplot2::geom_col(width=.55,fill=unname(cfg$binary_colors['Comparison']),colour='black',linewidth=.25)+
      ggplot2::geom_text(ggplot2::aes(label=count),hjust=-.2,size=5.7/(72.27/25.4),family='Arial')+
      ggplot2::scale_x_continuous(breaks=scales::breaks_pretty(4),expand=ggplot2::expansion(mult=c(0,.1)))+
      ggplot2::labs(title='Within-family associations',x='SGBs with BH q < 0.05',y=NULL)
    pn<-ggplot2::ggplot(plotdata,ggplot2::aes(N,model))+
      ggplot2::geom_col(width=.55,fill=unname(cfg$binary_colors['Reference']),colour='black',linewidth=.25)+
      ggplot2::geom_text(ggplot2::aes(label=N),hjust=-.2,size=5.7/(72.27/25.4),family='Arial')+
      ggplot2::scale_x_continuous(breaks=scales::breaks_pretty(4),expand=ggplot2::expansion(mult=c(0,.1)))+
      ggplot2::labs(title='Analysis samples',x='Participants',y=NULL)
    manifest[[2]]<-tf_export(list(plot=abc_wrap(list(pn,p)),width=183,height=95,panels=2,exemptions=list()),'Figure29_Model_comparison',out,cfg)
    legends<-c(legends,
      '**Figure29. Within-primary-family associations of personal hepatic steatosis with microbial composition.**',
      'The two panels display the same official LDM-omnibus P values against coefficients on relative-frequency (percentage-point) and arcsine-root scales, respectively. Each participant contributes one row; only disease-discordant primary families are included. Nuisance terms include family indicators and B1 covariates; permutations are restricted within families. The omnibus combines the scales through the official calibrated LDM procedure. Blue/orange: BH-supported lower/higher abundance with concordant coefficient directions on both scales; dark grey: scale-dependent directions; light grey: not BH-supported. No confidence intervals or fold changes are implied by the permutation statistic.',
      'Model comparison: numbers of BH-supported SGBs in estimable primary, same-generation and adult same-sample B1/B2 analyses. Different samples and model sizes prevent interpreting bar-height differences as formal heterogeneity tests. G1 has no informative same-generation discordant sets. G2 full B1 has only one residual degree of freedom and is reported as not reliably estimable, rather than removing covariates to force a result. Primary labels do not verify siblings, cohabitation or genetic matching.')
  }
  if(which=='C'){
    m<-readRDS(file.path(out,'private/B1_readout_plot_data.rds'));m$group<-factor(m$fatty_liver,0:1,c('No steatosis','Steatosis'))
    vars<-c('Favourable_richness','Unfavourable_richness','Favourable_abundance','Unfavourable_abundance')
    nf<-sum(!is.na(r$mapping$feature_id[r$mapping$set=='Favourable']));nu<-sum(!is.na(r$mapping$feature_id[r$mapping$set=='Unfavourable']))
    titles<-c(sprintf('Favourable richness (%d/50 mapped)',nf),sprintf('Unfavourable richness (%d/50 mapped)',nu),
      sprintf('Favourable abundance (%d/50 mapped)',nf),sprintf('Unfavourable abundance (%d/50 mapped)',nu))
    plots<-lapply(seq_along(vars),function(i){v<-vars[i];d<-m;if(i>2)d[[v]]<-d[[v]]*100
      ggplot2::ggplot(d,ggplot2::aes(group,.data[[v]],fill=group))+
        ggplot2::geom_boxplot(width=.5,colour='black',linewidth=.3,outlier.shape=NA)+
        ggplot2::geom_point(position=ggplot2::position_jitter(width=.14,height=0,seed=30),colour='#777777',size=.25,alpha=.15)+
        ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),levels(d$group)))+
        ggplot2::labs(title=titles[i],x=NULL,y=if(i<=2)'Detected published-set SGBs'else'Cumulative relative abundance (%)')})
    manifest[[1]]<-tf_export(list(plot=abc_wrap(plots,nrow=2),width=183,height=150,panels=4,exemptions=list()),'Figure30_External_sets_distributions',out,cfg)
    t<-r$results[r$results$model%in%c('B1','B2'),];t$set<-ifelse(grepl('^Favourable',t$readout),'Favourable','Unfavourable')
    t$label<-factor(paste(t$set,t$model),levels=rev(c('Favourable B1','Favourable B2','Unfavourable B1','Unfavourable B2')))
    plots<-lapply(c(TRUE,FALSE),function(rich){z<-t[grepl('richness',t$readout)==rich,]
      ggplot2::ggplot(z,ggplot2::aes(display_beta,label))+
        ggplot2::geom_vline(xintercept=0,colour='#999999',linewidth=.25,linetype=2)+
        ggplot2::geom_segment(ggplot2::aes(x=display_CI_low,xend=display_CI_high,yend=label),colour='black',linewidth=.35)+
        ggplot2::geom_point(ggplot2::aes(shape=model,fill=set),colour='black',size=1.9,stroke=.25)+
        ggplot2::scale_shape_manual(values=c(B1=21,B2=23))+
        ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),c('Favourable','Unfavourable')))+
        ggplot2::scale_x_continuous(breaks=scales::breaks_pretty(4))+
        ggplot2::labs(title=if(rich)'Richness contrasts'else'Cumulative abundance association',
          x=if(rich)'Adjusted SGB count difference (95% CI)'else'Steatosis probability difference (pp / 0.1 AST)',y=NULL)})
    manifest[[2]]<-tf_export(list(plot=abc_wrap(plots),width=183,height=100,panels=2,exemptions=list()),'Figure30_External_sets_associations',out,cfg)
    legends<-c(legends,
      '**Figure30. Externally predefined health-ranked SGB sets and personal hepatic steatosis.**',
      'Sets are the 50 lowest and 50 highest HEALTH_ranks from Asnicar et al. Nature 2026, frozen before YFS association testing and cross-checked against Supplementary Table 7. Exact SGB/group mappings and missing members are tabulated; absent IDs are not replaced by related species. These labels denote published associations, not universally beneficial or harmful organisms.',
      'Distributions: adopted B1 participants. Boxes show median/IQR, whiskers 1.5 IQR, grey dots individuals; blue no steatosis, orange steatosis. Counts use detected published members; cumulative abundance retains the full original microbiome denominator. Raw distributions are unadjusted.',
      'Associations: circles B1, diamonds adult B2; lines ordinary 95% robust family-GEE confidence intervals. Left: count outcome, steatosis predictor. Right: binary steatosis outcome, AST cumulative abundance predictor, preserving the original disease-model direction and Gaussian linear framework. Effects are percentage-point probability differences per 0.1 AST unit, not ORs; a linear probability model may yield predictions outside 0-1. Richness outcome direction and primary-family GEE are explicitly documented extensions. B1 includes age, sex, generation, batch, pool, season and antibiotics; B2 adds adult lifestyle and BMI. OLS source-style and same-adult-sample B1 results are in tables. Four readouts per model are BH-corrected; CI exclusion of zero alone does not establish BH significance.')
  }
  v3_table(ov_bind(manifest),'Figure_export_manifest_EN',out)
  writeLines(legends,file.path(out,'reports/Figure_legends_EN.md'),useBytes=TRUE)
}

abc_validate <- function(r,which,cfg,out) {
  checks<-list();add<-function(name,ok){checks[[length(checks)+1]]<<-data.frame(check=name,pass=isTRUE(ok));if(!isTRUE(ok))stop('Validation failed: ',name)}
  if(which=='A'){
    t<-r$results;add('Each genus/model/comparison unique',!anyDuplicated(paste(t$feature_id,t$model,t$comparison)))
    add('No invalid fit counted as discovery',!any(t$FDR_lt_005&!t$analysis_valid))
    add('Positive classified profile mass',all(r$coverage>0&r$coverage<=1+1e-10))
    add('Three independent-family community comparisons',nrow(r$community)==3)
    for(model in unique(t$model)){z<-t[t$model==model&t$eligible,];good<-z$analysis_valid&is.finite(z$P_raw);p<-rep(1,nrow(z));p[good]<-z$P_raw[good]
      add(paste(model,'joint genus BH exact'),max(abs(z$P_BH[good]-p.adjust(p,'BH')[good]))<1e-12)}
  }
  if(which=='B'){
    t<-r$results;add('Each SGB/model unique',!anyDuplicated(paste(t$feature_id,t$model)))
    add('Positive residual disease information',all(r$summary$within_family_disease_SS>0))
    add('Adequate residual df in estimated models',all(r$summary$residual_df[r$summary$status=='ESTIMATED']>10))
    for(model in unique(t$model)){z<-t[t$model==model,];add(paste(model,'BH exact'),max(abs(z$P_BH-p.adjust(z$P_raw,'BH')))<1e-12)}
    adult1<-readRDS(file.path(out,'private/Adult_B1_same_sample_analysis_design.rds'))$metadata
    adult2<-readRDS(file.path(out,'private/B2_analysis_design.rds'))$metadata
    add('Adult B1/B2 exact same IDs',identical(adult1$sample_id,adult2$sample_id))
    add('Adult models retain no minors',all(adult2$age>=18))
    for(model in r$summary$model){z<-readRDS(file.path(out,'private',paste0(model,'_analysis_design.rds')))
      add(paste(model,'each participant once'),!anyDuplicated(z$metadata$sample_id))
      add(paste(model,'all families disease-discordant'),all(tapply(z$metadata$fatty_liver,z$metadata$family_id,function(a)length(unique(a)))==2))}
    for(model in r$summary$model[r$summary$status=='ESTIMATED']){
      official<-readRDS(file.path(out,'private',paste0(model,'_official_LDM.rds')))
      z<-t[t$model==model,];p<-official$p.otu.omni;q<-official$q.otu.omni
      add(paste(model,'official P feature order exact'),identical(names(p),z$feature_id))
      add(paste(model,'official omnibus P exact'),max(abs(as.numeric(p)-z$P_raw))<1e-12)
      add(paste(model,'official native q exact'),max(abs(as.numeric(q)-z$q_official))<1e-12)
    }
  }
  if(which=='C'){
    t<-r$results;add('Exactly 100 frozen unique IDs',nrow(r$mapping)==100&&!anyDuplicated(r$mapping$SGB))
    add('Exactly 50 per external set',all(table(r$mapping$set)==50))
    add('Finite estimates and confidence intervals',all(is.finite(as.matrix(t[,c('beta','SE','CI_low','CI_high','P_raw')]))))
    for(model in unique(t$model)){z<-t[t$model==model,];add(paste(model,'four-readout BH exact'),nrow(z)==4&&max(abs(z$P_BH-p.adjust(z$P_raw,'BH')))<1e-12)}
    add('Cumulative disease-model direction preserved',all(t$outcome[grepl('abundance',t$readout)]=='Personal steatosis (0/1)'))
  }
  fl_write_csv(ov_bind(checks),file.path(out,'audit/Validation_checks.csv'))
}
