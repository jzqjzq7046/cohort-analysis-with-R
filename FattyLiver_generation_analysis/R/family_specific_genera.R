# Inert definitions; no participant data read or analysis runs while sourcing.
fg_protocol<-'family_genus_source_Fig1g_S5_rank_sum_20261029_v1'

fg_setup<-function(root){
  sys.source(file.path(root,'R/genus_generation_taxonomy.R'),envir=.GlobalEnv)
  cfg<-gg_setup(root)
  sys.source(file.path(root,'R/result_gallery.R'),envir=.GlobalEnv)
  for(p in c('readxl','matrixStats','ggrepel'))if(!requireNamespace(p,quietly=TRUE))stop('Missing project dependency: ',p)
  cfg
}

fg_log_bh<-function(lp){
  M<-length(lp);o<-order(lp);v<-lp[o]+log(M)-log(seq_len(M))
  q<-pmin(0,rev(cummin(rev(v))));a<-numeric(M);a[o]<-q;a
}

fg_source_audit<-function(path,out){
  raw<-suppressMessages(readxl::read_excel(path,sheet='Table_S5',col_names=FALSE,.name_repair='minimal'))
  keep<-vapply(raw[[4]],function(z)is.finite(suppressWarnings(as.numeric(z))),logical(1))
  t<-data.frame(genus=as.character(raw[[1]][keep]),family=as.character(raw[[2]][keep]),
    r=suppressWarnings(as.numeric(raw[[3]][keep])),P=suppressWarnings(as.numeric(raw[[4]][keep])))
  t<-t[is.finite(t$P)&is.finite(t$r),];M<-nrow(t)
  predicted<-2*pnorm(-abs(t$r)*sqrt(91));delta<-max(abs(predicted-t$P))
  g<-sort(unique(t$genus[t$r>0 & t$P<.05/M]))
  expected<-sort(paste0('g_',c('Anaerococcus','Corynebacterium','Oxalobacter','Pediococcus','Pseudomonas','Salmonella','Succinatimonas')))
  stopifnot(M==1824L,length(unique(t$genus))==76L,length(unique(t$family))==24L,
    identical(g,expected),delta<1e-10)
  rare<-suppressWarnings(as.numeric(raw[[8]]));valid<-is.finite(rare)
  a<-data.frame(item=c('Source tests','Source genera','Source families','Source N','Bonferroni raw-P cutoff',
    'Negative-log10 cutoff','Minimum reported prevalence count','r-to-P normal approximation maximum error',
    'Source seven enriched genera reproduced'),value=c(M,76,24,91,.05/M,-log10(.05/M),min(rare[valid]),delta,1))
  v3_table(a,'Original_article_S5_audit_EN',out)
  fl_write_csv(t,file.path(out,'audit/source_S5_tests.csv'))
  a
}

# Exact calibration for sparse pooled taxa: choose positive-carrier subsets and
# analytically count zero selections. No resampling or new asymptotic approximation.
fg_sparse_exact<-function(ranks2,observed,n1,max_carriers=20L){
  N<-length(ranks2);zero_rank<-attr(ranks2,'zero_rank');k<-attr(ranks2,'carriers')
  if(k>max_carriers)return(rep(NA_real_,length(n1)))
  pos<-ranks2[attr(ranks2,'positive')];n0<-N-k;ns<-sort(unique(n1));ans<-rep(NA_real_,length(n1))
  sums<-vector('list',min(k,max(ns))+1L);sums[[1]]<-c('0'=1)
  for(j in seq_len(length(sums)-1L)){
    v<-if(j==1L)pos else utils::combn(pos,j,FUN=sum,simplify=TRUE)
    sums[[j+1L]]<-table(v)
  }
  for(n in ns){
    score<-weight<-numeric()
    for(j in 0:min(k,n)){
      if(n-j>n0)next
      tab<-sums[[j+1L]];score<-c(score,as.numeric(names(tab))+(n-j)*zero_rank)
      weight<-c(weight,as.numeric(tab)*exp(lchoose(n0,n-j)-lchoose(N,n)))
    }
    d<-aggregate(weight,list(score=score),sum);d<-d[order(d$score),];prob<-d$x
    stopifnot(abs(sum(prob)-1)<1e-8);prob<-prob/sum(prob)
    j<-which(n1==n);at<-match(round(observed[j]),round(d$score))
    stopifnot(!anyNA(at));lower<-cumsum(prob);upper<-rev(cumsum(rev(prob)))
    ans[j]<-pmin(1,2*pmin(lower[at],upper[at]))
  }
  ans
}

fg_exact_floor<-function(ranks2,n1){
  N<-length(ranks2);sorted<-sort(ranks2);ns<-sort(unique(n1));v<-vapply(ns,function(n){
    low<-sorted[n];hi<-sorted[N-n+1L]
    pl<-lchoose(sum(sorted==low),n-sum(sorted<low))-lchoose(N,n)
    ph<-lchoose(sum(sorted==hi),n-sum(sorted>hi))-lchoose(N,n)
    min(0,log(2)+min(pl,ph))
  },numeric(1));v[match(n1,ns)]
}

# Exact zero-tie calibration when all nonzero observations are distinct.
# Conditional on j nonzero carriers, their rank sum has the ordinary Wilcoxon
# distribution on K positive observations. Mix over hypergeometric j.
fg_zero_tie_exact<-function(values,observed,n){
  N<-length(values);K<-sum(values>0);n0<-N-K
  if(anyDuplicated(values[values>0]))return(NA_real_)
  lower<-upper<-0
  for(j in max(0,n-n0):min(K,n)){
    w<-dhyper(j,K,n0,n);constant<-(n-j)*(n0+1)/2+j*n0+j*(j+1)/2
    q<-observed-constant
    if(j==0L||j==K){lo<-as.numeric(q>=0);hi<-as.numeric(q<=0)}else{
      lo<-pwilcox(floor(q),j,K-j);hi<-pwilcox(ceiling(q)-1,j,K-j,lower.tail=FALSE)
    }
    lower<-lower+w*lo;upper<-upper+w*hi
  }
  min(1,2*min(lower,upper))
}

fg_methods_exact_audit<-function(r,out){
  path<-file.path(out,'tables/Methods10_exact_calibration_EN.csv')
  if(file.exists(path))return(tf_read(path))
  a<-readRDS(file.path(out,'private/named_genus_input.rds'))
  fm<-readRDS(file.path(out,'private/family_label_map.rds'))
  t<-r$results[r$results$methods10_Bonferroni_support,]
  res<-t[,c('cohort','genus_id','genus','family_code','family_N','N','P_raw','P_Bonferroni_methods10')]
  res$exact_P<-NA_real_;res$positive_values_distinct<-FALSE;res$tested_methods_scope<-0L
  for(i in seq_len(nrow(t))){
    k<-if(t$cohort[i]=='All_ages')seq_len(nrow(a$metadata))else which(a$metadata$age>=4)
    v<-as.numeric(a$relative[t$genus_id[i],k]);raw<-fm$raw_family_id[match(t$family_code[i],fm$family_code)]
    within<-as.character(a$metadata$family_id[k])==raw
    stopifnot(sum(within)==t$family_N[i],length(v)==t$N[i])
    observed<-sum(rank(v,ties.method='average')[within])
    stopifnot(abs(observed-(t$W[i]+t$family_N[i]*(t$family_N[i]+1)/2))<1e-8)
    res$positive_values_distinct[i]<-!anyDuplicated(v[v>0])
    res$exact_P[i]<-fg_zero_tie_exact(v,observed,t$family_N[i])
    res$tested_methods_scope[i]<-sum(r$results$cohort==t$cohort[i]&r$results$prevalence_ge_10pct)
  }
  res$exact_P_Bonferroni_methods10<-pmin(1,res$exact_P*res$tested_methods_scope)
  res$exact_Bonferroni_support<-is.finite(res$exact_P_Bonferroni_methods10)&res$exact_P_Bonferroni_methods10<.05
  res$interpretation<-'Conditional individual-level calibration only; not covariate- or family-dependence-adjusted'
  v3_table(res,'Methods10_exact_calibration_EN',out);res
}

fg_fit<-function(A,m,label,out,family_map){
  stopifnot(identical(as.character(m$sample_id),colnames(A)),all(is.finite(A)),all(A>=0))
  ids<-rownames(A);G<-nrow(A);N<-ncol(A)
  fam<-factor(as.character(m$family_id),levels=family_map$raw_family_id);fam<-droplevels(fam)
  rawfam<-levels(fam);pub<-family_map$family_code[match(rawfam,family_map$raw_family_id)]
  n1<-as.integer(table(fam));F<-length(n1);M<-G*F
  R<-t(apply(A,1,rank,ties.method='average'));dimnames(R)<-dimnames(A)
  ranks<-t(rowsum(t(R),fam,reorder=FALSE));total<-t(rowsum(t(A),fam,reorder=FALSE))
  detected<-t(rowsum(t(1L*(A>0)),fam,reorder=FALSE))
  # rowsum(reorder=FALSE) uses first appearance, not factor order. Align explicitly.
  ranks<-ranks[,rawfam,drop=FALSE];total<-total[,rawfam,drop=FALSE];detected<-detected[,rawfam,drop=FALSE]
  tie<-vapply(seq_len(G),function(i){tt<-as.numeric(table(A[i,]));sum(tt^3-tt)},numeric(1))
  u<-sweep(ranks,2,n1*(n1+1)/2,'-');mu<-n1*(N-n1)/2
  variance<-outer((N+1-tie/(N*(N-1)))/12,n1*(N-n1))
  z<-sweep(u,2,mu,'-')/sqrt(variance);valid<-is.finite(z)&variance>0
  lp<-matrix(log(2)+pnorm(-abs(z),log.p=TRUE),nrow=G,ncol=F)
  lp[!valid]<-0;lp[lp>0]<-0
  family_mean<-sweep(total,2,n1,'/');rest_mean<-sweep(rowSums(A)-total,2,N-n1,'/')
  med<-matrix(NA_real_,G,F)
  for(j in seq_len(F))med[,j]<-matrixStats::rowMedians(A[,fam==rawfam[j],drop=FALSE])
  exact<-floor_lp<-matrix(NA_real_,G,F)
  carriers<-rowSums(A>0)
  for(i in seq_len(G)){
    r2<-2*R[i,];floor_lp[i,]<-fg_exact_floor(r2,n1)
    if(carriers[i]<=20L){
      attr(r2,'carriers')<-carriers[i];attr(r2,'positive')<-A[i,]>0
      attr(r2,'zero_rank')<-if(any(A[i,]==0))unique(r2[A[i,]==0]) else 0
      exact[i,]<-fg_sparse_exact(r2,2*ranks[i,],n1)
    }
  }
  # Real-data equivalence audit covers rare/common taxa, all observed family sizes.
  probe<-unique(c(seq_len(min(4,G)),order(carriers)[c(1,min(5,G),G)],which.max(rowMeans(A))))
  fj<-vapply(sort(unique(n1)),function(n)which(n1==n)[1],integer(1));audit<-list()
  for(i in probe)for(j in fj){
    w<-stats::wilcox.test(A[i,fam==rawfam[j]],A[i,fam!=rawfam[j]],alternative='two.sided',exact=FALSE,correct=FALSE)
    p<-exp(lp[i,j]);err<-abs(w$p.value-p)
    stopifnot(err<1e-10,abs(unname(w$statistic)-u[i,j])<1e-8)
    audit[[length(audit)+1L]]<-data.frame(cohort=label,genus_id=ids[i],family_code=pub[j],family_N=n1[j],
      reference_W=unname(w$statistic),fast_W=u[i,j],reference_P=w$p.value,fast_P=p,absolute_P_error=err)
  }
  fl_write_csv(ov_bind(audit),file.path(out,'audit',paste0(label,'_base_R_wilcoxon_equivalence.csv')))
  t<-data.frame(cohort=label,genus_id=rep(ids,F),genus=rep(sub('^g__','',ids),F),
    family_code=rep(pub,each=G),N=N,family_N=rep(n1,each=G),rest_N=rep(N-n1,each=G),
    pooled_detected_N=rep(carriers,F),pooled_prevalence=rep(carriers/N,F),
    family_detected_N=as.vector(detected),family_mean_RA_percent=100*as.vector(family_mean),
    rest_mean_RA_percent=100*as.vector(rest_mean),family_median_RA_percent=100*as.vector(med),
    W=as.vector(u),Z=as.vector(z),r_Z_over_sqrt_N=as.vector(z)/sqrt(N),
    rank_biserial=2*as.vector(u)/rep(n1*(N-n1),each=G)-1,
    P_raw=as.vector(exp(lp)),log_P_raw=as.vector(lp),negative_log10_P=-as.vector(lp)/log(10),
    P_underflow=as.vector(exp(lp))==0,estimable=as.vector(valid),
    exact_sparse_P=as.vector(exact),log_min_attainable_exact_P=as.vector(floor_lp),
    exact_audit_available=is.finite(as.vector(exact)))
  t$log_P_Bonferroni<-pmin(0,t$log_P_raw+log(M));t$P_Bonferroni<-exp(t$log_P_Bonferroni)
  t$log_P_BH<-fg_log_bh(t$log_P_raw);t$P_BH<-exp(t$log_P_BH)
  t$source_Bonferroni_support<-t$estimable & t$log_P_raw<log(.05/M)
  t$BH_support<-t$estimable&t$log_P_BH<log(.05)
  t$direction<-ifelse(t$Z>0,'Higher',ifelse(t$Z<0,'Lower','Equal'))
  t$source_threshold_P<-.05/M;t$source_threshold_negative_log10_P<--log10(.05/M);t$tests_in_family<-M
  t$exact_Bonferroni_support<-ifelse(t$exact_audit_available,t$exact_sparse_P<.05/M,NA)
  t$exact_threshold_unattainable<-t$log_min_attainable_exact_P>log(.05/M)
  t$source_only_sparse_failure<-t$source_Bonferroni_support&t$exact_audit_available&!t$exact_Bonferroni_support
  t$small_family_warning<-t$family_N<10L
  t$prevalence_ge_10pct<-t$pooled_prevalence>=.1
  j<-t$prevalence_ge_10pct;M10<-sum(j)
  t$P_Bonferroni_methods10<-t$P_BH_methods10<-NA_real_;t$methods10_Bonferroni_support<-FALSE
  t$P_Bonferroni_methods10[j]<-exp(pmin(0,t$log_P_raw[j]+log(M10)))
  t$P_BH_methods10[j]<-exp(fg_log_bh(t$log_P_raw[j]))
  t$methods10_Bonferroni_support[j]<-t$log_P_raw[j]<log(.05/M10)&t$estimable[j]
  stopifnot(nrow(t)==M,!anyNA(t$source_Bonferroni_support),all(t$P_raw>=0&t$P_raw<=1),
    max(abs(rowSums(total)-rowSums(A)))<1e-9)
  t
}

fg_build<-function(cfg,out){
  root<-cfg$project_root;p24<-lr_adopted(root,'24_Two_generation_reanalysis');p36<-lr_adopted(root,'36_Genus_generation_taxonomy')
  stage<-file.path(p24,'private/stage_prepare.rds');gen<-file.path(p36,'private/genus_relative_abundance.rds')
  prior<-file.path(p36,'private/analysis_result.rds')
  supp<-file.path(root,'literature_reviews/20261020_english_tables/sources/VallesColomer2022.xlsx')
  protected<-lr_hashes(c(stage,gen,prior,supp));old<-readRDS(prior);lr_check_hashes(old$protected)
  audit<-fg_source_audit(supp,out);x<-readRDS(stage);a<-readRDS(gen);m<-a$metadata
  stopifnot(identical(m$sample_id,x$metadata$sample_id),identical(m$family_id,x$metadata$family_id),
    identical(m$generation,x$metadata$generation),identical(colnames(a$relative),as.character(m$sample_id)),
    !anyDuplicated(m$sample_id),!anyNA(m$family_id),all(is.finite(m$age)),
    all(rowSums(table(m$family_id,m$generation)>0)>=2L))
  taxonomy<-old$taxonomy
  named<-taxonomy$feature_id[!taxonomy$unnamed_genome_genus_bin]
  A<-a$relative[match(named,rownames(a$relative)),,drop=FALSE]
  stopifnot(!anyNA(A),all(rowSums(A)>0),all(colSums(A)<=1+1e-10))
  fam<-sort(unique(as.character(m$family_id)));fm<-data.frame(raw_family_id=fam,family_code=sprintf('F%04d',seq_along(fam)))
  saveRDS(fm,file.path(out,'private/family_label_map.rds'))
  saveRDS(list(relative=A,metadata=m),file.path(out,'private/named_genus_input.rds'))
  family<-do.call(rbind,lapply(fam,function(f){d<-m[as.character(m$family_id)==f,]
    data.frame(family_code=fm$family_code[fm$raw_family_id==f],N=nrow(d),G0=sum(d$generation=='G0'),G1=sum(d$generation=='G1'),
      G2=sum(d$generation=='G2'),minors=sum(d$age<18),female=sum(d$sex=='Female'),age_min=min(d$age),age_max=max(d$age))
  }))
  v3_table(as.data.frame(table(family$N)), 'Observed_family_size_distribution_EN',out)
  v3_table(family,'Observed_family_composition_EN',out)
  filter<-taxonomy;filter$included_named_genus<-filter$feature_id%in%named
  filter$reason<-ifelse(filter$included_named_genus,'Uniquely aggregated named genus','Unnamed genome genus bin; excluded from source-style named-genus analysis')
  v3_table(filter,'Genus_filter_and_annotation_EN',out)
  results<-list();cohorts<-list(All_ages=seq_len(nrow(m)),Age_ge4=which(m$age>=4))
  for(label in names(cohorts)){
    j<-cohorts[[label]];B<-A[,j,drop=FALSE];B<-B[rowSums(B)>0,,drop=FALSE]
    v3_log('Figure1g family-genus ',label,': ',length(j),' people, ',nrow(B),' named genera, ',length(unique(m$family_id[j])),' labels')
    results[[label]]<-fg_fit(B,m[j,,drop=FALSE],label,out,fm)
  }
  t<-ov_bind(results)
  coverage<-data.frame(cohort=names(cohorts),participants=vapply(cohorts,length,integer(1)),
    families=vapply(cohorts,function(j)length(unique(m$family_id[j])),integer(1)),
    minors=vapply(cohorts,function(j)sum(m$age[j]<18),integer(1)),
    median_named_genus_mass_percent=vapply(cohorts,function(j)100*median(colSums(A[,j,drop=FALSE])),numeric(1)),
    removed_age_lt4=c(0,sum(m$age<4)))
  summary<-do.call(rbind,lapply(names(cohorts),function(label){d<-t[t$cohort==label,]
    data.frame(cohort=label,N=unique(d$N),families=length(unique(d$family_code)),genera=length(unique(d$genus_id)),
      tests=nrow(d),source_Bonferroni_cutoff=unique(d$source_threshold_P),source_negative_log10_cutoff=unique(d$source_threshold_negative_log10_P),
      source_higher_tests=sum(d$source_Bonferroni_support&d$Z>0),source_lower_tests=sum(d$source_Bonferroni_support&d$Z<0),
      source_higher_genera=length(unique(d$genus_id[d$source_Bonferroni_support&d$Z>0])),
      source_supported_families=length(unique(d$family_code[d$source_Bonferroni_support])),
      BH_supported_tests=sum(d$BH_support),sparse_exact_audited_tests=sum(d$exact_audit_available),
      source_significant_but_sparse_exact_fails=sum(d$source_only_sparse_failure),
      source_significant_exact_threshold_impossible=sum(d$source_Bonferroni_support&d$exact_threshold_unattainable),
      sparse_exact_Bonferroni_supported_tests=sum(d$exact_Bonferroni_support,na.rm=TRUE),
      methods10_genera=length(unique(d$genus_id[d$prevalence_ge_10pct])),methods10_Bonferroni_tests=sum(d$methods10_Bonferroni_support))
  }))
  gsum<-do.call(rbind,lapply(split(t,paste(t$cohort,t$genus_id,sep='|')),function(d){
    data.frame(cohort=d$cohort[1],genus_id=d$genus_id[1],genus=d$genus[1],pooled_prevalence=d$pooled_prevalence[1],
      higher_family_tests=sum(d$source_Bonferroni_support&d$Z>0),lower_family_tests=sum(d$source_Bonferroni_support&d$Z<0),
      min_log_P=min(d$log_P_raw),min_BH=min(d$P_BH),sparse_exact_available=all(d$exact_audit_available),
      source_only_sparse_failures=sum(d$source_only_sparse_failure))
  }))
  gs<-gsum[gsum$cohort=='All_ages'&gsum$higher_family_tests>0,];gs<-gs[order(gs$min_log_P,gs$genus_id),]
  selected<-head(gs,24L);selected$display_rank<-seq_len(nrow(selected))
  selected$selection_rule<-rep('At least one source-Bonferroni higher contrast; minimum raw P then genus ID; up to 24 display only',nrow(selected))
  means<-rowMeans(A);selected$mean_RA_percent<-100*means[match(selected$genus_id,rownames(A))]
  r<-list(protocol=fg_protocol,results=t,summary=summary,genus_summary=gsum,selected=selected,
    coverage=coverage,family_summary=family,source_audit=audit,protected=protected)
  fg_tables(r,out);lr_check_hashes(protected);r
}

fg_tables<-function(r,out){
  # The large full table is CSV only, avoiding an unusable million-cell HTML page.
  fl_write_csv(r$results,file.path(out,'tables/All_family_genus_Wilcoxon_tests_EN.csv'))
  v3_table(r$summary,'Family_genus_analysis_summary_EN',out)
  v3_table(r$coverage,'Analysis_population_and_coverage_EN',out)
  v3_table(r$genus_summary,'Genus_evidence_summary_EN',out)
  v3_table(r$selected,'Figure37_display_selection_EN',out)
  d<-r$results[r$results$source_Bonferroni_support,]
  fl_write_csv(d,file.path(out,'tables/Source_Bonferroni_supported_family_genus_tests_EN.csv'))
  fl_write_csv(r$results[r$results$exact_audit_available&r$results$source_Bonferroni_support,],
    file.path(out,'tables/Sparse_exact_calibration_of_source_signals_EN.csv'))
  v3_table(r$results[r$results$cohort=='All_ages'&r$results$genus_id%in%r$selected$genus_id&
      r$results$source_Bonferroni_support,], 'Figure37_supported_point_source_EN',out)
}

fg_plot<-function(r,cfg,out,stem='Figure37',title='Source-style family-versus-rest genus comparisons'){
  if(!nrow(r$selected))return(lr_empty('Family-associated genus patterns','No source-Bonferroni higher-family contrast'))
  s<-r$selected;d<-r$results[r$results$cohort=='All_ages'&r$results$genus_id%in%s$genus_id,]
  d$x<-match(d$genus_id,s$genus_id);d$mean_RA_percent<-s$mean_RA_percent[match(d$genus_id,s$genus_id)]
  d$state<-ifelse(d$source_Bonferroni_support,ifelse(d$Z>0,'Higher in family','Lower in family'),'Below threshold')
  # Deterministic offsets reveal repeated family comparisons; they do not change P.
  d$offset<-(as.integer(sub('^F','',d$family_code))%%23-11)/55;d$xplot<-d$x+d$offset
  top<-d[d$source_Bonferroni_support,];top<-top[order(top$log_P_raw,top$genus_id,top$family_code),]
  top<-top[!duplicated(top$genus_id),];top<-head(top,12L)
  top$label<-top$family_code
  mx<-max(d$negative_log10_P);cut<-unique(d$source_threshold_negative_log10_P)
  pal<-c('Higher in family'=unname(cfg$binary_colors['Comparison']),
    'Lower in family'=unname(cfg$binary_colors['Reference']),'Below threshold'='#BBBBBB')
  fill<-c('Higher in family'=unname(cfg$binary_colors['Comparison']),'Lower in family'='white','Below threshold'='white')
  p<-ggplot2::ggplot(d,ggplot2::aes(xplot,negative_log10_P))+
    ggplot2::geom_hline(yintercept=cut,linetype=2,linewidth=.3,colour='#555555')+
    ggplot2::geom_point(ggplot2::aes(size=mean_RA_percent,colour=state,fill=state),shape=21,stroke=.18,alpha=.58)+
    ggplot2::scale_colour_manual(values=pal,breaks=names(pal))+ggplot2::scale_fill_manual(values=fill,breaks=names(fill))+
    ggplot2::scale_size_area(name='Mean relative abundance (%)',max_size=2.5,breaks=scales::breaks_pretty(3))+
    ggplot2::scale_x_continuous(breaks=seq_len(nrow(s)),labels=s$genus,expand=ggplot2::expansion(add=.6))+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.02,.15)))+
    ggrepel::geom_text_repel(data=top,ggplot2::aes(label=label),seed=20261029,size=5.7/(72.27/25.4),
      colour='black',min.segment.length=0,segment.colour=NA,direction='y',nudge_y=mx*.035,
      box.padding=.4,point.padding=.15,max.overlaps=Inf,max.iter=20000)+
    ggplot2::labs(title=title,x=NULL,y='-log10(raw P)',colour=NULL,fill=NULL)+
    ggplot2::guides(size=ggplot2::guide_legend(order=2),colour=ggplot2::guide_legend(order=1),fill=ggplot2::guide_legend(order=1))+
    lr_theme()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=.5,size=5.7),
      legend.position='bottom',legend.box='vertical',legend.title=ggplot2::element_text(size=6),
      plot.margin=ggplot2::margin(6,8,5,5,'pt'))
  fl_write_csv(d,file.path(out,'tables',paste0(stem,'_all_displayed_point_source_EN.csv')))
  p
}

fg_methods10_view<-function(r){
  # A separately predeclared source-Methods scope, not an outcome-selected filter.
  a<-r;a$results<-r$results[r$results$prevalence_ge_10pct,]
  for(k in unique(a$results$cohort)){
    j<-a$results$cohort==k;M<-sum(j)
    a$results$source_threshold_P[j]<-.05/M
    a$results$source_threshold_negative_log10_P[j]<--log10(.05/M)
  }
  a$results$source_Bonferroni_support<-a$results$methods10_Bonferroni_support
  d<-a$results[a$results$cohort=='All_ages',]
  ids<-unique(d$genus_id[d$source_Bonferroni_support&d$Z>0])
  s<-data.frame(genus_id=ids,genus=sub('^g__','',ids))
  s$min_log_P<-vapply(ids,function(k)min(d$log_P_raw[d$genus_id==k]),numeric(1))
  s<-head(s[order(s$min_log_P,s$genus_id),],24L);s$display_rank<-seq_len(nrow(s))
  j<-match(s$genus_id,d$genus_id)
  s$mean_RA_percent<-(d$family_mean_RA_percent[j]*d$family_N[j]+d$rest_mean_RA_percent[j]*d$rest_N[j])/d$N[j]
  s$selection_rule<-rep('Prespecified >=10% prevalence scope; at least one Bonferroni higher contrast; minimum P then ID; display only',nrow(s))
  a$selected<-s;a
}

fg_figures<-function(r,cfg,out){
  p<-fg_plot(r,cfg,out)
  manifest<-lr_export(list(p),'Figure37_Family_associated_genera_source_style',out,cfg,width=183,height=150)
  a<-fg_methods10_view(r)
  v3_table(a$selected,'Figure37B_Methods10_display_selection_EN',out)
  d<-a$results[a$results$cohort=='All_ages'&a$results$source_Bonferroni_support,]
  v3_table(d,'Methods10_Bonferroni_supported_contrasts_EN',out)
  p2<-fg_plot(a,cfg,out,stem='Figure37B_Methods10',title='Family-versus-rest comparisons: prevalence >=10%')
  two<-lr_export(list(p2),'Figure37B_Family_genera_Methods10_scope',out,cfg,width=183,height=135)
  v3_table(ov_bind(list(manifest,two)),'Figure_export_manifest_EN',out)
}

fg_reports<-function(r,out){
  calibrated<-fg_methods_exact_audit(r,out)
  writeLines(c('# Figure 37 | Family-associated genus abundance patterns',
    sprintf('The source-style display includes %d selected genera from %d participants in %d primary-family labels. Every point represents one family-versus-rest comparison for a displayed genus. The complete testing family, including unshown genera, is retained in the result table.',nrow(r$selected),r$summary$N[1],r$summary$families[1]),
    'The vertical axis is the negative base-10 logarithm of the two-sided, tie-corrected asymptotic Wilcoxon rank-sum P value, without continuity correction. Bubble area is proportional to the cohort-average relative abundance of that genus. Filled orange and open blue circles denote higher and lower family ranks passing the global Bonferroni threshold; grey circles do not pass. The dashed line uses 0.05 divided by all eligible named-genus-by-family tests. Small deterministic horizontal offsets expose overlapping comparisons.',
    'The main display contains up to 24 genera with at least one higher-family Bonferroni contrast, ranked by minimum raw P then stable genus ID. Labels identify up to 12 genus peaks and use anonymous family codes. Selection affects display only. BH and the separately labelled >=10% prevalence scope appear in tables. The age>=4 source-range sensitivity is reported separately; the all-age main analysis retains minors.',
    'This applies the Figure 1g / Supplementary Table S5 framework of Valles-Colomer et al. Nature Microbiology 2022, doi:10.1038/s41564-021-01021-8. Relative abundance replaces source quantitative cell abundance, and a restrained direction palette replaces source family colours. Primary labels do not verify kinship or cohabitation. The test does not adjust age, sex or other covariates or dependence among the remaining families. Small observed family sizes and sparse taxa can make asymptotic P values severely anti-conservative. Exact sparse-taxon calibration and attainable-P flags must be inspected before biological interpretation. These results are exploratory, not confirmed heritability or microbial transmission.',
    '', '# Companion Figure 37B | Source-Methods prevalence scope',
    'The same source-style display and already calculated Wilcoxon statistics are shown for the separately predeclared genera detected in >=10% of participants. Bonferroni is recomputed across all genus-by-family tests in this scope, not only the plotted genera or significant tests. This follows the generic source-Methods prevalence rule, whereas Figure 37 follows the rare-genus-inclusive S5 testing scope. Colours, area encoding, anonymous family labels and main cohort are unchanged. Neither scope resolves small-family approximation or within-family dependence. All inference is reused, without a new test or an adjusted model.'),
    file.path(out,'reports/Figure_legend_EN.md'),useBytes=TRUE)
  writeLines(c('# Methods and source fidelity',
    'Valles-Colomer M, Bacigalupe R, Vieira-Silva S, et al. Variation and transmission of the human gut microbiota across multiple familial generations. Nature Microbiology. 2022;7:87–96. doi:10.1038/s41564-021-01021-8. https://www.nature.com/articles/s41564-021-01021-8',
    'Source locations: Figure 1g, Supplementary Table S5, Microbiome and metadata associations. The article does not release a complete analysis script or document continuity/exact settings. Our explicit asymptotic/no-continuity implementation is checked against base R. Z/sqrt(N) reproduces the published S5 r-to-P relation. The source S5 has 1,824 tests, including low-prevalence genera. Its 4.56 figure threshold matches 0.05/1824 and recovers all seven named enriched genera. This supports the reconstructed Bonferroni rule but does not establish an unreported original script.',
    'Data: completed expanded baseline and conservative genus aggregation. All nonzero named genera are included in the S5-style testing family; GGB groups, ambiguous assignments and unclassified mass are not renamed or redistributed. Overall prevalence >=10% is a separately corrected, predeclared subset because the generic source Methods contradict the low-prevalence S5 entries. Original source quantitative abundance, female-only verified pedigrees and >=4-year scope are not silently equated with our relative abundance, both-sex primary labels and all-age main cohort.',
    'Unit: each person contributes one abundance per genus, each family is compared with all remaining people. The same people recur in different family-versus-rest tests and members of the remainder can be related. Wilcoxon is a source-style descriptive test, not a family-cluster-robust model. Bonferroni covers every named genus x family comparison in each cohort; BH is also shown. Invalid constant tests count conservatively as P=1. No disease association selection, top-30 filter, covariate selection or re-normalization is used.',
    'Exact calibration: for each genus with <=20 detected people, enumerate the subsets of those carriers at each possible number selected into a family and sum rank scores. Weight each subset by choose(N_zero,n_family-n_carriers_selected)/choose(N,n_family); include all zero-only selections. This gives the exact conditional sampling distribution of the observed Wilcoxon rank sum, with a two-sided probability 2*min(lower-tail,upper-tail), capped at 1. No random sampling or outcome-dependent seed is involved. The exact diagnostic retains the full original Bonferroni denominator; missing audits on other taxa do not mean negative results. Minimum attainable exact P is computed from the combinatorial probability of extreme tied rank selections. These audits do not repair family dependence or demographic imbalance.',
    'Additional calibration of all source-Methods >=10%-prevalence Bonferroni signals: when positive abundances are distinct and only zeros are tied, condition on the hypergeometric number j of positive carriers in a family. The selected positive positions have a Wilcoxon rank-sum distribution among K positive observations. Sum its inclusive lower and upper tails with the hypergeometric weights, then take 2*min(tails), capped at 1. This exact diagnostic uses the full original >=10% testing denominator. It does not rerun or replace the source-style tests and is not a covariate- or dependence-adjusted analysis. Positive-value ties would leave the audit unavailable, rather than force an incorrect calculation.',
    'No absolute microbial load, adjusted disease effects, independent-family replication, strain sharing, inheritance or causal interpretation is claimed. Public family codes conceal raw primary labels; all original inputs and completed modules remain unchanged.'),
    file.path(out,'reports/Methods_and_references_EN.md'),useBytes=TRUE)
  lines<-c('# 37号：家系特异性菌属探索',
    '按原文Figure 1g / 补充表S5逐一家系与其余人群进行双侧Wilcoxon秩和检验。没有加入协变量；不是26/36号跨代丰度回归。',
    sprintf('全龄：%d人、%d个主家系、%d个已命名菌属，合计%d次比较。源表式渐近检验经全范围Bonferroni后有%d个较高比较，涉及%d个属。',r$summary$N[1],r$summary$families[1],r$summary$genera[1],r$summary$tests[1],r$summary$source_higher_tests[1],r$summary$source_higher_genera[1]),
    sprintf('本次源表式原始P门槛为%.6g，-log10(P)门槛为%.4f；没有照搬原图4.56。',r$summary$source_Bonferroni_cutoff[1],r$summary$source_negative_log10_cutoff[1]),
    sprintf('全龄稀有属精确校准覆盖%d项比较；其中%d项源式显著比较未通过精确概率的相同Bonferroni门槛，稀有属精确支持为%d项。',r$summary$sparse_exact_audited_tests[1],r$summary$source_significant_but_sparse_exact_fails[1],r$summary$sparse_exact_Bonferroni_supported_tests[1]),
    sprintf('正文>=10%%检出范围含%d属；全龄源式Bonferroni支持%d项。具体菌属与家系见Methods10_Bonferroni_supported_contrasts_EN。',r$summary$methods10_genera[1],r$summary$methods10_Bonferroni_tests[1]),
    '本队列各家系观察人数很少。小样本与大量并列零值会使渐近P严重偏小；不能将源式Bonferroni支持项直接当作可靠家系富集发现。请同时查看稀有属精确校准和最小可达P标记。',
    '主分析保留全部年龄；另列>=4岁原文范围敏感性。未要求三代均有成员，不分析MASLD。',
    '补充表S5实际包含低检出属；因此源表式全非零已命名属与正文>=10%检出范围分别报告并分别校正，不能混成同一显著名单。',
    sprintf('>=10%%范围的源式支持项另经精确条件秩和概率核查：全龄%d项中%d项仍过原完整检验范围的Bonferroni门槛；>=4岁%d项中%d项通过。精确核查也没有控制协变量或家系依赖。',
      sum(calibrated$cohort=='All_ages'),sum(calibrated$cohort=='All_ages'&calibrated$exact_Bonferroni_support),
      sum(calibrated$cohort=='Age_ge4'),sum(calibrated$cohort=='Age_ge4'&calibrated$exact_Bonferroni_support)),
    '图形为菌属横轴、-log10(P)纵轴，气泡面积按菌属全队列平均相对丰度，填充/空心表示较高/较低。最多24属仅用于主图；完整表保留所有检验。',
    '原文使用定量菌量；本次使用既有相对丰度及主家系标签，是已说明的数据适配。标签未核实直接亲缘或同住，不解释为传播或遗传。')
  supported<-calibrated[calibrated$cohort=='All_ages'&calibrated$exact_Bonferroni_support,]
  if(nrow(supported))lines<-c(lines,'', '## 通过追加精确校准的全龄探索信号',
    vapply(seq_len(nrow(supported)),function(i)sprintf('%s，匿名家系%s（观察%d人）：精确原始P=%.6g，完整>=10%%检验范围Bonferroni校正P=%.6g。',
      supported$genus[i],supported$family_code[i],supported$family_N[i],supported$exact_P[i],supported$exact_P_Bonferroni_methods10[i]),character(1)),
    '精确概率只修正小样本秩和分布，不控制年龄、性别、技术或同一家系相关性；因此仍是探索信号，不能视为独立家系验证或机制证据。')
  writeLines(lines,file.path(out,'reports/Results_summary_zh.md'),useBytes=TRUE)
}

fg_main<-function(args=character(),root=Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')){
  if(getRversion()<'4.5'){
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts/37_family_specific_genera.R')),args))
    if(status!=0L)stop('Project-R family-specific genus entry failed');return(invisible(NULL))
  }
  if(any(!args%in%c('--run','--check','--help','--figures-only','--resume-figures')))stop('Unknown mode')
  if('--help'%in%args){cat('37_family_specific_genera.R [--run | --check | --figures-only | --resume-figures]\n');return(invisible(NULL))}
  cfg<-fg_setup(root);cfg$output_module<-unname(cfg$output_modules['family_specific_genera'])
  cfg$run_id<-unname(cfg$run_ids['family_specific_genera']);cfg$run_id_explicit<-FALSE
  if('--check'%in%args){invisible(parse(file.path(root,'R/family_specific_genera.R')));cat('PASS: syntax/dependencies; no participant data read\n');return(invisible(NULL))}
  if('--figures-only'%in%args){
    out<-lr_adopted(root,cfg$output_module);r<-readRDS(file.path(out,'private/analysis_result.rds'));lr_check_hashes(r$protected)
    stamp<-format(Sys.time(),'%Y%m%d_%H%M%S');dest<-file.path(out,'audit',paste0('prior_figures_',stamp));dir.create(dest)
    file.copy(list.files(file.path(out,'figures'),full.names=TRUE),dest)
    fitfile<-file.path(out,'tables/All_family_genus_Wilcoxon_tests_EN.csv');fit_hash<-unname(tools::md5sum(fitfile))
    fg_figures(r,cfg,out);fg_reports(r,out);lr_check_hashes(r$protected)
    stopifnot(identical(fit_hash,unname(tools::md5sum(fitfile))))
    fl_publish_organized_results(cfg);return(invisible(out))
  }
  resume<-NULL
  if('--resume-figures'%in%args){
    dirs<-list.dirs(file.path(root,'outputs',cfg$output_module,'history'),recursive=FALSE,full.names=TRUE)
    dirs<-dirs[file.exists(file.path(dirs,'private/analysis_result.rds'))]
    if(!length(dirs))stop('No saved family-genus result to recover')
    resume<-dirs[which.max(file.info(file.path(dirs,'private/analysis_result.rds'))$mtime)]
  }
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(p in c('tables','figures','reports','audit','private'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
  writeLines(c('RUNNING',fg_protocol),file.path(out,'run_status.txt'))
  tryCatch({
    if(is.null(resume))r<-fg_build(cfg,out)else{
      r<-readRDS(file.path(resume,'private/analysis_result.rds'));lr_check_hashes(r$protected)
      for(p in c('tables','audit','private')){
        dst<-if(p=='audit')file.path(out,'audit/saved_analysis_audit')else file.path(out,p);dir.create(dst,recursive=TRUE,showWarnings=FALSE)
        files<-list.files(file.path(resume,p),full.names=TRUE);if(length(files))stopifnot(all(file.copy(files,dst,recursive=TRUE)))
      }
    }
    saveRDS(r,file.path(out,'private/analysis_result.rds'));fg_reports(r,out)
    fun<-c('fg_plot','fg_figures','lr_export','tf_export','nf_theme')
    src<-file.path(out,'audit/actual_plot_source.R')
    writeLines(unlist(lapply(fun,function(n)c(paste0(n,' <- '),deparse(get(n),width.cutoff=120L),''))),src)
    nf_qa_command(cfg,'validate_figure.py',c(shQuote(src),'--json'),file.path(out,'audit/source_preflight.json'))
    fg_figures(r,cfg,out);lr_check_hashes(r$protected)
    fl_write_csv(r$protected,file.path(out,'audit/protected_source_files.csv'))
    writeLines(capture.output(sessionInfo()),file.path(out,'audit/sessionInfo.txt'))
    writeLines(c('COMPLETED',fg_protocol),file.path(out,'run_status.txt'))
    target<-fl_output_complete_run(ctx);writeLines(basename(target),file.path(root,'outputs',cfg$output_module,'LATEST_RUN.txt'))
    v3_log('COMPLETED: ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
