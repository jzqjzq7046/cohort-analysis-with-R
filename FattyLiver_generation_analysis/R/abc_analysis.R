# Published methods reused with explicitly documented cohort adaptations.
abc_protocol <- 'ABC_Olivera_Zhu_Asnicar_20261027_v1'

abc_setup <- function(root) {
  if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
  .libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),.libPaths()))
  for(f in c('config.R','R/common.R','R/entry_helpers.R','R/confounders.R','R/nature_export_approved.R',
    'R/reanalysis_numerical_helpers.R','R/reanalysis_20261024.R','R/taxa_followup.R',
    'R/taxa_followup_figures.R','R/family_steatosis_microbiome.R','R/abc_figures.R','R/abc_ldm_acceleration.R'))
    sys.source(file.path(root,f),envir=.GlobalEnv)
  Sys.setenv(PYTHONPATH=paste(c(file.path(root,'.figure_qa_python'),Sys.getenv('PYTHONPATH')),collapse=.Platform$path.sep),
    OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
  stopifnot(all(vapply(c('Maaslin2','lmerTest','geepack','LDM','GUniFrac','readxl','digest','vegan'),requireNamespace,logical(1),quietly=TRUE)))
  fl_config(root)
}

abc_source <- function(root) {
  b<-file.path(root,'outputs/27_Family_steatosis_microbiome');id<-readLines(file.path(b,'LATEST_RUN.txt'),1)
  p<-file.path(b,'current',id);stopifnot(readLines(file.path(p,'run_status.txt'),1)=='COMPLETED')
  f<-file.path(p,'private/analysis_population.rds');x<-readRDS(f)
  stopifnot(identical(x$protocol,fp_protocol),!anyDuplicated(x$metadata$sample_id),
    identical(colnames(x$relative),as.character(x$metadata$sample_id)),max(abs(colSums(x$relative)-1))<1e-10)
  list(x=x,file=f,md5=unname(tools::md5sum(f)))
}

abc_design <- function(m,fields,set=NULL,age_spline=TRUE) {
  nz<-v3_nuisance(m,fields,age_spline=age_spline);Z<-nz$Z
  original<-colnames(Z)
  if(!is.null(set)){S<-model.matrix(~factor(set)-1);colnames(S)<-paste0('set_',seq_len(ncol(S)));Z<-cbind(S,Z);original<-c(colnames(S),original)}
  q<-qr(Z,tol=1e-9);keep<-sort(q$pivot[seq_len(q$rank)]);dropped<-setdiff(seq_len(ncol(Z)),keep);Z<-Z[,keep,drop=FALSE]
  # Orthogonal nuisance basis, same column space; removes aliases already absorbed by set indicators.
  Q<-qr.Q(qr(Z));stopifnot(ncol(Q)==ncol(Z),max(abs(crossprod(Q)-diag(ncol(Q))))<1e-8)
  list(Z=Q,mapping=data.frame(original=original,retained=seq_along(original)%in%keep),
    dropped=c(nz$dropped,original[dropped]),metadata=nz$metadata)
}

abc_gee <- function(y,e,Z,m) {
  d<-data.frame(y=y,exposure=e,Z,row.names=m$sample_id,check.names=FALSE)
  colnames(d)<-c('y','exposure',sprintf('z%03d',seq_len(ncol(Z))))
  f<-reformulate(colnames(d)[-1],response='y',intercept=FALSE)
  fit<-geepack::geeglm(f,data=d,id=as.integer(factor(m$family_id)),family=gaussian,corstr='independence',std.err='san.se')
  stopifnot(fit$geese$error==0);b<-unname(coef(fit)['exposure']);se<-sqrt(vcov(fit)['exposure','exposure'])
  stopifnot(all(is.finite(c(b,se))),se>0)
  list(row=data.frame(beta=b,SE=se,CI_low=b-qnorm(.975)*se,CI_high=b+qnorm(.975)*se,P_raw=2*pnorm(-abs(b/se))),fit=fit)
}

abc_genus <- function(x,out) {
  ft<-x$features[match(rownames(x$relative),x$features$feature_id),];stopifnot(!anyNA(ft$feature_id))
  gs<-regmatches(ft$taxonomy,gregexpr('g__[^|,; ]+',ft$taxonomy))
  gs<-lapply(gs,unique);ok<-vapply(gs,function(v)length(v)==1&&!grepl('unclassified|unknown|unidentified',v,ignore.case=TRUE),logical(1))
  genus<-vapply(gs,function(v)if(length(v)==1)v else NA_character_,character(1))
  map<-data.frame(feature_id=ft$feature_id,official_id=ft$official_id,genus=genus,unique_genus_annotations=lengths(gs),
    included=ok,genus_annotations=vapply(gs,paste,collapse=';',character(1)))
  v3_table(map,'00_SGB_to_genus_mapping_EN',out)
  G<-rowsum(x$relative[ok,,drop=FALSE],genus[ok],reorder=TRUE);fraction<-colSums(G)
  stopifnot(all(fraction>0),all(fraction<=1+1e-10),max(abs(fraction-colSums(x$relative[ok,,drop=FALSE])))<1e-10)
  saveRDS(data.frame(sample_id=colnames(G),classified_fraction=fraction),file.path(out,'private/genus_sample_coverage.rds'))
  v3_table(data.frame(SGBs_total=nrow(x$relative),SGBs_assigned=sum(ok),SGBs_excluded=sum(!ok),genus_bins=nrow(G),
    classified_fraction_median=median(fraction),classified_fraction_Q1=quantile(fraction,.25),
    classified_fraction_Q3=quantile(fraction,.75),classified_fraction_min=min(fraction)), '01_Taxonomy_coverage_EN',out)
  z<-x;z$relative<-G;z$features<-data.frame(feature_id=rownames(G),taxonomy=paste0(rownames(G),'|s__',sub('g__','',rownames(G))))
  list(x=z,G=G,fraction=fraction,mapping=map)
}

abc_A <- function(x,cfg,out,cache) {
  g<-abc_genus(x,out);sets<-fp_comparisons(x,out);effects<-alpha<-community<-ordinations<-list()
  for(key in names(sets$sets)){
    for(model in c('B1','Adult_B1_same_sample','B2')){
      m<-sets$sets[[key]][[model]];effects[[paste(key,model)]]<-fp_fit(g$x,m,model,key,cfg,out,file.path(cache,'genus'))
    }
    m<-sets$sets[[key]]$B1;m<-droplevels(m[order(m$family_id,m$sample_id),])
    A<-t(g$G[,match(m$sample_id,colnames(g$G)),drop=FALSE]);A<-A/rowSums(A)
    H<-vegan::diversity(A,index='shannon');Z<-abc_design(m,c(v3_P2,'n_older_observed'))$Z
    a<-abc_gee(H,m$older_any_case,Z,m);alpha[[key]]<-cbind(comparison=key,model='B1_family_GEE',N=nrow(m),families=length(unique(m$family_id)),a$row)
    set.seed(20261028L+match(key,names(sets$sets)));ix<-sort(vapply(split(seq_len(nrow(m)),m$family_id),function(ii)ii[sample.int(length(ii),1)],integer(1)))
    sm<-droplevels(m[ix,]);stopifnot(!anyDuplicated(sm$family_id));saveRDS(sm,file.path(out,'private',paste0(key,'_one_per_family.rds')))
    D<-vegan::vegdist(A[ix,,drop=FALSE],'bray');nz<-abc_design(sm,c(v3_P2,'n_older_observed'));Zs<-nz$Z
    e<-sm$older_any_case
    set.seed(20261128L+match(key,names(sets$sets)));perm<-vegan::adonis2(D~e,permutations=1000)
    bd<-vegan::betadisper(D,factor(e),type='median');set.seed(20261228L+match(key,names(sets$sets)))
    disp<-vegan::permutest(bd,permutations=1000)
    distance_metadata<-data.frame(e=e,row.names=sm$sample_id);distance_metadata$Zs<-I(Zs)
    dm<-GUniFrac::dmanova(D~Zs+e-1,data=distance_metadata,returnG=FALSE)
    saveRDS(list(adonis=perm,dispersion=disp,dmanova=dm),file.path(out,'audit',paste0(key,'_community_fits.rds')))
    dt<-dm$aov.tab;if(is.null(dt))dt<-dm$aov.table
    if(is.null(dt)){capture.output(str(dm),file=file.path(out,'audit/dmanova_structure.txt'));stop('Inspect dmanova output before extracting P')}
    di<-1L;stopifnot(rownames(dt)[di]=='Model(Adjusted)')
    pn<-grep('Pr|p.?value',colnames(dt),ignore.case=TRUE,value=TRUE)[1]
    rn<-grep('R2|R.sq',colnames(dt),ignore.case=TRUE,value=TRUE)[1]
    community[[key]]<-data.frame(comparison=key,N=nrow(sm),reference_N=sum(e==0),exposed_N=sum(e==1),
      Wilcoxon_Shannon_P=wilcox.test(H[ix]~e,exact=FALSE)$p.value,
      adonis_R2=perm$R2[1],adonis_P=perm$`Pr(>F)`[1],dispersion_P=disp$tab$`Pr(>F)`[1],
      dmanova_P=as.numeric(dt[di,pn]),dmanova_R2=if(!is.na(rn))as.numeric(dt[di,rn])else NA_real_,
      genus_classified_fraction_median=median(g$fraction[match(sm$sample_id,names(g$fraction))]))
    pc<-cmdscale(D,k=2,eig=TRUE,add=TRUE);ordinations[[key]]<-data.frame(comparison=key,sample_id=sm$sample_id,
      exposure=factor(e,0:1,c('No older case','Older case')),PC1=pc$points[,1],PC2=pc$points[,2])
    saveRDS(data.frame(m,shannon=H),file.path(out,'private',paste0(key,'_genus_shannon.rds')))
  }
  t<-fp_adjust(ov_bind(effects));t$taxonomic_level<-'Unambiguous genus / genome genus bin';t$genus<-sub('^g__','',t$feature_id)
  t$effect_unit<-'Log2 zero-replaced full-profile genus relative abundance; older case minus no observed older case'
  t$FDR_family<-paste0(t$model,'_all_eligible_genera_x_three_comparisons')
  a<-ov_bind(alpha);a$P_BH<-p.adjust(a$P_raw,'BH');c<-ov_bind(community)
  for(p in c('Wilcoxon_Shannon_P','adonis_P','dispersion_P','dmanova_P'))c[[paste0(p,'_BH')]]<-p.adjust(c[[p]],'BH')
  v3_table(t,'A_All_genus_associations_EN',out);v3_table(t[t$FDR_lt_005,],'A_BH_supported_genera_EN',out)
  v3_table(a,'A_Adjusted_genus_Shannon_EN',out);v3_table(c,'A_Genus_community_tests_EN',out)
  sm<-ov_bind(lapply(split(t,list(t$comparison,t$model),drop=TRUE),function(z)data.frame(comparison=z$comparison[1],model=z$model[1],
    N=z$N[1],tested=sum(z$eligible),valid=sum(z$analysis_valid),BH_lt_005=sum(z$FDR_lt_005))))
  v3_table(sm,'A_Genus_model_summary_EN',out)
  list(results=t,alpha=a,community=c,summary=sm,ordinations=ordinations,coverage=g$fraction)
}

abc_discordant <- function(m,fields,same_generation=NULL,adult=FALSE) {
  keep<-v3_cc(m,c('family_id','fatty_liver',fields));if(adult)keep<-keep&m$age>=18
  if(!is.null(same_generation))keep<-keep&as.character(m$generation)==same_generation
  m<-droplevels(m[keep,,drop=FALSE]);disc<-tapply(m$fatty_liver,m$family_id,function(v)length(unique(v))==2)
  m<-droplevels(m[as.character(m$family_id)%in%names(disc)[disc],,drop=FALSE]);m<-m[order(m$family_id,m$sample_id),]
  stopifnot(!anyDuplicated(m$sample_id),all(tapply(m$fatty_liver,m$family_id,function(v)length(unique(v)))==2));m
}

abc_B_fit <- function(x,m,fields,key,out,cache) {
  stopifnot(nrow(m)>20);nz<-abc_design(m,fields,as.character(m$family_id));Z<-nz$Z;e<-m$fatty_liver
  er<-drop(e-Z%*%crossprod(Z,e));df<-nrow(m)-ncol(Z)-1
  if(sum(er^2)<=1e-6||df<=10){
    saveRDS(list(metadata=m,Z=Z,disease_residual=er,design_mapping=nz$mapping),file.path(out,'private',paste0(key,'_analysis_design.rds')))
    return(list(results=NULL,summary=data.frame(model=key,N=nrow(m),cases=sum(e),families=nlevels(m$family_id),tested=0,
      global_P_frequency=NA_real_,global_P_AST=NA_real_,global_P_omnibus=NA_real_,BH_lt_005=NA_integer_,official_q_lt_005=NA_integer_,
      scale_disagreements=NA_integer_,residual_df=df,within_family_disease_SS=sum(er^2),status='NOT_RELIABLY_ESTIMABLE: residual df <=10 or exposure not identifiable')))
  }
  A<-t(x$relative[,match(m$sample_id,colnames(x$relative)),drop=FALSE]);nzobs<-colSums(A>0);pv<-colMeans(A>0)
  keep<-pv>=.1&nzobs>=5;A<-A[,keep,drop=FALSE];stopifnot(ncol(A)>0,all(rowSums(A)>0))
  signature<-digest::digest(list(abc_protocol,key,A,Z,e,LDM=as.character(packageVersion('LDM')),
    inference=c(max_permutations=10000,rejection_stop=100,FDR=.05,scale_rows=0,cores=1,omni3=0),
    acceleration=body(abc_ldm_acceleration)),algo='sha256')
  v3_table(data.frame(model=key,feature_id=rownames(x$relative),prevalence=pv,nonzero=nzobs,eligible=keep),paste0('B_Filter_',key,'_EN'),out)
  saveRDS(list(metadata=m,Z=Z,disease_residual=er,design_mapping=nz$mapping),file.path(out,'private',paste0(key,'_analysis_design.rds')))
  fl_write_csv(nz$mapping,file.path(out,'audit',paste0(key,'_nuisance_mapping.csv')))
  tf_log('B ',key,': ',nrow(m),' participants; ',nlevels(m$family_id),' discordant families; ',ncol(A),' SGBs; official LDM')
  cp<-file.path(cache,paste0(key,'.rds'))
  fit_cache<-tf_checkpoint(cp,signature,function(){
    Zc<-scale(Z,center=TRUE,scale=FALSE);qz<-qr(Zc,tol=1e-9);Zldm<-Zc[,sort(qz$pivot[seq_len(qz$rank)]),drop=FALSE]
    # Official package, set indicators AND restricted within-family permutations.
    dat<-data.frame(e=e,family_id=as.character(m$family_id),row.names=m$sample_id)
    dat$A<-I(A);dat$Zldm<-I(Zldm)
    accelerated<-abc_ldm_acceleration(m$family_id,function()LDM::ldm(A | Zldm ~ e,data=dat,cluster.id=family_id,
      perm.within.type='free',perm.between.type='none',n.perm.max=10000,n.rej.stop=100,
      seed=20261029L+match(key,c('B1','G0_same_generation','G2_same_generation','Adult_B1_same_sample','B2')),
      fdr.nominal=.05,scale.otu.table=FALSE,freq.scale.only=FALSE,test.omni3=FALSE,
      n.cores=1,verbose=TRUE))
    accelerated
  })
  ans<-fit_cache$result
  fl_write_csv(fit_cache$equivalence,file.path(out,'audit',paste0(key,'_official_statistic_equivalence.csv')))
  fl_write_csv(data.frame(LDM_version=as.character(packageVersion('LDM')),statistic_calls=fit_cache$calls,
    statistic_hash=fit_cache$official_statistic_hash,ldm_function_hash=fit_cache$official_ldm_hash),file.path(out,'audit',paste0(key,'_official_function_hashes.csv')))
  saveRDS(ans,file.path(out,'private',paste0(key,'_official_LDM.rds')))
  getp<-function(nm){a<-ans[[nm]];if(is.null(a))stop('Missing official LDM output: ',nm);v<-if(is.matrix(a))a[,ncol(a)]else a;stopifnot(length(v)==ncol(A));as.numeric(v)}
  # Descriptive adjusted coefficients; LDM omnibus P itself has no unique scale/direction.
  bf<-drop(crossprod(er,A))/sum(er^2);bt<-drop(crossprod(er,asin(sqrt(A))))/sum(er^2)
  t<-data.frame(model=key,feature_id=colnames(A),beta_relative=bf,beta_percentage_points=100*bf,beta_AST=bt,
    P_frequency=getp('p.otu.freq'),P_AST=getp('p.otu.tran'),P_raw=getp('p.otu.omni'),
    q_official=getp('q.otu.omni'),prevalence=pv[keep],N=nrow(m),cases=sum(e),families=nlevels(m$family_id),
    residual_df=df,within_family_disease_SS=sum(er^2),minors=sum(m$age<18),stringsAsFactors=FALSE)
  stopifnot(all(t$P_raw>0&t$P_raw<=1));t$P_BH<-p.adjust(t$P_raw,'BH');t$BH_lt_005<-t$P_BH<.05
  t$scale_direction_agrees<-sign(t$beta_relative)==sign(t$beta_AST)
  t$direction<-ifelse(!t$scale_direction_agrees,'Scale-dependent',ifelse(t$beta_relative<0,'Lower with steatosis','Higher with steatosis'))
  a<-match(t$feature_id,x$features$feature_id);t$SGB_id<-x$features$official_id[a];t$species_label<-ov_species_label(x$features$taxonomy[a])
  globals<-data.frame(model=key,N=nrow(m),cases=sum(e),families=nlevels(m$family_id),tested=ncol(A),
    global_P_frequency=as.numeric(ans$p.global.freq)[1],global_P_AST=as.numeric(ans$p.global.tran)[1],
    global_P_omnibus=as.numeric(ans$p.global.omni)[1],BH_lt_005=sum(t$BH_lt_005),official_q_lt_005=sum(t$q_official<.05),
    scale_disagreements=sum(!t$scale_direction_agrees),residual_df=df,within_family_disease_SS=sum(er^2),status='ESTIMATED')
  list(results=t,summary=globals)
}

abc_B <- function(x,cfg,out,cache) {
  m<-x$metadata;adult<-abc_discordant(m,v3_B2,adult=TRUE)
  specs<-list(B1=list(abc_discordant(m,v3_B1),v3_B1),
    G0_same_generation=list(abc_discordant(m,v3_B1,'G0'),v3_B1),
    G2_same_generation=list(abc_discordant(m,v3_B1,'G2'),v3_B1),
    Adult_B1_same_sample=list(adult,v3_B1),B2=list(adult,v3_B2))
  stopifnot(identical(specs$Adult_B1_same_sample[[1]]$sample_id,specs$B2[[1]]$sample_id))
  counts<-ov_bind(lapply(names(specs),function(k){z<-specs[[k]][[1]];data.frame(model=k,N=nrow(z),cases=sum(z$fatty_liver),
    controls=sum(z$fatty_liver==0),families=nlevels(z$family_id),G0=sum(z$generation=='G0'),G1=sum(z$generation=='G1'),G2=sum(z$generation=='G2'),minors=sum(z$age<18))}))
  v3_table(counts,'B_Discordant_family_feasibility_EN',out)
  overlap<-ov_bind(lapply(names(specs),function(k){z<-specs[[k]][[1]];ov_bind(lapply(split(z,list(z$generation,z$fatty_liver),drop=TRUE),function(a)
    data.frame(model=k,generation=as.character(a$generation[1]),steatosis=a$fatty_liver[1],N=nrow(a),age_min=min(a$age),age_Q1=quantile(a$age,.25),age_median=median(a$age),age_Q3=quantile(a$age,.75),age_max=max(a$age))))}))
  v3_table(overlap,'B_Age_generation_overlap_EN',out)
  fits<-lapply(names(specs),function(k)abc_B_fit(x,specs[[k]][[1]],specs[[k]][[2]],k,out,cache))
  t<-ov_bind(lapply(fits,`[[`,'results'));s<-ov_bind(lapply(fits,`[[`,'summary'));s$global_P_BH<-p.adjust(s$global_P_omnibus,'BH')
  v3_table(t,'B_All_within_family_LDM_SGB_results_EN',out);v3_table(t[t$BH_lt_005,],'B_BH_supported_SGBs_EN',out);v3_table(s,'B_Global_and_taxon_summary_EN',out)
  list(results=t,summary=s,feasibility=counts,overlap=overlap)
}

abc_C_sets <- function(x,cfg,out) {
  path<-file.path(cfg$project_root,'literature_reviews/ABC_primary_20261027/Asnicar2026_Supplementary_Tables.xlsx')
  tab<-suppressMessages(readxl::read_excel(path,sheet='S5',skip=3,.name_repair='unique'))
  rank<-data.frame(SGB=tab$SGB,health_rank=as.numeric(tab$HEALTH_ranks));stopifnot(nrow(rank)==661,!anyNA(rank),!anyDuplicated(rank$SGB))
  rank<-rank[order(rank$health_rank,rank$SGB),];sel<-rbind(transform(head(rank,50),set='Favourable'),transform(tail(rank,50),set='Unfavourable'))
  check<-readxl::read_excel(path,sheet='S7',skip=3,col_names=FALSE,.name_repair='unique')[[1]];check<-as.character(check);check<-check[grepl('^SGB[0-9]+(_group)?$',check)]
  stopifnot(length(check)==100,setequal(sel$SGB,check),!anyDuplicated(sel$SGB))
  ft<-x$features[match(rownames(x$relative),x$features$feature_id),];stopifnot(!anyDuplicated(ft$official_id[!is.na(ft$official_id)]))
  j<-match(sel$SGB,ft$official_id);sel$feature_id<-ft$feature_id[j];sel$mapping<-ifelse(is.na(j),'Absent from processed feature manifest','Exact official SGB/group identifier')
  sel$observed_prevalence<-ifelse(is.na(j),0,rowMeans(x$relative>0)[j]);sel$mean_relative_abundance<-ifelse(is.na(j),0,rowMeans(x$relative)[j])
  v3_table(sel,'C_Frozen_external_100_SGB_mapping_EN',out)
  v3_table(data.frame(set=c('Favourable','Unfavourable'),published_members=50,mapped_members=sapply(c('Favourable','Unfavourable'),function(s)sum(!is.na(sel$feature_id[sel$set==s])))), 'C_External_set_coverage_EN',out)
  d<-data.frame(sample_id=colnames(x$relative))
  for(s in c('Favourable','Unfavourable')){ids<-na.omit(sel$feature_id[sel$set==s]);A<-x$relative[ids,,drop=FALSE]
    d[[paste0(s,'_richness')]]<-colSums(A>0);d[[paste0(s,'_abundance')]]<-colSums(A)}
  stopifnot(all(d$Favourable_richness<=50),all(d$Unfavourable_richness<=50),all(d$Favourable_abundance+d$Unfavourable_abundance<=1+1e-10))
  saveRDS(d,file.path(out,'private/external_readouts.rds'));list(readouts=d,mapping=sel,path=path)
}

abc_C <- function(x,cfg,out,cache) {
  xx<-readRDS(tf_source(cfg$project_root)$files[1]);sets<-abc_C_sets(xx,cfg,out);m<-xx$metadata
  a<-match(m$sample_id,sets$readouts$sample_id);stopifnot(!anyNA(a));m<-cbind(m,sets$readouts[a,-1])
  b1<-v3_cc(m,c('fatty_liver','family_id',v3_B1));b2<-m$age>=18&v3_cc(m,c('fatty_liver','family_id',v3_B2))
  src<-m$age>=18&v3_cc(m,c('fatty_liver','family_id','age','sex','BMI18'))
  specs<-list(B1=list(b1,v3_B1,FALSE),Adult_B1_same_sample=list(b2,v3_B1,FALSE),B2=list(b2,v3_B2,FALSE),
    Source_age_sex_BMI_family_GEE=list(src,c('age','sex','BMI18'),FALSE),Source_age_sex_BMI_OLS=list(src,c('age','sex','BMI18'),TRUE))
  rows<-list();readouts<-c('Favourable_richness','Unfavourable_richness','Favourable_abundance','Unfavourable_abundance')
  for(key in names(specs)){
    sp<-specs[[key]];z<-droplevels(m[sp[[1]],]);z<-z[order(z$family_id,z$sample_id),];Z<-abc_design(z,sp[[2]],age_spline=!grepl('^Source_',key))$Z
    for(v in readouts){rich<-grepl('richness',v);response<-if(rich)z[[v]]else z$fatty_liver
      exposure<-if(rich)z$fatty_liver else asin(sqrt(z[[v]]))
      if(sp[[3]]){dd<-data.frame(y=response,exposure=exposure,Z);fit<-lm(y~.-1,data=dd);co<-coef(summary(fit))['exposure',];b<-unname(co[1]);se<-unname(co[2]);crit<-qt(.975,df.residual(fit))
        rr<-data.frame(beta=b,SE=se,CI_low=b-crit*se,CI_high=b+crit*se,P_raw=unname(co[4]))
      }else{ans<-abc_gee(response,exposure,Z,z);rr<-ans$row;fit<-ans$fit}
      saveRDS(fit,file.path(out,'private',paste0(key,'_',v,'_fit.rds')))
      rows[[paste(key,v)]]<-cbind(model=key,readout=v,N=nrow(z),cases=sum(z$fatty_liver),families=nlevels(z$family_id),minors=sum(z$age<18),
        outcome=if(rich)'Detected SGB count'else'Personal steatosis (0/1)',
        predictor=if(rich)'Personal steatosis (0/1)'else'AST cumulative relative abundance',rr)
    }
  }
  t<-ov_bind(rows);t$P_BH<-ave(t$P_raw,t$model,FUN=function(p)p.adjust(p,'BH'));t$BH_lt_005<-t$P_BH<.05
  t$effect_unit<-ifelse(grepl('richness',t$readout),'SGB count: steatosis minus no steatosis','Steatosis probability difference per 1 AST unit; linear, not logistic')
  t$display_beta<-ifelse(grepl('richness',t$readout),t$beta,t$beta*10)
  t$display_CI_low<-ifelse(grepl('richness',t$readout),t$CI_low,t$CI_low*10);t$display_CI_high<-ifelse(grepl('richness',t$readout),t$CI_high,t$CI_high*10)
  t$display_unit<-ifelse(grepl('richness',t$readout),'SGB count difference','Steatosis probability difference (percentage points) per 0.1 AST unit')
  v3_table(t,'C_External_health_set_associations_EN',out)
  desc<-ov_bind(lapply(readouts,function(v)ov_bind(lapply(0:1,function(e){z<-m[b1&m$fatty_liver==e,];a<-z[[v]];data.frame(readout=v,steatosis=e,N=nrow(z),median=median(a),Q1=quantile(a,.25),Q3=quantile(a,.75),mean=mean(a),SD=sd(a))}))))
  v3_table(desc,'C_External_set_distributions_EN',out);saveRDS(m[b1,],file.path(out,'private/B1_readout_plot_data.rds'))
  fl_write_csv(data.frame(path=sets$path,md5=unname(tools::md5sum(sets$path))),file.path(out,'audit/frozen_external_source.csv'))
  list(results=t,summary=desc,mapping=sets$mapping)
}

abc_reports <- function(r,which,cfg,out) {
  plan<-readLines(file.path(cfg$project_root,'analysis_plan/adopted_ABC_20261027.md'),warn=FALSE,encoding='UTF-8')
  writeLines(plan,file.path(out,'reports/Methods_and_references_EN.md'),useBytes=TRUE)
  tab<-if(which=='A')r$summary else if(which=='B')r$summary else r$results
  zh<-c(paste0('# ',which,' 分析结果'),'',if(which=='A')'在自身无超声脂肪肝的年轻成员中，分析观察年长一代脂肪肝与属水平菌群。'else if(which=='B')
    '仅在同时有患病与未患病成员的主家系内比较个人脂肪肝与菌群，使用官方LDM：家系固定效应＋家系内置换。'else
    '固定使用Asnicar原文健康排名最有利／最不利各50种菌；检验个人超声脂肪肝与菌集合的关系。',
    '',capture.output(print(tab,row.names=FALSE)), '',
    '既往27号SGB阴性结果保留。本分析与原文的人群、疾病、测序类型和协变量存在明确适配，不声称是原文全流程逐字复现。',
    '主家系编号未核实直接亲缘或共同居住；不推断传播、因果、中介或未来发病风险。',
    '原始RDS保留；仅汇总表和图片公开至结果目录。完整方法与参考文献见Methods_and_references_EN.md。')
  writeLines(zh,file.path(out,'reports/Results_summary_zh.md'),useBytes=TRUE)
  capture.output(sessionInfo(),file=file.path(out,'audit/sessionInfo.txt'))
}

abc_publish <- function(out,cfg,which) {
  if(file.exists(file.path(cfg$project_root,'config/results_viewer.json'))){
    gallery<-new.env(parent=baseenv());sys.source(file.path(cfg$project_root,'R/result_gallery.R'),envir=gallery)
    return(gallery$fl_publish_organized_results(cfg))
  }
  viewer<-file.path(cfg$project_root,'Results_表格与图片/ABC_followup_20261027',paste0(which,'_',cfg$output_module))
  dir.create(viewer,recursive=TRUE,showWarnings=FALSE)
  files<-unlist(lapply(c('tables','figures','reports'),function(p)list.files(file.path(out,p),recursive=TRUE,full.names=TRUE)))
  rel<-substring(files,nchar(out)+2);dest<-file.path(viewer,rel)
  for(i in seq_along(files)){dir.create(dirname(dest[i]),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(files[i],dest[i],overwrite=TRUE))}
  stopifnot(identical(unname(tools::md5sum(files)),unname(tools::md5sum(dest))))
  writeLines(c('<!doctype html><meta charset="utf-8"><style>body{font:15px Arial;margin:30px}img{max-width:100%;width:1100px}</style>',
    paste0('<h1>',which,' analysis</h1>'),'<ul>',paste0('<li><a href="',gsub('\\','/',rel,fixed=TRUE),'">',basename(rel),'</a></li>'),'</ul>',
    paste0('<img src="',gsub('\\','/',rel[grepl('\\.png$',rel)],fixed=TRUE),'">')),file.path(viewer,'index.html'),useBytes=TRUE)
  parent<-dirname(viewer);dirs<-list.dirs(parent,recursive=FALSE,full.names=FALSE)
  dirs<-dirs[file.exists(file.path(parent,dirs,'index.html'))]
  writeLines(c('<!doctype html><meta charset="utf-8"><style>body{font:16px Arial;margin:35px;line-height:1.8}</style>',
    '<h1>ABC follow-up: hepatic steatosis and the microbiome</h1>',
    '<p>A: unaffected younger members, genus/community analysis. B: within-family personal-steatosis contrast. C: externally predefined health-ranked SGB sets.</p>',
    '<p>Published methods with documented cohort adaptations. Original module-27 negative results are retained.</p><ul>',
    paste0('<li><a href="',dirs,'/index.html">',dirs,'</a></li>'),'</ul>'),file.path(parent,'index.html'),useBytes=TRUE)
  viewer
}

abc_main <- function(which,args=character(),root=Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')) {
  stopifnot(which%in%c('A','B','C'))
  if(which=='C'){sys.source(file.path(root,'R/external_health_ranking.R'),envir=.GlobalEnv);return(ehr_main(args,root))}
  if(getRversion()<'4.5'){entry<-c(A='28_family_genus_community.R',B='29_within_family_steatosis.R',C='30_external_health_sets.R')[[which]]
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts',entry)),vapply(args,shQuote,character(1))))
    if(status!=0)stop('Project R ABC entry failed');return(invisible(NULL))}
  cfg<-abc_setup(root);key<-c(A='abc_genus',B='abc_within',C='abc_external')[[which]]
  cfg$output_module<-unname(cfg$output_modules[key]);cfg$run_id<-unname(cfg$run_ids[key]);cfg$run_id_explicit<-FALSE
  if('--check'%in%args){cat('ABC ',which,': dependencies and source check\n');abc_source(root);return(invisible(NULL))}
  if(any(!args%in%c('--run','--figures-only')))stop('Supported modes: --check --run --figures-only')
  module<-file.path(root,'outputs',cfg$output_module);cache<-file.path(module,'_checkpoints/private');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  if('--figures-only'%in%args){id<-readLines(file.path(module,'LATEST_RUN.txt'),1);out<-file.path(module,'current',id);r<-readRDS(file.path(out,'private/analysis_result.rds'))
    abc_figures(r,which,cfg,out);abc_publish(out,cfg,which);return(invisible(out))}
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(p in c('private','audit','tables','figures','reports'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
  writeLines(c('RUNNING',abc_protocol),file.path(out,'run_status.txt'))
  tryCatch({src<-abc_source(root);fl_write_csv(data.frame(path=src$file,md5=src$md5),file.path(out,'audit/input_fingerprints.csv'))
    r<-switch(which,A=abc_A(src$x,cfg,out,cache),B=abc_B(src$x,cfg,out,cache),C=abc_C(src$x,cfg,out,cache))
    saveRDS(r,file.path(out,'private/analysis_result.rds'));abc_reports(r,which,cfg,out);abc_validate(r,which,cfg,out);abc_figures(r,which,cfg,out)
    stopifnot(unname(tools::md5sum(src$file))==src$md5)
    writeLines(c('COMPLETED',abc_protocol),file.path(out,'run_status.txt'));target<-fl_output_complete_run(ctx)
    writeLines(cfg$run_id,file.path(module,'LATEST_RUN.txt'));abc_publish(target,cfg,which);tf_log('COMPLETED ABC ',which,': ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
