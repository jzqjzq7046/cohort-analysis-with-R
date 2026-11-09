# Revised modules 25/26. Inert until tf_main(). Module 24 is read-only.
# 25: official MaAsLin2 LOG model effect comparisons (Ma 2024; Mei 2024).
# 26: Gaussian independence GEE, clustered by primary family (Liang & Zeger 1986).
# The cross-generation exposure definition is a cohort-specific adaptation.
tf_protocol <- 'taxa_effect_comparison_and_family_regression_v2_20261026'
tf_log <- function(...) {cat(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),paste0(...,collapse=''),'\n');flush.console()}
tf_read <- function(path) read.csv(path,check.names=FALSE,fileEncoding='UTF-8-BOM',stringsAsFactors=FALSE)
tf_bool <- function(z) !is.na(z)&as.logical(z)

tf_source <- function(root) {
  base<-file.path(root,'outputs/24_Two_generation_reanalysis')
  id<-readLines(file.path(base,'LATEST_RUN.txt'),n=1,warn=FALSE)
  stopifnot(length(id)==1L,grepl('^[A-Za-z0-9_.-]+$',id))
  src<-file.path(base,'current',id)
  if(!file.exists(file.path(src,'run_status.txt'))||readLines(file.path(src,'run_status.txt'),n=1)!='COMPLETED')
    stop('The adopted two-generation source must be COMPLETED; no old-cohort fallback is allowed.')
  files<-c(file.path(src,'private/stage_prepare.rds'),
    file.path(src,'tables',paste0('Table6_SGB_',c('B1','Adult_B1_same_sample','B2'),'_all_results_EN.csv')))
  if(any(!file.exists(files)))stop('Missing adopted source files: ',paste(files[!file.exists(files)],collapse='; '))
  list(path=src,run_id=id,files=files,md5=tools::md5sum(files))
}

tf_adjust <- function(t,family) {
  stopifnot(!anyDuplicated(t$feature_id))
  tested<-!grepl('^FILTERED',t$status);valid<-tested&tf_bool(t$analysis_valid)&is.finite(t$P_raw)&t$P_raw>=0&t$P_raw<=1
  p<-rep(1,sum(tested));p[valid[tested]]<-t$P_raw[valid]
  t$P_BH<-NA_real_;q<-p.adjust(p,method='BH');t$P_BH[which(tested)[valid[tested]]]<-q[valid[tested]]
  t$FDR_family<-family;t$FDR_tests<-sum(tested)
  t$FDR_lt_005<-valid&!is.na(t$P_BH)&t$P_BH<.05
  t$direction<-ifelse(!valid,'Not estimable',ifelse(t$beta_log2<0,'Lower','Higher'))
  t
}

tf_baseline <- function(src) {
  nms<-c('B1','Adult_B1_same_sample','B2')
  setNames(lapply(nms,function(nm)tf_adjust(tf_read(file.path(src$path,'tables',paste0('Table6_SGB_',nm,'_all_results_EN.csv'))),paste0('Disease_',nm,'_all_tested_SGBs'))),nms)
}

tf_load <- function(src) {
  x<-readRDS(src$files[1]);m<-x$metadata
  stopifnot(identical(x$protocol,v3_protocol),!anyDuplicated(m$sample_id),
    identical(colnames(x$relative),as.character(m$sample_id)),!anyNA(m$family_id),
    max(abs(colSums(x$relative)-1))<1e-10,all(is.finite(x$relative)),all(x$relative>=0),
    all(rowSums(table(m$family_id,m$generation)>0)>=2))
  x
}

tf_checkpoint <- function(path,signature,compute) {
  if(file.exists(path)){z<-readRDS(path);if(identical(z$signature,signature))return(z$value)}
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  z<-compute()
  tmp<-paste0(path,'.tmp');saveRDS(list(signature=signature,value=z),tmp)
  if(file.exists(path))unlink(path)
  if(!file.rename(tmp,path))stop('Cannot adopt checkpoint: ',path)
  z
}

tf_official <- function(d,ids,model,cfg,cache,signature) {
  pieces<-list();batches<-split(ids,ceiling(seq_along(ids)/20L))
  for(i in seq_along(batches)){
    keys<-batches[[i]];cp<-file.path(cache,sprintf('batch_%03d.rds',i))
    piece<-tf_checkpoint(cp,paste(signature,paste(keys,collapse=';')),function(){
      tf_log(model,': batch ',i,'/',length(batches),' (',length(keys),' SGBs)')
      A<-d$relative[,keys,drop=FALSE];expected<-vapply(keys,function(k)get('LOG',asNamespace('Maaslin2'))(A[,k]),numeric(nrow(A)))
      rownames(expected)<-rownames(A);colnames(expected)<-keys
      folder<-file.path(cache,sprintf('official_%03d',i));dir.create(folder,recursive=TRUE,showWarnings=FALSE)
      console<-file(file.path(folder,'console.log'),'wt');sink(console)
      ans<-tryCatch(Maaslin2::Maaslin2(input_data=as.data.frame(A),input_metadata=d$design,output=file.path(folder,'fit'),
        fixed_effects=c(d$bases,d$nuisance),random_effects='family_id',normalization='NONE',transform='LOG',standardize=FALSE,
        min_prevalence=0,min_abundance=0,min_variance=0,analysis_method='LM',cores=1,correction='BH',
        plot_heatmap=FALSE,plot_scatter=FALSE,save_models=TRUE),finally={sink();close(console)})
      # Per-batch package q-values are intentionally discarded.
      ov_extract(ans,keys,d,model,expected)
    })
    pieces[[i]]<-piece
  }
  list(effects=ov_bind(lapply(pieces,`[[`,'effects')),diagnostics=ov_bind(lapply(pieces,`[[`,'diagnostics')))
}

tf_no_antibiotics <- function(x,src,cfg,out,cache) {
  mask<-v3_micro_samples(x)$B1$mask & !is.na(x$metadata$abx_recent)&x$metadata$abx_recent==0
  d<-v3_design(x,mask,v3_B1);ids<-d$filter$feature_id[d$filter$retained]
  if(!length(ids))stop('No eligible SGBs in the antibiotic exclusion sample.')
  sig<-tf_signature(src,cfg,'no_antibiotics',list(samples=d$metadata$sample_id,filter=d$filter))
  e<-tf_official(d,ids,'No_recent_antibiotics',cfg,file.path(cache,'no_antibiotics'),sig)
  t<-ov_expand(e$effects,d,x$candidates,'No_recent_antibiotics')
  f<-x$features[match(t$feature_id,x$features$feature_id),,drop=FALSE]
  result<-data.frame(model='No_recent_antibiotics',feature_id=t$feature_id,SGB_id=sub('^t__','',t$feature_id),
    species_label=ov_species_label(f$taxonomy),N=nrow(d$metadata),cases=sum(d$metadata$fatty_liver),
    families=length(unique(d$metadata$family_id)),beta_log2=t$estimate,SE=t$SE,CI_low=t$CI_low,CI_high=t$CI_high,
    P_raw=t$p_raw,analysis_valid=t$analysis_valid,status=t$status,singular=t$singular)
  result<-tf_adjust(result,'Disease_no_antibiotics_all_tested_SGBs')
  v3_table(result,'Antibiotic_exclusion_all_SGBs_EN',out)
  fl_write_csv(d$mapping,file.path(out,'audit/antibiotic_exclusion_design.csv'))
  fl_write_csv(d$filter,file.path(out,'audit/antibiotic_exclusion_filter.csv'))
  fl_write_csv(e$diagnostics,file.path(out,'audit/antibiotic_exclusion_diagnostics.csv'))
  result
}

tf_signature <- function(src,cfg,kind,extra=NULL) digest::digest(list(protocol=tf_protocol,
  parent_MD5=unname(src$md5),kind=kind,extra=extra,
  source_MD5=unname(tools::md5sum(file.path(cfg$project_root,c('R/taxa_followup.R',
    'R/reanalysis_20261024.R','R/reanalysis_numerical_helpers.R')))),
  Maaslin2=as.character(packageVersion('Maaslin2'))),algo='sha256')

tf_effect_units <- function(t) {
  t$abundance_ratio<-2^t$beta_log2;t$ratio_CI_low<-2^t$CI_low;t$ratio_CI_high<-2^t$CI_high
  t$effect_unit<-'Log2 zero-replaced relative abundance: steatosis minus no steatosis'
  t$ratio_unit<-'Ratio on the zero-replaced relative-abundance scale; not absolute microbial load'
  t
}

tf_compare <- function(a,b,label) {
  z<-merge(a[,c('feature_id','beta_log2','P_BH','analysis_valid')],b[,c('feature_id','beta_log2','P_BH','analysis_valid')],by='feature_id',suffixes=c('_x','_y'))
  z<-z[tf_bool(z$analysis_valid_x)&tf_bool(z$analysis_valid_y)&is.finite(z$beta_log2_x)&is.finite(z$beta_log2_y),]
  z$comparison<-label;z$same_direction<-sign(z$beta_log2_x)==sign(z$beta_log2_y)
  summary<-data.frame(comparison=label,common_valid_SGBs=nrow(z),same_direction=sum(z$same_direction),
    same_direction_percent=100*mean(z$same_direction),Spearman_rho=cor(z$beta_log2_x,z$beta_log2_y,method='spearman'))
  list(data=z,summary=summary)
}

tf_robustness <- function(x,src,cfg,out,cache) {
  base<-lapply(tf_baseline(src),tf_effect_units)
  anti<-tf_effect_units(tf_no_antibiotics(x,src,cfg,out,cache))
  tables<-c(base,list(No_recent_antibiotics=anti))
  for(nm in names(tables))v3_table(tables[[nm]],paste0('Disease_',nm,'_all_SGBs_EN'),out)
  adult<-tf_compare(base$Adult_B1_same_sample,base$B2,'Adult B1 versus adult B2: same participants')
  antibiotics<-tf_compare(base$B1,anti,'Overall B1 versus antibiotic exclusion: different samples')
  comparisons<-ov_bind(list(adult$data,antibiotics$data));comparison_summary<-ov_bind(list(adult$summary,antibiotics$summary))
  v3_table(comparisons,'Effect_comparison_all_common_SGBs_EN',out)
  v3_table(comparison_summary,'Effect_comparison_summary_EN',out)
  summary<-ov_bind(lapply(names(tables),function(nm){t<-tables[[nm]];data.frame(model=nm,N=t$N[1],cases=t$cases[1],
    tested=max(t$FDR_tests),valid=sum(tf_bool(t$analysis_valid)),raw_P_lt_005=sum(t$P_raw<.05&tf_bool(t$analysis_valid),na.rm=TRUE),
    BH_lt_005=sum(t$FDR_lt_005),BH_lower=sum(t$FDR_lt_005&t$beta_log2<0,na.rm=TRUE),BH_higher=sum(t$FDR_lt_005&t$beta_log2>0,na.rm=TRUE))}))
  v3_table(summary,'Robustness_model_summary_EN',out)
  list(tables=tables,comparisons=comparisons,comparison_summary=comparison_summary,summary=summary)
}

# Frisch-Waugh-Lovell plus family-cluster HC0 sandwich. Runtime equivalence
# against geepack::geeglm Gaussian/independence/san.se is mandatory below.
tf_cluster_regression <- function(X,Y,W,family_id) {
  stopifnot(nrow(X)==nrow(Y),nrow(W)==nrow(X),length(family_id)==nrow(X),all(is.finite(X)),all(is.finite(Y)),all(is.finite(W)))
  qw<-qr(W,tol=1e-9);W<-W[,sort(qw$pivot[seq_len(qw$rank)]),drop=FALSE];qw<-qr(W)
  XR<-qr.resid(qw,X);YR<-qr.resid(qw,Y);ss<-colSums(XR^2)
  valid<-is.finite(ss)&ss>1e-10;beta<-rep(NA_real_,ncol(X));beta[valid]<-colSums(XR[,valid,drop=FALSE]*YR[,valid,drop=FALSE])/ss[valid]
  errors<-YR[,valid,drop=FALSE]-sweep(XR[,valid,drop=FALSE],2,beta[valid],'*')
  scores<-rowsum(XR[,valid,drop=FALSE]*errors,group=family_id,reorder=FALSE)
  SE<-rep(NA_real_,ncol(X));SE[valid]<-sqrt(colSums(scores^2))/ss[valid]
  valid<-valid&is.finite(beta)&is.finite(SE)&SE>0
  list(beta=beta,SE=SE,P=2*pnorm(-abs(beta/SE)),valid=valid,W=W,
    residual_SD=sqrt(colSums(YR^2)/nrow(Y)),predictor_residual_SS=ss)
}

tf_verify_gee <- function(X,Y,fit,family_id,comparison,model,out) {
  W<-fit$W;colnames(W)<-paste0('w',seq_len(ncol(W)))
  checks<-lapply(head(which(fit$valid),3L),function(j){
    md<-data.frame(y=Y[,j],x=X[,j],W,cluster=as.integer(factor(family_id)))
    formula<-reformulate(c('x',colnames(W)),response='y',intercept=FALSE)
    ref<-geepack::geeglm(formula,data=md,id=cluster,family=gaussian(),corstr='independence',std.err='san.se')
    stopifnot(ref$geese$error==0L)
    beta<-unname(coef(ref)['x']);SE<-sqrt(vcov(ref)['x','x']);P<-2*pnorm(-abs(beta/SE))
    err<-max(abs(c(beta-fit$beta[j],SE-fit$SE[j],P-fit$P[j])))
    if(!is.finite(err)||err>1e-6)stop('Family GEE equivalence failed: ',comparison,' ',model,' ',colnames(X)[j],' error=',err)
    data.frame(comparison=comparison,model=model,feature_id=colnames(X)[j],beta_reference=beta,SE_reference=SE,
      beta_vectorized=fit$beta[j],SE_vectorized=fit$SE[j],maximum_error=err,status='PASS')
  })
  stopifnot(length(checks)>0)
  fl_write_csv(ov_bind(checks),file.path(out,'audit',paste0('GEE_backend_equivalence_',comparison,'_',model,'.csv')))
}

tf_family_analysis <- function(x,src,cfg,out,cache) {
  ids<-x$candidates;A<-x$relative[ids,,drop=FALSE]
  replacement<-apply(A,1,function(z)if(any(z>0))min(z[z>0])/2 else NA_real_)
  stopifnot(all(is.finite(replacement)),all(replacement>0))
  Z<-log2(A);zero<-which(A==0,arr.ind=TRUE);Z[zero]<-log2(replacement[zero[,1]])
  fl_write_csv(data.frame(feature_id=ids,zero_replacement=replacement),file.path(out,'audit/fixed_zero_replacements.csv'))
  all<-list();samples<-list()
  for(g in v3_pairs){comparison<-paste(g,collapse='-');m<-v3_scope(x$metadata,g)
    m<-v3_scope(m[v3_cc(m,v3_P2),,drop=FALSE],g);m<-m[order(m$family_id,m$sample_id),,drop=FALSE]
    a<-which(m$generation==g[1]);b<-which(m$generation==g[2]);edge<-which(outer(as.character(m$family_id[a]),as.character(m$family_id[b]),'=='),arr.ind=TRUE)
    edge<-edge[order(m$family_id[a[edge[,1]]],m$sample_id[b[edge[,2]]],m$sample_id[a[edge[,1]]]),,drop=FALSE]
    older<-a[edge[,1]];younger<-b[edge[,2]];fam<-as.character(m$family_id[older]);stopifnot(all(fam==as.character(m$family_id[younger])))
    index<-match(m$sample_id,colnames(A));pvA<-rowMeans(A[,index[a],drop=FALSE]>0);pvB<-rowMeans(A[,index[b],drop=FALSE]>0)
    eligible<-pvA>=.10&pvB>=.10;kept<-ids[eligible]
    X<-t(Z[kept,index[older],drop=FALSE]);Y<-t(Z[kept,index[younger],drop=FALSE])
    oa<-v3_nuisance(m[a,,drop=FALSE],v3_P2);ob<-v3_nuisance(m[b,,drop=FALSE],v3_P2)
    WA<-oa$Z[edge[,1],,drop=FALSE];WB<-ob$Z[edge[,2],,drop=FALSE]
    colnames(WA)<-paste0('older_',colnames(WA));colnames(WB)<-paste0('younger_',colnames(WB))
    W2<-cbind(WA,WB)
    fl_write_csv(data.frame(endpoint=c(rep('older',ncol(WA)),rep('younger',ncol(WB))),column=c(colnames(WA),colnames(WB))),file.path(out,'audit',paste0('Endpoint_covariates_',comparison,'.csv')))
    saveRDS(data.frame(older_id=m$sample_id[older],younger_id=m$sample_id[younger],family_id=fam),file.path(out,'private',paste0('Observed_pairs_',comparison,'.rds')))
    fl_write_csv(data.frame(feature_id=ids,older_prevalence=pvA,younger_prevalence=pvB,eligible=eligible),file.path(out,'audit',paste0('Filter_',comparison,'.csv')))
    samples[[comparison]]<-data.frame(comparison=comparison,participants=nrow(m),older_N=length(a),younger_N=length(b),families=length(unique(fam)),observed_within_family_pairs=nrow(edge),eligible_SGBs=length(kept))
    for(model in c('P1','P2')){
      W<-if(model=='P1')matrix(1,nrow(edge),1,dimnames=list(NULL,'Intercept'))else W2
      tf_log('Family abundance regression ',comparison,' ',model,': ',length(kept),' SGBs; ',nrow(edge),' observed pairs; ',length(unique(fam)),' families')
      fit<-tf_cluster_regression(X,Y,W,fam);tf_verify_gee(X,Y,fit,fam,comparison,model,out)
      j<-match(ids,kept);beta<-fit$beta[j];SE<-fit$SE[j];valid<-eligible&!is.na(j)&tf_bool(fit$valid[j])
      t<-data.frame(comparison=comparison,model=model,feature_id=ids,SGB_id=sub('^t__','',ids),eligible=eligible,
        older_prevalence=pvA,younger_prevalence=pvB,N=nrow(m),older_N=length(a),younger_N=length(b),families=length(unique(fam)),within_family_pairs=nrow(edge),
        beta_log2_per_log2=beta,SE=SE,CI_low=beta-qnorm(.975)*SE,CI_high=beta+qnorm(.975)*SE,P_raw=fit$P[j],
        analysis_valid=valid,status=ifelse(!eligible,'FILTERED_PREVALENCE',ifelse(valid,'ESTIMATED','NOT_ESTIMABLE')),
        nuisance_rank=ncol(fit$W),predictor_residual_SS=fit$predictor_residual_SS[j],
        younger_ratio_per_older_doubling=2^beta,ratio_CI_low=2^(beta-qnorm(.975)*SE),ratio_CI_high=2^(beta+qnorm(.975)*SE))
      t$species_label<-ov_species_label(x$features$taxonomy[match(ids,x$features$feature_id)])
      all[[paste(comparison,model)]]<-t
    }
  }
  t<-ov_bind(all);t$P_BH<-NA_real_;t$FDR_tests<-0L
  for(model in c('P1','P2')){
    k<-which(t$model==model&t$eligible);good<-tf_bool(t$analysis_valid[k])&is.finite(t$P_raw[k]);p<-rep(1,length(k));p[good]<-t$P_raw[k][good]
    q<-p.adjust(p,'BH');t$P_BH[k[good]]<-q[good];t$FDR_tests[t$model==model]<-length(k)
  }
  t$FDR_family<-paste0(t$model,'_all_eligible_SGBs_x_three_comparisons');t$FDR_lt_005<-is.finite(t$P_BH)&t$P_BH<.05
  t$effect_unit<-'Younger log2 zero-replaced abundance per one-unit older log2 abundance; family-cluster robust GEE'
  v3_table(t,'Taxon_cross_generation_GEE_all_results_EN',out)
  sample_table<-ov_bind(samples);v3_table(sample_table,'Family_pair_sample_sizes_EN',out)
  summary<-ov_bind(lapply(split(t,list(t$model,t$comparison),drop=TRUE),function(z)data.frame(model=z$model[1],comparison=z$comparison[1],
    eligible=sum(z$eligible),valid=sum(z$analysis_valid),BH_positive=sum(z$FDR_lt_005&z$beta_log2_per_log2>0,na.rm=TRUE),
    BH_negative=sum(z$FDR_lt_005&z$beta_log2_per_log2<0,na.rm=TRUE),raw_P_lt_005=sum(z$analysis_valid&z$P_raw<.05,na.rm=TRUE))))
  v3_table(summary,'Taxon_concordance_summary_EN',out)
  list(results=t,samples=sample_table,summary=summary)
}

tf_reports <- function(result,analysis,src,cfg,out) {
  common<-c('# Revised taxon follow-up analyses','',paste('Protocol:',tf_protocol),paste('Adopted module-24 source:',src$run_id),
    'Outcome: ultrasound-defined hepatic steatosis. MASLD is not analysed. Only the two relevant observed generations are required; minors remain eligible.',
    'Primary-family labels are observed grouping variables; direct parent-child and grandparent-grandchild relationships have not been verified.',
    'Every tested eligible feature is retained in the full tables, including unsuccessful estimates. BH correction uses all eligible tests, with invalid estimates conservatively counted as P=1; their displayed P/q remain unavailable.',
    'Display selection is separate from analysis filtering and does not define a new discovery or validation set.')
  detail<-if(analysis=='taxa_robustness')c('','## Module 25: disease association effect comparisons',
    'Read the completed overall B1 and same-adult-sample B1/B2 official MaAsLin2 1.24.1 results from module 24. No overall or adult models are refitted just for figures.',
    'Fit official MaAsLin2 LOG models after excluding recent antibiotics: log2 zero-replaced relative abundance ~ steatosis + B1 nuisance terms + (1 | primary family). Zeros use the official half-smallest-positive replacement within each model; positives are unchanged. Ratios 2^beta refer to this zero-replaced abundance scale, not absolute bacterial load.',
    'B1: generation, generation-specific age splines, sex, preparation batch, sequencing pool, season sine/cosine, recent antibiotics. B2: B1 plus smoking, alcohol, diet score, physical activity, total energy and continuous adult BMI. Adult B1 and B2 use exactly the same complete adult sample. The antibiotic-excluded model drops the now constant antibiotic column.',
    'Use the existing module-24 candidate union and model-specific prevalence >=10%. Correct across all tested species separately for each declared model, never within batches.',
    'Figure25a/b compare all common valid coefficients, following the sensitivity coefficient-comparison presentation in Mei et al. 2024, Extended Data Figs. 4c and 5. Spearman rho is descriptive; these are same-cohort sensitivity analyses, not independent replication.',
    'The companion overall-abundance figure displays at most 16 valid overall B1 BH-supported SGBs, ranked by BH q, raw P and ID. The full tables contain all candidates. No selected-feature bootstrap or four-model-intersection discovery rule is used.',
    'Reference adaptation: Ma et al. reduced/full covariate strategy and MaAsLin2; Mei et al. sensitivity effect comparisons. Their exact populations, thresholds and covariates differ from this approved cohort protocol.')else c('','## Module 26: single-taxon cross-generation abundance concordance',
    'For each G0-G1, G1-G2 and G0-G2 comparison, enumerate each observed within-primary-family cross-generation dyad once. Individuals may recur in dyads. No between-family dyad distance is calculated.',
    'Exposure: older member log2 zero-replaced relative abundance of the SGB. Outcome: younger member log2 zero-replaced relative abundance of the same SGB. Fix per-SGB zero replacement at half its minimum positive abundance in the full adopted baseline, identically across both generation roles and all models.',
    'P1: younger abundance ~ older abundance. P2 adds both members age splines, sex, preparation batch, pool, season sine/cosine and recent antibiotics. P1 and P2 use identical P2-complete endpoints and observed pairs. Each SGB requires >=10% detection in EACH endpoint generation. No disease association P values enter filtering.',
    'Gaussian identity-link GEE with independence working correlation, clustering by primary family and the uncorrected family-cluster sandwich variance. Cluster identifiers are contiguous. All dyads enter with equal weight; therefore the estimand describes observed dyads, and families contributing more pairs have more weight.',
    'Vectorized Frisch-Waugh-Lovell estimates and HC0 cluster scores are checked against actual geepack::geeglm Gaussian independence robust estimates for the first three valid features in every comparison/model; maximum absolute beta/SE/P disagreement must be <1e-6.',
    'A positive beta indicates abundance concordance. 2^beta is the younger abundance ratio per doubling of older abundance on the zero-replaced relative-abundance scale; it is not an OR, absolute agreement, heritability, causal transmission or verified parent-child inheritance. Zero replacement and compositional constraints can influence estimates.',
    'BH pools every eligible SGB x all three comparisons, separately for P1 and P2. The plotted 16 SGBs are selected for display by minimum adjusted P2 q across comparisons, then ID; all eligible SGBs were analysed.',
    'Method provenance: the estimating equations and cluster sandwich follow Liang & Zeger; geepack provides the reference implementation. Applying this established regression to older/younger SGB abundance is our cohort-specific adaptation. Familial microbiome publications motivate the question but did not use this exact exposure/formula. It is not presented as verbatim replication of Valles-Colomer or Gacesa.')
  refs<-c('','## References',
    'Ma W, et al. Gut microbiome composition and metabolic activity in women with diverticulitis. Nature Communications. 2024;15:3612. doi:10.1038/s41467-024-47859-4.',
    'Mei Z, et al. Strain-specific gut microbial signatures in type 2 diabetes identified in a cross-cohort analysis of 8,117 metagenomes. Nature Medicine. 2024;30:2265-2276. doi:10.1038/s41591-024-03067-7.',
    'Mallick H, et al. Multivariable association discovery in population-scale meta-omics studies. PLoS Computational Biology. 2021;17:e1009442. doi:10.1371/journal.pcbi.1009442.',
    'Liang KY, Zeger SL. Longitudinal data analysis using generalized linear models. Biometrika. 1986;73:13-22. doi:10.1093/biomet/73.1.13.',
    'Halekoh U, Hojsgaard S, Yan J. The R package geepack for generalized estimating equations. Journal of Statistical Software. 2006;15(2):1-11. doi:10.18637/jss.v015.i02.',
    'Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. JRSS B. 1995;57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.')
  writeLines(c(common,detail,refs),file.path(out,'reports/Methods_and_references_EN.md'),useBytes=TRUE)
}

tf_publish <- function(target,cfg,analysis) {
  if(file.exists(file.path(cfg$project_root,'config/results_viewer.json'))){
    gallery<-new.env(parent=baseenv());sys.source(file.path(cfg$project_root,'R/result_gallery.R'),envir=gallery)
    return(gallery$fl_publish_organized_results(cfg))
  }
  viewer<-file.path(cfg$project_root,'Results_表格与图片','Followup_taxa_revised_20261026',cfg$output_module,basename(target))
  dir.create(viewer,recursive=TRUE,showWarnings=FALSE)
  public<-list.files(target,recursive=TRUE,full.names=FALSE);public<-gsub('\\','/',public,fixed=TRUE)
  public<-public[grepl('^(tables|figures|reports)/',public)]
  records<-lapply(public,function(p){s<-file.path(target,p);d<-file.path(viewer,p);dir.create(dirname(d),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(s,d,overwrite=TRUE))stop('Cannot publish ',p)
    hash<-unname(tools::md5sum(s));stopifnot(unname(tools::md5sum(d))==hash);data.frame(file=p,md5=hash)})
  fl_write_csv(ov_bind(records),file.path(viewer,'file_manifest.csv'))
  html<-c('<!doctype html><meta charset="utf-8"><style>body{font:15px Arial;margin:28px;line-height:1.6}img{max-width:100%;width:1100px}</style>',
    paste0('<h1>',if(analysis=='taxa_robustness')'Taxon association: sensitivity effect comparisons'else'Single-taxon cross-generation abundance concordance','</h1>'),
    '<ul>',paste0('<li><a href="',public,'">',public,'</a></li>'),'</ul>',paste0('<img src="',public[grepl('\\.png$',public)],'">'))
  writeLines(html,file.path(viewer,'index.html'),useBytes=TRUE);tf_log('Viewing copies: ',viewer)
  invisible(viewer)
}

tf_main <- function(analysis,args=character(),root=Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')) {
  stopifnot(analysis%in%c('taxa_robustness','taxon_family_sharing'))
  entry<-if(analysis=='taxa_robustness')'25_taxa_robustness.R'else'26_taxon_family_sharing.R'
  if(getRversion()<'4.5'){
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts',entry)),vapply(args,shQuote,character(1))))
    if(status!=0)stop('Project-R follow-up failed');return(invisible(NULL))
  }
  .libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),file.path(R.home(),'library'),.libPaths()))
  for(f in c('config.R','R/common.R','R/entry_helpers.R','R/nature_export_approved.R','R/reanalysis_numerical_helpers.R','R/reanalysis_20261024.R','R/taxa_followup_figures.R'))sys.source(file.path(root,f),envir=.GlobalEnv)
  cfg<-fl_config(root);cfg$output_module<-unname(cfg$output_modules[analysis]);cfg$run_id<-unname(cfg$run_ids[analysis]);cfg$run_id_explicit<-FALSE
  deps<-c('Maaslin2','lme4','lmerTest','geepack','digest','ggplot2','patchwork','svglite','jsonlite')
  stopifnot(all(vapply(deps,requireNamespace,logical(1),quietly=TRUE)),as.character(packageVersion('Maaslin2'))=='1.24.1')
  if(!identical(as.character(Matrix::Matrix.Version()$abi),as.character(lme4:::.Matrix.abi.build.version)))stop('Matrix ABI mismatch; use project R --vanilla')
  src<-tf_source(root)
  if('--check'%in%args){cat('Static source/package check PASS. No fits or plots.\n');return(invisible(NULL))}
  if(any(!args%in%c('--run','--figures-only')))stop('Supported: --check, --run, --figures-only')
  module<-file.path(root,'outputs',cfg$output_module)
  if('--figures-only'%in%args){
    id<-readLines(file.path(module,'LATEST_RUN.txt'),1);out<-file.path(module,'current',id);cfg<-readRDS(file.path(out,'private/run_config.rds'))
    result<-readRDS(file.path(out,'private/analysis_result.rds'))
    Sys.setenv(PYTHONPATH=paste(c(file.path(root,'.figure_qa_python'),Sys.getenv('PYTHONPATH')),collapse=.Platform$path.sep))
    tf_figures(result,analysis,cfg,out);tf_publish(out,cfg,analysis);return(invisible(out))
  }
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(p in c('tables','figures','reports','audit','private'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
  sig<-tf_signature(src,cfg,analysis);cache<-file.path(module,'_checkpoints','private',substr(sig,1,24));dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  writeLines(c('RUNNING',tf_protocol,analysis),file.path(out,'run_status.txt'))
  saveRDS(cfg,file.path(out,'private/run_config.rds'));saveRDS(src,file.path(out,'private/source_pointer.rds'))
  fl_write_csv(data.frame(path=src$files,md5=unname(src$md5)),file.path(out,'audit/adopted_source_fingerprints.csv'))
  Sys.setenv(PYTHONPATH=paste(c(file.path(root,'.figure_qa_python'),Sys.getenv('PYTHONPATH')),collapse=.Platform$path.sep),OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4')
  tf_log('START ',analysis,'; ',cfg$run_id)
  tryCatch({
    x<-tf_load(src)
    fl_write_csv(x$features[match(x$candidates,x$features$feature_id),,drop=FALSE],file.path(out,'tables/SGB_taxonomy_mapping_EN.csv'))
    result<-tf_checkpoint(file.path(out,'private/analysis_checkpoint.rds'),sig,function(){
      if(analysis=='taxa_robustness')tf_robustness(x,src,cfg,out,cache)else tf_family_analysis(x,src,cfg,out,cache)})
    saveRDS(result,file.path(out,'private/analysis_result.rds'));tf_reports(result,analysis,src,cfg,out)
    tf_figures(result,analysis,cfg,out)
    stopifnot(identical(unname(src$md5),unname(tools::md5sum(src$files))))
    writeLines(c('COMPLETED',tf_protocol,analysis),file.path(out,'run_status.txt'))
    target<-fl_output_complete_run(ctx);writeLines(cfg$run_id,file.path(module,'LATEST_RUN.txt'));tf_publish(target,cfg,analysis)
    tf_log('COMPLETED: ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}

