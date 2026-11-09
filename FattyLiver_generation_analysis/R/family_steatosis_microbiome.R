# Official MaAsLin2 LOG framework (Mallick 2021); unaffected-family-background
# comparison follows Olivera 2025. Disease/labels/SGBs/covariates are explicit
# cohort applications, not an exact reproduction of Crohn's disease biology.
fp_protocol <- 'unaffected_younger_observed_older_steatosis_Maaslin2_20261027_v1'

fp_prepare <- function(src,cfg,out) {
  x<-readRDS(src$files[1]);stopifnot(identical(x$protocol,v3_protocol))
  path<-file.path(cfg$project_root,'data/processed/processed_microbiome_v1/microbiome_processed.rds')
  manifest<-tf_read(file.path(dirname(path),'dataset_manifest.csv'))
  stopifnot(unname(tools::md5sum(path))==manifest$md5[manifest$file=='microbiome_processed.rds'])
  d<-readRDS(path);m<-d$metadata;cl<-x$clinical
  stopifnot(!anyDuplicated(m$sample_id),!anyDuplicated(m$participant_id),!anyDuplicated(cl$participant_id),
    identical(colnames(d$full$counts),as.character(m$sample_id)))
  k<-match(as.character(m$participant_id),cl$participant_id)
  # The clinical cohort includes all families with >=2 clinically observed
  # generations; no older stool measurement is required in this module.
  keep<-!is.na(k);m<-m[keep,,drop=FALSE];k<-k[keep]
  m$family_id<-fl_primary_family_key(m$family_id)
  stopifnot(!anyNA(m$family_id),all(m$family_id==cl$family_id[k]),
    all(as.character(m$generation)==as.character(cl$generation[k])),all(m$fatty_liver==cl$fatty_liver[k]))
  cd<-fl_cv_coldata(cfg$microbiome_file);j<-match(m$sample_id,rownames(cd));stopifnot(!anyNA(j),!anyDuplicated(rownames(cd)))
  m$abx_recent<-as.numeric(cd$ANTBV18[j]);m$smoking<-as.numeric(cd$smoke18[j]);m$physical_activity<-as.numeric(cd$PAI18[j])
  stopifnot(all(is.na(m$abx_recent)|m$abx_recent%in%0:1),all(is.na(m$smoking)|m$smoking%in%0:1))
  m$physical_activity[is.na(m$age)|m$age<18]<-NA_real_
  month<-as.integer(format(cd$ULOSTE_PVM18[j],'%m'));m$season_sin<-sin(2*pi*(month-1)/12);m$season_cos<-cos(2*pi*(month-1)/12)
  m$alcohol<-m$alkodpd18;m$diet_score<-m$dietscore_g18
  ns<-fl_cv_coldata('E:/YFS/workspace/microbiome/tse_insulin.rds');stopifnot(!anyDuplicated(ns$SampleID))
  nj<-match(m$sample_id,as.character(ns$SampleID));kcal<-as.numeric(ns$ENERCK18[nj]);stopifnot(all(is.na(kcal)|kcal>0));m$energy_log<-log10(kcal)
  for(v in c('generation','sex','prep_batch','seq_pool','family_id'))m[[v]]<-factor(m[[v]])
  m$generation<-factor(m$generation,levels=c('G0','G1','G2'))
  counts<-d$full$counts[,match(m$sample_id,colnames(d$full$counts)),drop=FALSE]
  stopifnot(all(counts>=0),all(is.finite(counts)),all(counts==floor(counts)),all(colSums(counts)>0))
  relative<-sweep(counts,2,colSums(counts),'/');stopifnot(max(abs(colSums(relative)-1))<1e-10)
  # Verify the overlapping adopted participants, but do not inherit its
  # two-microbiome-generation restriction or its preselected candidate union.
  a<-match(x$metadata$sample_id,m$sample_id);stopifnot(!anyNA(a))
  for(v in unique(c(v3_B2,'fatty_liver'))){
    aa<-m[[v]][a];bb<-x$metadata[[v]]
    if(is.factor(aa)||is.factor(bb))stopifnot(identical(as.character(aa),as.character(bb))) else
      stopifnot(isTRUE(all.equal(as.numeric(aa),as.numeric(bb),tolerance=1e-10)))
  }
  prepared<-list(metadata=m,clinical=cl,relative=relative,features=d$full$features,
    prior_microbiome_ids=as.character(x$metadata$sample_id),protocol=fp_protocol)
  v3_table(data.frame(cohort=c('Adopted clinical two-generation cohort','All stool profiles linked to clinical two-generation families',
    'Prior two-microbiome-generation baseline','Linked younger stool profiles added before comparison-specific restrictions'),
    N=c(nrow(cl),nrow(m),nrow(x$metadata),sum(!m$sample_id%in%x$metadata$sample_id))), '00_Source_population_EN',out)
  fl_write_csv(data.frame(path=c(src$files[1],path,cfg$microbiome_file,'E:/YFS/workspace/microbiome/tse_insulin.rds'),
    md5=unname(tools::md5sum(c(src$files[1],path,cfg$microbiome_file,'E:/YFS/workspace/microbiome/tse_insulin.rds')))),file.path(out,'audit/input_fingerprints.csv'))
  prepared
}

fp_comparisons <- function(x,out) {
  sets<-flow<-desc<-missing<-list()
  for(g in v3_pairs){key<-paste(g,collapse='-');older<-x$clinical[x$clinical$generation==g[1],,drop=FALSE]
    blocks<-split(seq_len(nrow(older)),as.character(older$family_id))
    senior<-ov_bind(lapply(names(blocks),function(f){z<-older[blocks[[f]],,drop=FALSE];data.frame(family_id=f,
      n_older_observed=nrow(z),n_older_cases=sum(z$fatty_liver),older_any_case=as.integer(any(z$fatty_liver==1)))}))
    young<-droplevels(x$metadata[x$metadata$generation==g[2],,drop=FALSE]);k<-match(as.character(young$family_id),senior$family_id)
    linked<-!is.na(k);m<-young[linked,,drop=FALSE];k<-k[linked]
    for(v in c('n_older_observed','n_older_cases','older_any_case'))m[[v]]<-senior[[v]][k]
    m$comparison<-key;stopifnot(!anyDuplicated(m$sample_id),all(m$n_older_observed>=1),all(m$older_any_case==as.integer(m$n_older_cases>0)))
    healthy<-droplevels(m[m$fatty_liver==0,,drop=FALSE]);core<-v3_cc(healthy,c(v3_P2,'n_older_observed'))
    adult<-!is.na(healthy$age)&healthy$age>=18&v3_cc(healthy,c(setdiff(v3_B2,'generation'),'n_older_observed'))
    sets[[key]]<-list(eligible=m,unaffected=healthy,B1=droplevels(healthy[core,,drop=FALSE]),
      Adult_B1_same_sample=droplevels(healthy[adult,,drop=FALSE]),B2=droplevels(healthy[adult,,drop=FALSE]))
    flow[[key]]<-data.frame(comparison=key,stage=c('Younger stool profiles in clinically eligible families',
      'Has observed clinical older-generation member','Younger has personal steatosis: excluded','Younger without personal steatosis',
      'All-age B1 complete','Adult B1/B2 common complete'),
      N=c(nrow(young),nrow(m),sum(m$fatty_liver==1),nrow(healthy),sum(core),sum(adult)))
    for(model in c('B1','Adult_B1_same_sample','B2'))for(e in 0:1){z<-sets[[key]][[model]];z<-z[z$older_any_case==e,,drop=FALSE]
      ages<-quantile(z$age,c(.25,.5,.75),na.rm=TRUE,names=FALSE)
      desc[[paste(key,model,e)]]<-data.frame(comparison=key,model=model,older_any_case=e,
        exposure=if(e)'At least one observed older case'else'No observed older case',N=nrow(z),families=length(unique(z$family_id)),
        minors=sum(z$age<18),female=sum(z$sex=='Female'),age_median=ages[2],age_Q1=ages[1],age_Q3=ages[3],
        mean_observed_older=mean(z$n_older_observed),min_observed_older=min(z$n_older_observed),max_observed_older=max(z$n_older_observed),
        added_beyond_old_microbiome_scope=sum(!z$sample_id%in%x$prior_microbiome_ids))
    }
    for(v in unique(c(v3_P2,setdiff(v3_B2,'generation'))))missing[[paste(key,v)]]<-data.frame(comparison=key,variable=v,
      unaffected_N=nrow(healthy),missing=sum(is.na(healthy[[v]])),adult_N=sum(healthy$age>=18,na.rm=TRUE),
      adult_missing=sum(is.na(healthy[[v]])&healthy$age>=18,na.rm=TRUE))
    saveRDS(m,file.path(out,'private',paste0('Linked_younger_',key,'.rds')))
  }
  v3_table(ov_bind(flow),'01_Participant_flow_EN',out);v3_table(ov_bind(desc),'02_Analysis_groups_EN',out)
  v3_table(ov_bind(missing),'03_Covariate_missingness_EN',out)
  list(sets=sets,flow=ov_bind(flow),groups=ov_bind(desc))
}

fp_design <- function(x,m,model) {
  m<-droplevels(m[order(m$family_id,m$sample_id),,drop=FALSE]);stopifnot(all(m$fatty_liver==0),!anyDuplicated(m$sample_id))
  fields<-c(if(model=='B2')setdiff(v3_B2,'generation')else v3_P2,'n_older_observed')
  stopifnot(nrow(m)>20,length(unique(m$older_any_case))==2,!anyNA(m[,fields,drop=FALSE]))
  nz<-v3_nuisance(m,fields);m<-nz$metadata;Z<-nz$Z[,colnames(nz$Z)!='(Intercept)',drop=FALSE]
  centers<-colMeans(Z);scales<-apply(Z,2,sd);stopifnot(all(is.finite(scales)&scales>0))
  map<-data.frame(column=sprintf('nuis_%03d',seq_len(ncol(Z))),original=colnames(Z),center=centers,scale=scales)
  Z<-sweep(sweep(Z,2,centers,'-'),2,scales,'/');colnames(Z)<-map$column
  design<-data.frame(family_id=factor(m$family_id),older_case=factor(m$older_any_case,levels=0:1,labels=c('No','Yes')),
    Z,row.names=m$sample_id,check.names=FALSE)
  mm<-model.matrix(reformulate(c('older_case',map$column)),design);stopifnot(qr(mm)$rank==ncol(mm),nrow(mm)>ncol(mm)+5)
  A<-t(x$relative[,match(m$sample_id,colnames(x$relative)),drop=FALSE]);pv<-colMeans(A>0)
  random<-any(table(m$family_id)>1L)
  list(metadata=m,design=design,relative=A,mapping=map,nuisance=map$column,random=random,dropped=nz$dropped,
    formula=paste('LOG abundance ~ older_case +',paste(map$column,collapse=' + '),if(random)'+ (1 | family_id)'else''),
    filter=data.frame(feature_id=colnames(A),prevalence=pv,retained=pv>=.10))
}

fp_extract <- function(ans,ids,d,expected,model,key) {
  effects<-diagnostics<-list()
  for(id in ids){fit<-ans$fits[[id]];messages<-character();err<-'';singular<-NA;icc<-hmin<-NA_real_;opt<-0L;cs<-NULL
    if(inherits(fit,'merMod')){
      messages<-as.character(fit@optinfo$conv$lme4$messages);opt<-unlist(fit@optinfo$conv$opt)
      H<-fit@optinfo$derivs$Hessian;if(is.matrix(H)&&all(is.finite(H)))hmin<-min(eigen(H,symmetric=TRUE,only.values=TRUE)$values)
      singular<-lme4::isSingular(fit,tol=1e-4)
      data<-model.frame(fit);formula<-formula(fit);na.action<-stats::na.fail
      ft<-tryCatch(withCallingHandlers(lmerTest::as_lmerModLmerTest(fit),warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);NULL})
      if(!is.null(ft))cs<-coef(summary(ft))
      if(!is.null(cs))df<-cs['older_caseYes','df']else df<-NA_real_
      vv<-as.data.frame(lme4::VarCorr(fit));vf<-sum(vv$vcov[vv$grp!='Residual']);icc<-vf/(vf+sigma(fit)^2)
    }else if(inherits(fit,'glm')||inherits(fit,'lm')){cs<-coef(summary(fit));df<-df.residual(fit)}else{err<-'Official MaAsLin2 did not return a fitted model';df<-NA_real_}
    b<-se<-p<-lo<-hi<-NA_real_;finite<-FALSE;matched<-FALSE;response_error<-NA_real_
    if(!is.null(cs)&&'older_caseYes'%in%rownames(cs)){
      mf<-model.frame(fit);at<-match(rownames(mf),rownames(expected));stopifnot(!anyNA(at))
      response_error<-max(abs(as.numeric(model.response(mf))-expected[at,id]));stopifnot(response_error<1e-10)
      b<-cs['older_caseYes','Estimate'];se<-cs['older_caseYes','Std. Error'];p<-cs['older_caseYes','Pr(>|t|)']
      finite<-all(is.finite(c(b,se,p,df)))&&se>0&&df>0&&p>=0&&p<=1
      if(finite){critical<-qt(.975,df);lo<-b-critical*se;hi<-b+critical*se}
      rr<-ans$results;rr<-rr[rr$feature==id&rr$metadata=='older_case'&as.character(rr$value)=='Yes',,drop=FALSE]
      if(nrow(rr)==1L&&finite)matched<-isTRUE(all.equal(unname(b),rr$coef[1],tolerance=1e-7))&&isTRUE(all.equal(unname(p),rr$pval[1],tolerance=1e-7))
    }
    conv<-messages[!grepl('boundary.*singular|singular fit',messages,ignore.case=TRUE)]
    bad<-nzchar(err)||length(conv)>0||any(opt!=0)||(!is.na(hmin)&&hmin< -1e-6)
    state<-if(!finite)'NOT_ESTIMABLE'else if(!matched)'OFFICIAL_RESULT_MISMATCH'else if(bad)'UNRELIABLE_CONVERGENCE'else if(isTRUE(singular))'ESTIMATED_SINGULAR'else'ESTIMATED'
    rv<-if(!is.null(cs))residuals(fit)else NA_real_;sr<-sd(rv);skew<-if(is.finite(sr)&&sr>0)mean((rv-mean(rv))^3)/sr^3 else NA_real_
    effects[[id]]<-data.frame(comparison=key,model=model,feature_id=id,beta_log2=b,SE=se,df=df,CI_low=lo,CI_high=hi,P_raw=p,
      analysis_valid=state%in%c('ESTIMATED','ESTIMATED_SINGULAR'),status=state,singular=singular)
    diagnostics[[id]]<-data.frame(comparison=key,model=model,feature_id=id,status=state,random_intercept=d$random,
      family_ICC=icc,residual_skewness=skew,minimum_Hessian_eigenvalue=hmin,optimizer_code=paste(opt,collapse=';'),
      official_beta_P_match=matched,response_max_error=response_error,warnings=paste(unique(messages),collapse='; '),error=err)
  }
  list(effects=ov_bind(effects),diagnostics=ov_bind(diagnostics))
}

fp_fit <- function(x,m,model,key,cfg,out,cache) {
  d<-fp_design(x,m,model);ids<-d$filter$feature_id[d$filter$retained]
  prefix<-paste(key,model,sep='_');stopifnot(length(ids)>0)
  fl_write_csv(d$mapping,file.path(out,'audit',paste0(prefix,'_design.csv')))
  fl_write_csv(d$filter,file.path(out,'audit',paste0(prefix,'_filter.csv')))
  saveRDS(list(metadata=d$metadata,design=d$design,filter=d$filter,formula=d$formula,random=d$random,dropped=d$dropped),file.path(out,'private',paste0(prefix,'_design.rds')))
  sig<-digest::digest(list(protocol=fp_protocol,A=d$relative[,ids,drop=FALSE],design=d$design,random=d$random,
    Maaslin2=as.character(packageVersion('Maaslin2')),extract=body(fp_extract)),algo='sha256')
  pieces<-list();batches<-split(ids,ceiling(seq_along(ids)/25L));cache<-file.path(cache,prefix,substr(sig,1,20))
  tf_log(key,' ',model,': ',nrow(m),' unaffected younger members, ',sum(m$older_any_case),' exposed; ',length(ids),' SGBs; ',if(d$random)'family LMM'else'LM: one younger per family')
  for(i in seq_along(batches)){keys<-batches[[i]]
    pieces[[i]]<-tf_checkpoint(file.path(cache,sprintf('batch_%03d.rds',i)),paste(sig,paste(keys,collapse=';')),function(){
      tf_log(prefix,': official MaAsLin2 batch ',i,'/',length(batches))
      A<-d$relative[,keys,drop=FALSE];expected<-vapply(keys,function(k)get('LOG',asNamespace('Maaslin2'))(A[,k]),numeric(nrow(A)))
      dimnames(expected)<-dimnames(A);folder<-file.path(cache,sprintf('official_%03d',i));dir.create(folder,recursive=TRUE,showWarnings=FALSE)
      console<-file(file.path(folder,'console.log'),'wt');sink(console)
      ans<-tryCatch(Maaslin2::Maaslin2(input_data=as.data.frame(A),input_metadata=d$design,output=file.path(folder,'fit'),
        fixed_effects=c('older_case',d$nuisance),random_effects=if(d$random)'family_id'else NULL,
        normalization='NONE',transform='LOG',standardize=FALSE,min_prevalence=0,min_abundance=0,min_variance=0,
        analysis_method='LM',cores=1,correction='BH',plot_heatmap=FALSE,plot_scatter=FALSE,save_models=TRUE),finally={sink();close(console)})
      z<-fp_extract(ans,keys,d,expected,model,key)
      saveRDS(list(feature_ids=keys,package_results=ans$results),file.path(folder,'official_results.rds'))
      z
    })
  }
  e<-ov_bind(lapply(pieces,`[[`,'effects'));dg<-ov_bind(lapply(pieces,`[[`,'diagnostics'))
  fl_write_csv(dg,file.path(out,'audit',paste0(prefix,'_diagnostics.csv')))
  j<-match(d$filter$feature_id,e$feature_id);t<-e[j,,drop=FALSE];t$feature_id<-d$filter$feature_id;t$comparison<-key;t$model<-model
  t$status[is.na(j)]<-'FILTERED_PREVALENCE_LT_10_PERCENT';t$analysis_valid[is.na(j)]<-FALSE
  t$eligible<-d$filter$retained;t$prevalence<-d$filter$prevalence
  t$N<-nrow(m);t$exposed_N<-sum(m$older_any_case);t$reference_N<-sum(m$older_any_case==0);t$families<-length(unique(m$family_id));t$minors<-sum(m$age<18)
  t$random_intercept<-d$random;t$SGB_id<-sub('^t__','',t$feature_id);t$species_label<-ov_species_label(x$features$taxonomy[match(t$feature_id,x$features$feature_id)])
  t$abundance_ratio<-2^t$beta_log2;t$ratio_CI_low<-2^t$CI_low;t$ratio_CI_high<-2^t$CI_high
  t$effect_unit<-'Log2 zero-replaced relative abundance: observed older case minus no observed older case, among unaffected younger members'
  rownames(t)<-NULL;t
}

fp_adjust <- function(t) {
  t$P_BH<-t$P_BH_within_comparison<-NA_real_;t$FDR_tests<-0L
  for(model in unique(t$model)){
    k<-which(t$model==model&t$eligible);good<-tf_bool(t$analysis_valid[k])&is.finite(t$P_raw[k]);p<-rep(1,length(k));p[good]<-t$P_raw[k][good]
    q<-p.adjust(p,'BH');t$P_BH[k[good]]<-q[good];t$FDR_tests[t$model==model]<-length(k)
    for(key in unique(t$comparison)){a<-which(t$model==model&t$comparison==key&t$eligible);valid<-tf_bool(t$analysis_valid[a])&is.finite(t$P_raw[a]);pp<-rep(1,length(a));pp[valid]<-t$P_raw[a][valid]
      qq<-p.adjust(pp,'BH');t$P_BH_within_comparison[a[valid]]<-qq[valid]}
  }
  t$FDR_family<-paste0(t$model,'_all_eligible_SGBs_x_three_comparisons');t$FDR_lt_005<-is.finite(t$P_BH)&t$P_BH<.05
  t$direction<-ifelse(!tf_bool(t$analysis_valid),'Not estimable',ifelse(t$beta_log2<0,'Lower with older case','Higher with older case'))
  t
}

fp_analyse <- function(x,sets,cfg,out,cache) {
  all<-list()
  for(key in names(sets$sets))for(model in c('B1','Adult_B1_same_sample','B2'))
    all[[paste(key,model)]]<-fp_fit(x,sets$sets[[key]][[model]],model,key,cfg,out,cache)
  t<-fp_adjust(ov_bind(all));v3_table(t,'04_All_SGB_family_exposure_results_EN',out)
  summary<-ov_bind(lapply(split(t,list(t$comparison,t$model),drop=TRUE),function(z)data.frame(comparison=z$comparison[1],model=z$model[1],
    N=z$N[1],exposed_N=z$exposed_N[1],reference_N=z$reference_N[1],families=z$families[1],minors=z$minors[1],
    random_intercept=z$random_intercept[1],eligible_SGBs=sum(z$eligible),valid=sum(z$analysis_valid),
    unreliable=sum(z$eligible&!z$analysis_valid),singular=sum(z$analysis_valid&z$singular,na.rm=TRUE),
    raw_P_lt_005=sum(z$analysis_valid&z$P_raw<.05,na.rm=TRUE),BH_joint_lt_005=sum(z$FDR_lt_005),
    BH_lower=sum(z$FDR_lt_005&z$beta_log2<0,na.rm=TRUE),BH_higher=sum(z$FDR_lt_005&z$beta_log2>0,na.rm=TRUE))))
  v3_table(summary,'05_Model_summary_EN',out)
  selected<-t[t$FDR_lt_005,,drop=FALSE];selected<-selected[order(selected$model,selected$P_BH,selected$P_raw),,drop=FALSE]
  v3_table(selected,'06_BH_supported_associations_EN',out)
  # Existing disease/concordance results are contextual annotations only.
  disease<-tf_read(file.path(tf_source(cfg$project_root)$path,'tables/Table6_SGB_B1_all_results_EN.csv'))
  ann<-t[t$model=='B1',,drop=FALSE];a<-match(ann$feature_id,disease$feature_id)
  ann$personal_steatosis_beta_25<-disease$beta_log2[a];ann$personal_steatosis_BH_25<-disease$P_BH[a]
  p26<-file.path(cfg$project_root,'outputs/26_Taxon_family_sharing');id26<-readLines(file.path(p26,'LATEST_RUN.txt'),1)
  c26<-tf_read(file.path(p26,'current',id26,'tables/Taxon_cross_generation_GEE_all_results_EN.csv'));c26<-c26[c26$model=='P2',]
  a<-match(paste(ann$feature_id,ann$comparison),paste(c26$feature_id,c26$comparison))
  ann$concordance_beta_26<-c26$beta_log2_per_log2[a];ann$concordance_BH_26<-c26$P_BH[a]
  ann$annotation_note<-'Same-cohort contextual evidence; not independent validation or mediation'
  v3_table(ann,'07_Context_with_existing_disease_and_concordance_results_EN',out)
  for(key in names(sets$sets))stopifnot(identical(sets$sets[[key]]$Adult_B1_same_sample$sample_id,sets$sets[[key]]$B2$sample_id))
  list(results=t,summary=summary,groups=sets$groups,flow=sets$flow)
}

fp_reports <- function(r,cfg,out) {
  methods<-c('# Methods and references','',paste('Protocol:',fp_protocol),
    'Question: among younger participants without personal ultrasound steatosis, is at least one observed older-generation steatosis case associated with their microbial abundances?',
    'Only the relevant two clinically observed generations are required. Older members need ultrasound records, not stool samples. Younger members each contribute one row per comparison. Only the primary family label is used. Missing older examinations are not disease-negative observations. Minors remain eligible.',
    'Primary exposure: any steatosis among the observed older members; the reference has no cases among those observed. Ultrasound RASVAMAKSA18 codes 2/3 indicate steatosis, code 1 indicates absence. No MASLD analysis.',
    'Official MaAsLin2 1.24.1; LM analysis_method, relative abundances, normalization NONE, transform LOG, standardize FALSE. Relative abundances are calculated from the full processed SGB matrix before filtering, not renormalized after filtering. The official LOG function replaces zeros by half the within-model minimum positive abundance of that feature, then takes log2. Candidate SGBs are not selected by earlier disease P values; retain prevalence >=10% in the actual analysis sample.',
    'B1 covariates: younger age natural spline (up to 3 df), sex, library preparation batch, sequencing pool, season sine/cosine, recent antibiotics, and observed older-member count if variable. Generation is constant within each comparison and is not entered. Continuous nuisance columns are scaled only for conditioning; exposure and microbial outcome are not SD-standardized.',
    'Use a family random intercept when families contribute repeated younger members. If every family contributes one younger member, use the official ordinary Gaussian LM. Random-effect boundary/singularity is reported, not automatically treated as optimizer failure.',
    'Adult extension: use one common complete sample age >=18 for adult B1 and B2; B2 adds smoking, alcohol, diet score, physical activity, log10 total energy and continuous BMI. No new imputation is performed.',
    'Official beta/P must match the extracted fit. LMM CI uses Satterthwaite t degrees of freedom; ordinary LM CI uses residual t degrees of freedom. Record optimizer warnings, Hessian, singularity, residual skewness, family ICC and response-order checks. Invalid estimates are counted conservatively as P=1 in multiplicity, and have no reported q.',
    'Main BH correction pools all eligible SGB-by-three-comparison tests within each model set (B1, adult B1, B2 separately). Also provide within-comparison BH as a labelled secondary source-style result, never substitute it for the joint main q. Per-computational-batch package q values are discarded.',
    '2^beta is the adjusted ratio of geometric means on the official zero-replaced relative-abundance scale. It is not an OR, absolute microbial load, transmission probability, genetic heritability or a future disease-risk estimate.',
    '', '## Published method and cohort adaptations',
    'Olivera et al. used unaffected verified first-degree relatives, CD multiplex versus simplex family burden, 16S genus relative abundances, MaAsLin2 with demographic adjustment/family clustering and BH q<0.05. Their age/sex/race/family-size adjustment is not claimed to equal our B1/B2. Our observed older steatosis exposure, primary labels, metagenomic SGBs, covariate set and joint three-comparison correction are explicitly YFS applications.',
    'Lee et al. support the healthy-offspring/parental-disease-exposure approach and >=10% prevalence; their contrast was parental diagnosis during/before versus after the perinatal period and their microbiome model adjusted offspring age/sex. We do not claim to reproduce that exposure timing.',
    'Mallick et al. provide the official LOG multivariable association framework. The existing approved cohort covariates and family grouping are reused rather than introducing an unpublished new abundance estimator. The current analysis is not a literal replication of an entire Crohn study.',
    'The deferred inflammation/clinical-marker module is not run. No new three-group disease-stage, score, prediction or mediation model is included.',
    '', '## References',
    'Olivera PA, Martinez-Lozano H, Leibovitzh H, et al. Healthy First-Degree Relatives From Multiplex Families vs Simplex Families Have Higher Subclinical Intestinal Inflammation, a Distinct Fecal Microbial Signature, and Harbor a Higher Risk of Developing Crohn\'s Disease. Gastroenterology. 2025;168(1):99-110.e2. doi:10.1053/j.gastro.2024.08.031. https://doi.org/10.1053/j.gastro.2024.08.031',
    'Lee SH, Bushra M, Qiu L, et al. Early Life Exposure to Parental Crohn\'s Disease Is Associated With Offspring\'s Gut Microbiome, Gut Permeability, and Increased Risk of Future Crohn\'s Disease. Gastroenterology. 2025;168(2):385-388.e3. doi:10.1053/j.gastro.2024.09.033. https://doi.org/10.1053/j.gastro.2024.09.033',
    'Mallick H, Rahnavard A, McIver LJ, et al. Multivariable association discovery in population-scale meta-omics studies. PLOS Computational Biology. 2021;17(11):e1009442. doi:10.1371/journal.pcbi.1009442. https://doi.org/10.1371/journal.pcbi.1009442')
  writeLines(methods,file.path(out,'reports/Methods_and_references_EN.md'),useBytes=TRUE)
  s<-r$summary[r$summary$model=='B1',];s<-s[match(c('G0-G1','G1-G2','G0-G2'),s$comparison),]
  lines<-c('# 27号：年长成员脂肪肝与未患病年轻成员菌群','',
    '仅运行第一项。年轻成员均无超声脂肪肝；暴露为同一主家系已观察年长成员至少一人有脂肪肝。年长成员不要求有菌群样本。',
    '| 比较 | 未患病年轻成员 | 暴露组 | 参照组 | 家系 | 合格SGB | 主BH q<0.05 | 较低 | 较高 |',
    '|---|---:|---:|---:|---:|---:|---:|---:|---:|')
  for(i in seq_len(nrow(s)))lines<-c(lines,sprintf('| %s | %d | %d | %d | %d | %d | %d | %d | %d |',s$comparison[i],s$N[i],s$exposed_N[i],s$reference_N[i],s$families[i],s$eligible_SGBs[i],s$BH_joint_lt_005[i],s$BH_lower[i],s$BH_higher[i]))
  if(sum(s$BH_joint_lt_005)==0)lines<-c(lines,'',
    '三组全龄主分析均未检出统一BH校正后显著的菌。目前没有充分统计证据表明：在自身无脂肪肝的年轻成员中，观察到的年长一代脂肪肝与特定菌丰度有关。这不等于证明不存在关联，也不否定此前个人脂肪肝与菌群的关联或同一家系跨代丰度相关。',
    paste0('主模型原始P<0.05项数按G0-G1/G1-G2/G0-G2为：',paste(s$raw_P_lt_005,collapse='/'),'；最小统一BH q为',format(min(r$results$P_BH[r$results$model=='B1'],na.rm=TRUE),digits=4),'。原始P支持项不能表述为多重校正后的发现。'))
  lines<-c(lines,'','## 成人同样本补充结果','',
    '| 比较 | 模型 | 人数 | 暴露组 | 参照组 | 合格SGB | 主BH q<0.05 |',
    '|---|---|---:|---:|---:|---:|---:|')
  a<-r$summary[r$summary$model!='B1',]
  for(i in seq_len(nrow(a)))lines<-c(lines,sprintf('| %s | %s | %d | %d | %d | %d | %d |',a$comparison[i],a$model[i],a$N[i],a$exposed_N[i],a$reference_N[i],a$eligible_SGBs[i],a$BH_joint_lt_005[i]))
  lines<-c(lines,'','主调整：年轻成员年龄样条、性别、建库批次、测序池、季节、近期抗生素，以及非恒定的年长成员观察人数。成人扩展在同一完整样本再加入吸烟、饮酒、饮食评分、体力活动、总能量、BMI。',
    '所有合格菌均检验。主BH按每套模型的全部菌×三个代际比较统一校正；各比较内BH另列，仅供透明比较。',
    '丰度比为零值替换后的相对丰度尺度，不是绝对菌量或OR。主家系未核实直接亲子关系；不据此推断菌群传播、发病前变化或中介作用。',
    '原文采用已核实亲属、克罗恩病负担和属水平数据；这里沿用官方MaAsLin2分析方式，疾病、主家系标签、SGB及调整集为本队列适配。完整来源见Methods_and_references_EN.md。')
  writeLines(lines,file.path(out,'reports/Results_summary_zh.md'),useBytes=TRUE)
  capture.output(sessionInfo(),file=file.path(out,'audit/sessionInfo.txt'))
}

fp_publish <- function(out,cfg) {
  if(file.exists(file.path(cfg$project_root,'config/results_viewer.json'))){
    gallery<-new.env(parent=baseenv());sys.source(file.path(cfg$project_root,'R/result_gallery.R'),envir=gallery)
    return(gallery$fl_publish_organized_results(cfg))
  }
  viewer<-file.path(cfg$project_root,'Results_表格与图片/Family_steatosis_microbiome_20261027');dir.create(viewer,recursive=TRUE,showWarnings=FALSE)
  files<-unlist(lapply(c('tables','figures','reports'),function(p)list.files(file.path(out,p),recursive=TRUE,full.names=TRUE)))
  rel<-substring(files,nchar(out)+2L);dest<-file.path(viewer,rel)
  for(i in seq_along(files)){dir.create(dirname(dest[i]),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(files[i],dest[i],overwrite=TRUE))}
  stopifnot(identical(unname(tools::md5sum(files)),unname(tools::md5sum(dest))))
  fl_write_csv(data.frame(source=files,viewing_copy=dest,md5=unname(tools::md5sum(files))),file.path(out,'audit/viewer_copy_manifest.csv'))
  writeLines(c('<!doctype html><meta charset="utf-8"><style>body{font:15px Arial;margin:28px;line-height:1.6}img{max-width:100%;width:1100px}</style>',
    '<h1>Observed older steatosis and the microbiome of unaffected younger members</h1>',
    '<p>Official MaAsLin2; primary-family labels; ultrasound steatosis; relevant two generations.</p>',
    '<ul>',paste0('<li><a href="',gsub('\\','/',rel,fixed=TRUE),'">',basename(rel),'</a></li>'),'</ul>',
    paste0('<img src="',gsub('\\','/',rel[grepl('\\.png$',rel)],fixed=TRUE),'">')),file.path(viewer,'index.html'),useBytes=TRUE)
  invisible(viewer)
}

fp_main <- function(args=character(),root=Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')) {
  if(getRversion()<'4.5'){
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts/27_family_steatosis_microbiome.R')),vapply(args,shQuote,character(1))))
    if(status!=0)stop('Project-R family-exposure analysis failed');return(invisible(NULL))
  }
  .libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),file.path(R.home(),'library'),.libPaths()))
  for(f in c('config.R','R/common.R','R/entry_helpers.R','R/confounders.R','R/nature_export_approved.R','R/reanalysis_numerical_helpers.R',
    'R/reanalysis_20261024.R','R/taxa_followup.R','R/taxa_followup_figures.R','R/family_steatosis_microbiome_figures.R',
    'R/family_steatosis_microbiome_diagnostics.R'))sys.source(file.path(root,f),envir=.GlobalEnv)
  cfg<-fl_config(root);cfg$output_module<-unname(cfg$output_modules['family_steatosis_microbiome']);cfg$run_id<-unname(cfg$run_ids['family_steatosis_microbiome']);cfg$run_id_explicit<-FALSE
  deps<-c('Maaslin2','lme4','lmerTest','digest','ggplot2','patchwork','svglite','jsonlite');stopifnot(all(vapply(deps,requireNamespace,logical(1),quietly=TRUE)),as.character(packageVersion('Maaslin2'))=='1.24.1')
  if(!identical(as.character(Matrix::Matrix.Version()$abi),as.character(lme4:::.Matrix.abi.build.version)))stop('Matrix ABI mismatch')
  src<-tf_source(root)
  if('--check'%in%args){cat('PASS: source, package and R/Matrix ABI checks. No model fits.\n');return(invisible(NULL))}
  if(any(!args%in%c('--run','--prepare-only','--figures-only')))stop('Supported: --check, --run, --prepare-only, --figures-only')
  module<-file.path(root,'outputs',cfg$output_module);cache<-file.path(module,'_checkpoints/private');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  Sys.setenv(PYTHONPATH=paste(c(file.path(root,'.figure_qa_python'),Sys.getenv('PYTHONPATH')),collapse=.Platform$path.sep),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
  if('--figures-only'%in%args){id<-readLines(file.path(module,'LATEST_RUN.txt'),1);out<-file.path(module,'current',id)
    r<-readRDS(file.path(out,'private/analysis_result.rds'));fp_validate(r,cfg,out);fp_diagnostics(r,cfg,out);fp_figures(r,cfg,out);fp_publish(out,cfg);return(invisible(out))}
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(p in c('tables','figures','reports','audit','private'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
  writeLines(c('RUNNING',fp_protocol),file.path(out,'run_status.txt'));saveRDS(cfg,file.path(out,'private/run_config.rds'))
  tf_log('START module 27; output ',out)
  tryCatch({
    sig<-digest::digest(list(protocol=fp_protocol,source_MD5=src$md5[1],prepare=body(fp_prepare)),algo='sha256')
    x<-tf_checkpoint(file.path(cache,'expanded_population.rds'),sig,function()fp_prepare(src,cfg,out))
    # Reproduce source audit even when resuming preparation from checkpoint.
    v3_table(data.frame(cohort=c('Adopted clinical two-generation cohort','All stool profiles linked to clinical two-generation families',
      'Prior two-microbiome-generation baseline','Linked stool profiles added before comparison-specific restrictions'),
      N=c(nrow(x$clinical),nrow(x$metadata),length(x$prior_microbiome_ids),sum(!x$metadata$sample_id%in%x$prior_microbiome_ids))),
      '00_Source_population_EN',out)
    source_paths<-c(src$files[1],file.path(root,'data/processed/processed_microbiome_v1/microbiome_processed.rds'),
      cfg$microbiome_file,'E:/YFS/workspace/microbiome/tse_insulin.rds')
    fl_write_csv(data.frame(path=source_paths,md5=unname(tools::md5sum(source_paths))),file.path(out,'audit/input_fingerprints.csv'))
    saveRDS(x,file.path(out,'private/analysis_population.rds'));sets<-fp_comparisons(x,out)
    fl_write_csv(sets$groups,file.path(out,'audit/pre_fit_group_counts.csv'))
    if('--prepare-only'%in%args){writeLines(c('PREPARED_ONLY',fp_protocol),file.path(out,'run_status.txt'));tf_log('PREPARE ONLY: ',out);return(invisible(out))}
    r<-fp_analyse(x,sets,cfg,out,cache);saveRDS(r,file.path(out,'private/analysis_result.rds'))
    fp_reports(r,cfg,out);fp_validate(r,cfg,out);fp_diagnostics(r,cfg,out);fp_figures(r,cfg,out)
    stopifnot(identical(unname(src$md5),unname(tools::md5sum(src$files))))
    writeLines(c('COMPLETED',fp_protocol),file.path(out,'run_status.txt'));target<-fl_output_complete_run(ctx)
    writeLines(cfg$run_id,file.path(module,'LATEST_RUN.txt'));viewer<-fp_publish(target,cfg)
    tf_log('COMPLETED: ',target);tf_log('Viewer: ',viewer);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
