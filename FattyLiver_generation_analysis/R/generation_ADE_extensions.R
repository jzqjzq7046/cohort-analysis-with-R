# Inert definitions: user-authorized A/D/E extension. Original data are read-only.
ade_protocol <- 'generation_function_strain_v1_20261101'
ade_root <- function() Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
ade_log <- function(...) {
  line<-paste(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),paste0(...,collapse=''))
  cat(line,'\n');flush.console();p<-getOption('fatty_liver.ade_log')
  if(is.character(p)&&length(p)==1L)cat(line,'\n',file=p,append=TRUE)
}
ade_setup <- function(root) {
  if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
  .libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),file.path(R.home(),'library')))
  for(f in c('config.R','R/common.R','R/entry_helpers.R','R/nature_export_approved.R',
    'R/reanalysis_numerical_helpers.R','R/reanalysis_20261024.R','R/taxa_followup.R','R/taxa_followup_figures.R'))
    sys.source(file.path(root,f),envir=.GlobalEnv)
  pk<-c('Maaslin2','lme4','lmerTest','ggplot2','patchwork','svglite','jsonlite','digest','SummarizedExperiment','TreeSummarizedExperiment')
  if(!all(vapply(pk,requireNamespace,logical(1),quietly=TRUE)))stop('Missing project dependency')
  stopifnot(as.character(packageVersion('Maaslin2'))=='1.24.1',
    identical(as.character(Matrix::Matrix.Version()$abi),as.character(lme4:::.Matrix.abi.build.version)))
  Sys.setenv(PYTHONPATH=paste(c(file.path(root,'.figure_qa_python'),Sys.getenv('PYTHONPATH')),collapse=.Platform$path.sep),
    OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
  fl_config(root)
}
ade_table <- function(x,name,out) {
  if(!ncol(x))x<-data.frame(status=character())
  v3_table(x,name,out)
}
ade_labels <- function(f) {
  f$SGB_id<-sub('^t__','',f$feature_id)
  # Same conservative binomial rule as the adopted module 40.
  rows<-lapply(seq_len(nrow(f)),function(i){
    tx<-f$taxonomy[i];if(is.na(tx))tx<-''
    tok<-unlist(strsplit(tx,'[|,]'));ss<-sort(unique(sub('^s__','',tok[grepl('^s__',tok)])))
    ss<-ss[nzchar(ss)];readable<-gsub('_',' ',ss,fixed=TRUE)
    named<-length(ss)==1L&&grepl('^(Candidatus )?[A-Z][a-z]+ [a-z][a-z-]+$',readable)&&
      !grepl('\\b(bacterium|archaeon|unclassified|uncultured|unidentified|unknown|sp|spp)\\b',readable,ignore.case=TRUE,perl=TRUE)
    nm<-if(named)readable else NA_character_
    data.frame(feature_id=f$feature_id[i],feature_label=if(named)paste0(nm,' (',f$SGB_id[i],')')else f$SGB_id[i],
      scientific_name=nm,annotation_candidates=paste(readable,collapse='; '),
      annotation_status=if(length(ss)>1)'Conflicting'else if(named)'Named'else'Unresolved',taxonomy=tx)
  });ov_bind(rows)
}
ade_fingerprint <- function(paths) data.frame(path=normalizePath(paths,winslash='/',mustWork=TRUE),md5=unname(tools::md5sum(paths)))
ade_protect <- function(h) stopifnot(all(file.exists(h$path)),identical(unname(tools::md5sum(h$path)),h$md5))

ade_import_paths <- function(cfg,x,out) {
  p<-'E:/YFS/documents-export-2026-08-31/tse_humann_pathabundance_h5se'
  files<-file.path(p,c('se.rds','assays.h5'));stopifnot(all(file.exists(files)))
  protected<-ade_fingerprint(files)
  # SE row/column identities can be read without HDF5Array. Numeric import is
  # delegated to the existing local h5py; no metadata is inferred from position.
  z<-readRDS(files[1]);stopifnot(inherits(z,'SummarizedExperiment'))
  rn<-rownames(z);cn<-colnames(z);rd<-as.data.frame(SummarizedExperiment::rowData(z))
  stopifnot(length(rn)==nrow(z),length(cn)==ncol(z),!anyDuplicated(rn),!anyDuplicated(cn))
  ii<-which(!grepl('|',rn,fixed=TRUE));stopifnot(length(ii)>1L)
  desc<-list(dataset='assay001',n_samples=length(cn),n_source_rows=length(rn),zero_based_row_indices=ii-1L)
  mf<-file.path(out,'private/humann_row_indices.json');bin<-file.path(out,'private/humann_community.bin')
  jsonlite::write_json(desc,mf,auto_unbox=TRUE,pretty=TRUE)
  status<-system2(cfg$python,c(shQuote(file.path(cfg$project_root,'tools/read_humann_h5.py')),
    shQuote(files[2]),shQuote(mf),shQuote(bin)))
  if(status!=0L)stop('Existing HUMAnN HDF5 numeric import failed')
  con<-file(bin,'rb');on.exit(close(con),add=TRUE)
  a<-matrix(readBin(con,'double',n=length(ii)*length(cn),size=8,endian='little'),nrow=length(ii),
    dimnames=list(rn[ii],cn));stopifnot(length(a)==length(ii)*length(cn),all(is.finite(a)),all(a>=0))
  j<-match(x$metadata$sample_id,cn);stopifnot(!anyNA(j))
  cd<-as.data.frame(SummarizedExperiment::colData(z))
  # Cross-check sample mapping with the preserved participant identifier.
  v<-if('tutknom18'%in%names(cd))'tutknom18'else if('TUTKNOM18'%in%names(cd))'TUTKNOM18'else NULL
  if(!is.null(v))stopifnot(all(as.character(cd[[v]][j])==as.character(x$metadata$participant_id)))
  a<-a[,j,drop=FALSE];denom<-colSums(a);stopifnot(all(is.finite(denom)),all(denom>0))
  fid<-if('Gene_id'%in%names(rd))as.character(rd$Gene_id[ii])else sub(':.*','',rn[ii])
  stopifnot(!anyDuplicated(fid),all(nzchar(fid)))
  special<-toupper(fid)%in%c('UNMAPPED','UNINTEGRATED','UNGROUPED','READS_UNMAPPED')
  biological<-!special;ra<-sweep(a,2,denom,'/');rownames(ra)<-fid
  features<-data.frame(feature_id=fid,feature_label=fid,
    pathway_name=if('Gene'%in%names(rd))as.character(rd$Gene[ii])else sub('^[^:]+: ?','',rn[ii]),
    original_label=rn[ii],special=special,stringsAsFactors=FALSE)
  coverage<-data.frame(source_rows=length(rn),community_rows=length(ii),biological_pathways=sum(biological),
    special_rows=sum(special),source_samples=length(cn),matched_baseline_samples=length(j),
    identity_crosscheck=if(is.null(v))'Sample names only'else paste0('Sample names and ',v),
    denominator='Sum of all supplied unstratified community rows before filtering; special rows included when supplied',
    source_sum_min=min(denom),source_sum_median=median(denom),source_sum_max=max(denom),
    baseline_RA_sum_max_error=max(abs(colSums(ra)-1)))
  ade_table(features,'Pathway_annotation_and_filter_scope_EN',out);ade_table(coverage,'HUMAnN_data_availability_EN',out)
  saveRDS(list(metadata=x$metadata,relative=ra,features=features,candidates=fid[biological],denominator=denom,
    source_files=files,protected=protected),file.path(out,'private/functional_input.rds'))
  x$relative<-ra;x$candidates<-fid[biological];x$features<-features
  list(x=x,protected=protected,coverage=coverage)
}

ade_extract <- function(ans,keys,expected,design,model,mode,feature_ids) {
  effects<-joint<-diagnostics<-list()
  terms<-if(mode=='generation')paste0('case_G',0:2)else'case_overall'
  for(k in seq_along(keys)){
    key<-keys[k];id<-feature_ids[k];fit<-ans$fits[[key]];warn<-character();err<-'';ft<-cs<-NULL
    hmin<-NA_real_;singular<-NA;opt<-0;icc<-skew<-NA_real_
    if(inherits(fit,'merMod')){
      mf<-stats::model.frame(fit);at<-match(rownames(mf),rownames(expected))
      if(anyNA(at)||max(abs(stats::model.response(mf)-expected[at,key]))>1e-10)stop('Response/order mismatch: ',id)
      warn<-as.character(fit@optinfo$conv$lme4$messages);opt<-unlist(fit@optinfo$conv$opt)
      H<-fit@optinfo$derivs$Hessian
      if(is.matrix(H)&&all(is.finite(H)))hmin<-min(eigen(H,symmetric=TRUE,only.values=TRUE)$values)
      singular<-lme4::isSingular(fit,tol=1e-4)
      data<-mf;formula<-stats::formula(fit);na.action<-stats::na.fail
      ft<-tryCatch(withCallingHandlers(if(inherits(fit,'lmerModLmerTest'))fit else lmerTest::as_lmerModLmerTest(fit),
        warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);NULL})
      if(!is.null(ft))cs<-coef(summary(ft))
      vv<-as.data.frame(lme4::VarCorr(fit));vf<-sum(vv$vcov[vv$grp!='Residual']);icc<-vf/(vf+sigma(fit)^2)
      rv<-residuals(fit);if(sd(rv)>0)skew<-mean((rv-mean(rv))^3)/sd(rv)^3
    }else err<-'Official package did not return a mixed model'
    warn<-unique(warn[nzchar(warn)]);cw<-warn[!grepl('boundary.*singular|singular fit',warn,ignore.case=TRUE)]
    bad<-nzchar(err)||length(cw)>0||any(opt!=0)||(!is.na(hmin)&&hmin< -1e-6)
    state<-if(bad)'UNRELIABLE_CONVERGENCE'else if(isTRUE(singular))'ESTIMATED_SINGULAR'else'ESTIMATED'
    official_match<-TRUE;finite_all<-TRUE
    for(term in terms){
      b<-se<-df<-p<-lo<-hi<-NA_real_;ok<-FALSE;check<-FALSE
      if(!is.null(cs)&&term%in%rownames(cs)){
        b<-cs[term,'Estimate'];se<-cs[term,'Std. Error'];df<-cs[term,'df'];p<-cs[term,'Pr(>|t|)']
        ok<-all(is.finite(c(b,se,df,p)))&&se>0&&df>0&&p>=0&&p<=1
        if(ok){lo<-b-qt(.975,df)*se;hi<-b+qt(.975,df)*se}
        rr<-ans$results;rr<-rr[rr$feature==key&rr$metadata==term,,drop=FALSE]
        check<-nrow(rr)==1&&ok&&abs(b-rr$coef[1])<1e-7&&abs(p-rr$pval[1])<1e-7
        if(!check)official_match<-FALSE
      }
      finite_all<-finite_all&&ok
      effects[[length(effects)+1L]]<-data.frame(model=model,feature_id=id,generation=if(mode=='generation')sub('case_','',term)else'Overall',
        beta_log2=b,SE=se,df=df,CI_low=lo,CI_high=hi,abundance_ratio=2^b,ratio_CI_low=2^lo,ratio_CI_high=2^hi,
        P_raw=p,analysis_valid=ok&&!bad,official_result_match=if(exists('check',inherits=FALSE))check else FALSE,
        status=if(!ok)'NOT_ESTIMABLE'else state,singular=singular,zero_replacement=min(design$relative[design$relative[,id]>0,id])/2)
    }
    if(!official_match&&finite_all)stop('Official coefficient/P mismatch: ',id)
    jp<-Fv<-dfn<-dfd<-NA_real_
    if(mode=='generation'&&!is.null(ft)&&finite_all){
      b<-lme4::fixef(ft);L<-matrix(0,2,length(b),dimnames=list(NULL,names(b)))
      L[1,terms[1:2]]<-c(1,-1);L[2,terms[2:3]]<-c(1,-1)
      jt<-tryCatch(lmerTest::contestMD(ft,L,ddf='Satterthwaite'),error=function(e){err<<-paste(err,conditionMessage(e));NULL})
      if(!is.null(jt)){jp<-jt[['Pr(>F)']][1];Fv<-jt[['F value']][1];dfn<-jt[['NumDF']][1];dfd<-jt[['DenDF']][1]}
      joint[[length(joint)+1L]]<-data.frame(model=model,feature_id=id,F_statistic=Fv,numerator_df=dfn,denominator_df=dfd,
        P_raw=jp,analysis_valid=!bad&&all(is.finite(c(jp,Fv,dfn,dfd)))&&dfn==2&&dfd>0,
        status=if(!is.finite(jp))'NOT_ESTIMABLE'else state,singular=singular)
    }else if(mode=='generation')joint[[length(joint)+1L]]<-data.frame(model=model,feature_id=id,
      F_statistic=NA_real_,numerator_df=2,denominator_df=NA_real_,P_raw=NA_real_,analysis_valid=FALSE,status=state,singular=singular)
    diagnostics[[length(diagnostics)+1L]]<-data.frame(model=model,feature_id=id,status=state,singular=singular,
      optimizer_code=paste(opt,collapse=';'),minimum_Hessian_eigenvalue=hmin,family_ICC=icc,residual_skewness=skew,
      warnings=paste(warn,collapse='; '),error=err,official_result_match=official_match)
  }
  list(effects=ov_bind(effects),heterogeneity=ov_bind(joint),diagnostics=ov_bind(diagnostics))
}

ade_bh <- function(t,family,total=nrow(t)) {
  if(!nrow(t))return(t)
  good<-tf_bool(t$analysis_valid)&is.finite(t$P_raw)&t$P_raw>=0&t$P_raw<=1
  pv<-rep(1,nrow(t));pv[good]<-t$P_raw[good]
  t$P_BH<-NA_real_;qq<-p.adjust(pv,method='BH',n=total);t$P_BH[good]<-qq[good]
  t$FDR_lt_005<-good&!is.na(t$P_BH)&t$P_BH<.05;t$FDR_family<-family;t$FDR_tests<-total;t
}

ade_fit_set <- function(x,sc,model,mode,cfg,out,cache,probe=FALSE) {
  d<-v3_design(x,sc$mask,sc$fields);pv<-d$filter$prevalence
  vari<-vapply(seq_len(ncol(d$relative)),function(j)length(unique(d$relative[,j]))>1L,logical(1))
  d$filter$retained<-d$filter$retained&vari;d$filter$nonzero_variance<-vari
  ids<-d$filter$feature_id[d$filter$retained];if(probe)ids<-head(ids,2)
  md<-d$design;md$fatty_liver<-NULL
  terms<-if(mode=='generation')paste0('case_G',0:2)else'case_overall'
  if(mode=='generation')for(g in c('G0','G1','G2'))md[[paste0('case_',g)]]<-as.integer(d$metadata$generation==g)*d$metadata$fatty_liver else
    md$case_overall<-d$metadata$fatty_liver
  mm<-model.matrix(reformulate(c(terms,d$nuisance)),md);stopifnot(qr(mm)$rank==ncol(mm))
  mapping<-data.frame(model=model,feature_id=ids,package_feature=sprintf('f%05d',seq_along(ids)))
  signature<-digest::digest(list(protocol=ade_protocol,relative=d$relative[,ids,drop=FALSE],design=md,
    mode=mode,model=model,extract=body(ade_extract),package=as.character(packageVersion('Maaslin2'))),algo='sha256')
  cp<-file.path(cache,substr(signature,1,14));dir.create(cp,recursive=TRUE,showWarnings=FALSE)
  batches<-split(seq_along(ids),ceiling(seq_along(ids)/20L));pieces<-progress<-list()
  started<-Sys.time()
  for(i in seq_along(batches)){
    rows<-batches[[i]];keys<-mapping$package_feature[rows];bi<-ids[rows]
    f<-file.path(cp,sprintf('b%03d.rds',i));reuse<-file.exists(f)
    if(reuse)piece<-readRDS(f)else{
      ade_log(model,' ',mode,': batch ',i,'/',length(batches),' (',length(rows),' features; N=',nrow(md),')')
      A<-d$relative[,bi,drop=FALSE];colnames(A)<-keys
      expected<-vapply(keys,function(k)get('LOG',asNamespace('Maaslin2'))(A[,k]),numeric(nrow(A)))
      rownames(expected)<-rownames(A);colnames(expected)<-keys
      dir<-file.path(cp,sprintf('m%03d',i));dir.create(dir,recursive=TRUE,showWarnings=FALSE)
      con<-file(file.path(dir,'console.log'),'wt');sink(con)
      ans<-tryCatch(Maaslin2::Maaslin2(input_data=as.data.frame(A),input_metadata=md,output=file.path(dir,'fit'),
        fixed_effects=c(terms,d$nuisance),random_effects='family_id',normalization='NONE',transform='LOG',standardize=FALSE,
        min_prevalence=0,min_abundance=0,min_variance=0,analysis_method='LM',cores=1,correction='BH',
        plot_heatmap=FALSE,plot_scatter=FALSE,save_models=TRUE),finally={sink();close(con)})
      if(is.null(ans$fits))stop('Official package did not return saved mixed-model fits')
      piece<-ade_extract(ans,keys,expected,d,model,mode,bi)
      saveRDS(piece,f);rm(ans,expected);gc(FALSE)
    }
    pieces[[i]]<-piece
    progress[[i]]<-data.frame(model=model,mode=mode,batch=i,total_batches=length(batches),
      completed_features=sum(lengths(batches[seq_len(i)])),eligible_features=length(ids),reused=reuse,
      elapsed_minutes=as.numeric(difftime(Sys.time(),started,units='mins')))
    fl_write_csv(ov_bind(progress),file.path(out,'audit',paste0('progress_',model,'_',mode,'.csv')))
  }
  eff<-ov_bind(lapply(pieces,`[[`,'effects'));het<-ov_bind(lapply(pieces,`[[`,'heterogeneity'))
  eff<-ade_bh(eff,paste0(model,'_',mode,'_all_feature_by_generation_effects'),nrow(eff))
  het<-ade_bh(het,paste0(model,'_all_eligible_feature_joint_generation_tests'),nrow(het))
  for(v in c('N','cases','families')){
    if(nrow(het))het[[v]]<-switch(v,N=nrow(d$metadata),cases=sum(d$metadata$fatty_liver),families=length(unique(d$metadata$family_id)))
    eff[[v]]<-vapply(eff$generation,function(g){ii<-if(g=='Overall')rep(TRUE,nrow(d$metadata))else d$metadata$generation==g
      switch(v,N=sum(ii),cases=sum(d$metadata$fatty_liver[ii]),families=length(unique(d$metadata$family_id[ii])))},numeric(1))
  }
  ade_table(d$filter,paste0('Feature_filter_',model,'_',mode,'_EN'),out)
  fl_write_csv(d$mapping,file.path(out,'audit',paste0('design_',model,'_',mode,'.csv')))
  fl_write_csv(mapping,file.path(out,'audit',paste0('feature_map_',model,'_',mode,'.csv')))
  counts<-as.data.frame(table(generation=d$metadata$generation,steatosis=d$metadata$fatty_liver));counts$model<-model
  list(effects=eff,heterogeneity=het,diagnostics=ov_bind(lapply(pieces,`[[`,'diagnostics')),
    counts=counts,filter=d$filter,signature=signature,sample_ids=as.character(d$metadata$sample_id),mode=mode,model=model)
}

ade_analyze <- function(kind,x,cfg,out,cache,probe=FALSE) {
  sc<-v3_micro_samples(x);fits<-list()
  for(nm in names(sc)){
    fits[[nm]]<-ade_fit_set(x,sc[[nm]],nm,if(kind=='A')'generation'else'overall',cfg,out,cache,probe)
    if(probe)break
  }
  if(kind=='D'&&!probe)fits[['B1_generation']]<-ade_fit_set(x,sc$B1,'B1','generation',cfg,out,cache,FALSE)
  eff<-ov_bind(lapply(fits,`[[`,'effects'));het<-ov_bind(lapply(fits,`[[`,'heterogeneity'))
  anno<-if(kind=='A')ade_labels(x$features[match(x$candidates,x$features$feature_id),,drop=FALSE])else x$features
  for(nm in c('effects','heterogeneity')){
    z<-if(nm=='effects')eff else het
    if(nrow(z)){at<-match(z$feature_id,anno$feature_id);for(v in setdiff(names(anno),'feature_id'))z[[v]]<-anno[[v]][at]}
    if(nm=='effects')eff<-z else het<-z
  }
  r<-list(protocol=ade_protocol,kind=kind,fits=fits,effects=eff,heterogeneity=het,annotation=anno,
    diagnostics=ov_bind(lapply(fits,`[[`,'diagnostics')),counts=ov_bind(lapply(fits,`[[`,'counts')))
  ade_table(eff,'All_feature_effects_EN',out);ade_table(het,'All_joint_generation_tests_EN',out)
  ade_table(r$diagnostics,'All_model_diagnostics_EN',out);ade_table(anno,'Feature_annotations_EN',out)
  ade_table(r$counts,'Analysis_sample_counts_EN',out)
  summary<-ov_bind(lapply(fits,function(f){h<-f$heterogeneity;e<-f$effects
    data.frame(model=f$model,mode=f$mode,participants=length(f$sample_ids),eligible_features=sum(f$filter$retained),
      valid_effects=sum(tf_bool(e$analysis_valid)),effect_tests=nrow(e),significant_effects=sum(e$FDR_lt_005),
      valid_heterogeneity_tests=if(nrow(h))sum(tf_bool(h$analysis_valid))else NA_integer_,
      significant_heterogeneity=if(nrow(h))sum(h$FDR_lt_005)else NA_integer_)}))
  r$summary<-summary;ade_table(summary,'Analysis_model_summary_EN',out)
  if(all(c('Adult_B1_same_sample','B2')%in%names(fits)))stopifnot(identical(fits$Adult_B1_same_sample$sample_ids,fits$B2$sample_ids))
  r
}

ade_main <- function(kind,args=character(),root=ade_root()) {
  stopifnot(kind%in%c('A','D','E'))
  if(getRversion()<'4.5'){
    entry<-c(A='41_A_generation_heterogeneity.R',D='42_D_functional_pathways.R',E='43_E_strain_sharing.R')[[kind]]
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts',entry)),args))
    if(status!=0)stop('Project-R ADE entry failed');return(invisible(NULL))
  }
  cfg<-ade_setup(root);key<-c(A='generation_heterogeneity',D='functional_pathways',E='strain_sharing')[[kind]]
  cfg$output_module<-unname(cfg$output_modules[key]);cfg$run_id<-unname(cfg$run_ids[key]);cfg$run_id_explicit<-FALSE
  if(any(!args%in%c('--run','--probe','--figures-only','--check')))stop('Supported: --run, --probe, --figures-only, --check')
  src<-tf_source(root);x<-tf_load(src)
  if('--check'%in%args){cat('PASS: source scope, original metadata, required packages and Matrix ABI\n');return(invisible(TRUE))}
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(p in c('private','audit','tables','figures','reports'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
  writeLines(c('RUNNING',ade_protocol,kind),file.path(out,'run_status.txt'))
  old_log<-getOption('fatty_liver.ade_log');options(fatty_liver.ade_log=file.path(out,'audit/execution.log'))
  on.exit(options(fatty_liver.ade_log=old_log),add=TRUE)
  tryCatch({
    ade_log('START ',kind,' ',cfg$run_id)
    protected<-ade_fingerprint(src$files);fl_write_csv(protected,file.path(out,'audit/parent_fingerprints.csv'))
    if(kind=='E'){
      sys.source(file.path(root,'R/strain_sharing_extension.R'),envir=.GlobalEnv)
      result<-ade_strain_main(x,cfg,out)
    }else if('--figures-only'%in%args){
      module<-file.path(root,'outputs',cfg$output_module)
      rds<-list.files(file.path(module,'current'),pattern='^analysis_result.rds$',recursive=TRUE,full.names=TRUE)
      if(!length(rds))rds<-list.files(file.path(module,'history'),pattern='^analysis_result.rds$',recursive=TRUE,full.names=TRUE)
      if(!length(rds))stop('No saved analysis result for figure recovery')
      file<-rds[which.max(file.info(rds)$mtime)];result<-readRDS(file)
      protected<-result$protected
      previous<-dirname(dirname(file))
      for(folder in c('tables','reports','audit')){
        fs<-list.files(file.path(previous,folder),full.names=TRUE,all.files=TRUE,no..=TRUE)
        if(length(fs)&&!all(file.copy(fs,file.path(out,folder),recursive=TRUE,overwrite=TRUE,copy.date=TRUE)))
          stop('Cannot preserve saved numerical tables/reports during figure recovery')
      }
      for(name in c('functional_input.rds','humann_row_indices.json','humann_community.json')){
        p<-file.path(previous,'private',name)
        if(file.exists(p)&&!file.copy(p,file.path(out,'private',name),overwrite=TRUE))stop('Cannot preserve functional import audit')
      }
    }else{
      if(kind=='D'){imp<-ade_import_paths(cfg,x,out);x<-imp$x;protected<-ov_bind(list(protected,imp$protected))}
      cache<-file.path(root,'outputs',cfg$output_module,'_checkpoints/private');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
      result<-ade_analyze(kind,x,cfg,out,cache,'--probe'%in%args);result$protected<-protected
    }
    saveRDS(result,file.path(out,'private/analysis_result.rds'));saveRDS(cfg,file.path(out,'private/run_config.rds'))
    if('--probe'%in%args){
      ade_log('PROBE COMPLETE; not adopted or published')
      writeLines(c('PROBE_COMPLETE',ade_protocol),file.path(out,'run_status.txt'));return(invisible(result))
    }
    if(kind!='E'){
      sys.source(file.path(root,'R/generation_ADE_figures.R'),envir=.GlobalEnv)
      ade_reports(result,cfg,out);ade_figures(result,cfg,out)
    }
    ade_protect(protected);capture.output(sessionInfo(),file=file.path(out,'audit/sessionInfo.txt'))
    writeLines(c('COMPLETED',ade_protocol,if(kind=='E'&&identical(result$status,'INPUTS_REQUIRED'))'Data availability audit complete; strain-sharing inference NOT RUN: missing sequence-based inputs'else'Analysis and rendered QA completed'),file.path(out,'run_status.txt'))
    # Serialize only publication/gallery refresh when independent A and D runs
    # finish close together. No model or source-data lock is needed.
    lock<-file.path(root,'inventory/ADE_20261101/publication.lock')
    dir.create(dirname(lock),recursive=TRUE,showWarnings=FALSE);start_wait<-Sys.time()
    while(!dir.create(lock,showWarnings=FALSE)){
      if(as.numeric(difftime(Sys.time(),start_wait,units='secs'))>240)stop('Another ADE run is still publishing; rerun with --figures-only')
      Sys.sleep(.5)
    }
    final<-tryCatch({
      target<-fl_output_complete_run(ctx)
      writeLines(cfg$run_id,file.path(root,'outputs',cfg$output_module,'LATEST_RUN.txt'));target
    },finally={
      resolved<-normalizePath(lock,winslash='/',mustWork=TRUE)
      base<-normalizePath(file.path(root,'inventory/ADE_20261101'),winslash='/',mustWork=TRUE)
      stopifnot(identical(tolower(dirname(resolved)),tolower(base)),basename(resolved)=='publication.lock')
      unlink(resolved,recursive=TRUE)
    })
    # Publication moves _running to current; append to the adopted log path.
    options(fatty_liver.ade_log=file.path(final,'audit/execution.log'))
    ade_log('COMPLETED: ',final);invisible(final)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
