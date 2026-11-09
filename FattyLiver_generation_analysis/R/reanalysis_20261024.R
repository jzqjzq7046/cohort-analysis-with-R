# Active protocol: user-approved C1/C2/B1/B2/P1/P2, 2026-10-24.
# Ultrasound steatosis only. Primary labels are not verified pedigrees/households.
v3_protocol <- 'two_generation_C1_C2_B1_B2_P1_P2_20261024'
v3_pairs <- list(c('G0','G1'),c('G1','G2'),c('G0','G2'))
v3_B1 <- c('generation','age','sex','prep_batch','seq_pool','season_sin','season_cos','abx_recent')
v3_B2 <- c(v3_B1,'smoking','alcohol','diet_score','physical_activity','energy_log','BMI18')
v3_P2 <- setdiff(v3_B1,'generation')
v3_log <- function(...) {cat(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),paste0(...,collapse=''),'\n');flush.console()}
v3_cc <- function(m,fields) {j<-complete.cases(m[,fields,drop=FALSE]);for(v in fields)if(is.numeric(m[[v]]))j<-j&is.finite(m[[v]]);j}
v3_scope <- function(m,g=NULL) {
 m$family_id<-fl_primary_family_key(m$family_id)
 valid<-!is.na(m$family_id)&m$generation%in%c('G0','G1','G2')&!is.na(m$fatty_liver)
 z<-m[valid,,drop=FALSE];tab<-table(z$family_id,factor(z$generation,levels=c('G0','G1','G2')))
 families<-if(is.null(g))rownames(tab)[rowSums(tab>0)>=2]else rownames(tab)[rowSums(tab[,g,drop=FALSE]>0)==2]
 keep<-valid&m$family_id%in%families
 if(!is.null(g))keep<-keep&m$generation%in%g
 droplevels(m[keep,,drop=FALSE])
}
v3_table <- function(z,name,out) {
 p<-file.path(out,'tables',paste0(name,'.csv'));dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);fl_write_csv(z,p)
 if(nrow(z)<=150L) {
  esc<-function(x){x<-gsub('&','&amp;',as.character(x),fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);x[is.na(x)]<-'';x}
  lines<-c('<!doctype html><meta charset="utf-8"><style>body{font:13px Arial;margin:25px}table{border-collapse:collapse}th,td{padding:6px 9px;border-bottom:1px solid #bbb;text-align:left}th{border-top:2px solid black;border-bottom:1px solid black}tr:last-child td{border-bottom:2px solid black}</style>',paste0('<h3>',name,'</h3><table><tr>',paste0('<th>',esc(names(z)),'</th>',collapse=''),'</tr>'))
  for(i in seq_len(nrow(z)))lines<-c(lines,paste0('<tr>',paste0('<td>',esc(unlist(z[i,,drop=FALSE])),'</td>',collapse=''),'</tr>'))
  writeLines(c(lines,'</table>'),sub('csv$','html',p),useBytes=TRUE)
 }
 invisible(z)
}
v3_counts <- function(m,label) do.call(rbind,lapply(c('Overall','G0','G1','G2','Minor','Adult'),function(g){
 j<-if(g=='Overall')rep(TRUE,nrow(m))else if(g=='Minor')!is.na(m$age)&m$age<18 else if(g=='Adult')!is.na(m$age)&m$age>=18 else m$generation==g
 z<-m[j,,drop=FALSE];data.frame(cohort=label,group=g,N=nrow(z),families=length(unique(z$family_id)),cases=sum(z$fatty_liver),prevalence_percent=if(nrow(z))100*mean(z$fatty_liver)else NA_real_)
}))
v3_prepare <- function(cfg,out) {
 v3_log('Preparing expanded clinical and microbiome cohorts')
 fields<-unique(c(setdiff(fl_requested_fields(),c('HOMA18','INSU18')),'alue18','kaumaa18','PAI18'))
 csv<-file.path(out,'private','clinical_selected.csv')
 if(!file.exists(csv)){
  status<-system2(cfg$python,c(shQuote(file.path(cfg$project_root,'scripts/export_sas_for_r.py')),'--source',shQuote(cfg$clinical_file),'--out',shQuote(csv),'--fields',shQuote(paste(fields,collapse=','))))
  if(status!=0)stop('Clinical SAS export failed')
 }
 raw<-read.csv(csv,colClasses='character',fileEncoding='UTF-8-BOM',check.names=FALSE)
 cl<-data.frame(participant_id=fl_text(fl_column(raw,'tutknom18',TRUE)),family_id=fl_primary_family_key(fl_column(raw,'tutkno80perhe',TRUE)),generation=toupper(fl_text(fl_column(raw,'sukupolvi',TRUE))))
 stopifnot(!anyDuplicated(cl$participant_id),!anyNA(cl$participant_id))
 cl$sample_id<-paste0('clinical_',cl$participant_id)
 cl$fat_code<-fl_number(fl_column(raw,'RASVAMAKSA18',TRUE));cl$fatty_liver<-ifelse(cl$fat_code%in%1:3,as.integer(cl$fat_code%in%c(2,3)),NA_integer_)
 cl$age<-fl_age(fl_date(fl_column(raw,'syntymapvm')),fl_date(fl_column(raw,'MAKSAPVM18')))
 for(v in c('ika_status18','TARKKAIKA18')){a<-fl_number(fl_column(raw,v));j<-is.na(cl$age)&is.finite(a)&a>=0&a<=120;cl$age[j]<-a[j]}
 sp<-fl_number(fl_column(raw,'sp',TRUE));stopifnot(all(is.na(sp)|sp%in%c(1,2)))
 cl$sex<-ifelse(is.na(sp),NA_character_,ifelse(sp==1,'Female','Male'))
 for(v in c('BMI18','alkodpd18','dietscore_g18','PITUUS18','PAINO18','WAIST18','syst18','DKVAI18','GLUK18','HbA1c18','TOTKOL18','TRIGLY18','HDLKOL18','LDLKOL18','CRP18','ALAT18','ASAT18','S_GT18'))cl[[v]]<-fl_number(fl_column(raw,v))
 cl$age_group<-ifelse(is.na(cl$age),'Age unavailable',ifelse(cl$age<18,'Minor','Adult'));cl$alcohol<-cl$alkodpd18
 all_cl<-cl[cl$generation%in%c('G0','G1','G2')&!is.na(cl$fatty_liver),,drop=FALSE]
 cl<-v3_scope(all_cl)
 # Unrestricted historical processed package contains all 4,623 ultrasound-eligible microbiome samples.
 folder<-file.path(cfg$project_root,'data/processed/processed_microbiome_v1')
 path<-file.path(folder,'microbiome_processed.rds');mf<-read.csv(file.path(folder,'dataset_manifest.csv'),fileEncoding='UTF-8-BOM')
 stopifnot(unname(tools::md5sum(path))==mf$md5[mf$file=='microbiome_processed.rds'])
 d<-readRDS(path);m<-d$metadata
 stopifnot(identical(colnames(d$full$counts),m$sample_id),!anyDuplicated(m$sample_id),all(d$full$counts>=0),all(d$full$counts==floor(d$full$counts)))
 k<-match(as.character(m$participant_id),all_cl$participant_id)
 stopifnot(!anyNA(k),all(m$fat_code==all_cl$fat_code[k]),all(as.character(m$generation)==all_cl$generation[k]),all(as.character(m$sex)==all_cl$sex[k]),max(abs(m$age-all_cl$age[k]),na.rm=TRUE)<1e-6)
 source_cd<-fl_cv_coldata(cfg$microbiome_file);j<-match(m$sample_id,rownames(source_cd));stopifnot(!anyNA(j),!anyDuplicated(rownames(source_cd)))
 cd<-source_cd[j,,drop=FALSE]
 m$family_id<-fl_primary_family_key(m$family_id)
 stopifnot(all(is.na(m$family_id)&is.na(all_cl$family_id[k])|m$family_id==all_cl$family_id[k],na.rm=TRUE))
 m$abx_recent<-as.numeric(cd$ANTBV18);m$smoking<-as.numeric(cd$smoke18)
 stopifnot(all(is.na(m$abx_recent)|m$abx_recent%in%0:1),all(is.na(m$smoking)|m$smoking%in%0:1))
 m$alcohol<-m$alkodpd18;m$diet_score<-m$dietscore_g18;m$physical_activity<-as.numeric(cd$PAI18);m$physical_activity[is.na(m$age)|m$age<18]<-NA_real_
 month<-as.integer(format(cd$ULOSTE_PVM18,'%m'));m$season_sin<-sin(2*pi*(month-1)/12);m$season_cos<-cos(2*pi*(month-1)/12)
 ns<-fl_cv_coldata('E:/YFS/workspace/microbiome/tse_insulin.rds');stopifnot(!anyDuplicated(ns$SampleID))
 nj<-match(m$sample_id,as.character(ns$SampleID));kcal<-as.numeric(ns$ENERCK18[nj]);stopifnot(all(is.na(kcal)|kcal>0));m$energy_log<-log10(kcal)
 fl_write_csv(data.frame(source='ENERCK18',microbiome_samples=nrow(m),unmatched_source_IDs=sum(is.na(nj)),action='Retain as missing; only adult B2 complete-case eligibility is affected'),file.path(out,'audit/nutrient_linkage.csv'))
 for(v in c('WAIST18','syst18','DKVAI18','GLUK18','HbA1c18','TRIGLY18','HDLKOL18','ALAT18','ASAT18'))m[[v]]<-all_cl[[v]][k]
 m<-v3_scope(m);ix<-match(m$sample_id,colnames(d$full$counts));counts<-d$full$counts[,ix,drop=FALSE]
 counts<-counts[rowSums(counts)>0,,drop=FALSE];relative<-sweep(counts,2,colSums(counts),'/')
 stopifnot(max(abs(colSums(relative)-1))<1e-10,max(abs(colSums(counts)-m$matrix_total))==0)
 m$shannon<--colSums(ifelse(relative>0,relative*log(pmax(relative,.Machine$double.xmin)),0));m$observed_sgbs<-colSums(counts>0)
 for(v in c('generation','sex','prep_batch','seq_pool'))m[[v]]<-factor(m[[v]])
 m$generation<-factor(m$generation,levels=c('G0','G1','G2'));m$family_id<-factor(m$family_id)
 ft<-d$full$features[match(rownames(counts),d$full$features$feature_id),,drop=FALSE];stopifnot(!anyNA(ft$feature_id))
 # Candidate union is rebuilt on the expanded baseline, never inherited from the old 319-family analysis.
 pv<-sapply(c('G0','G1','G2'),function(g)rowMeans(relative[,m$generation==g,drop=FALSE]>0))
 candidate<-apply(pv>=.10,1,any)
 filter<-data.frame(feature_id=rownames(counts),prevalence_overall=rowMeans(relative>0),pv,candidate=candidate)
 v3_table(filter,'00_Baseline_feature_filter_EN',out)
 paths<-c(cfg$clinical_file,path,file.path(folder,'dataset_manifest.csv'),cfg$microbiome_file,'E:/YFS/workspace/microbiome/tse_insulin.rds')
 fingerprints<-data.frame(path=paths,md5=unname(tools::md5sum(paths)));fl_write_csv(fingerprints,file.path(out,'audit/input_fingerprints.csv'))
 keep_B1<-v3_cc(m,c('family_id','fatty_liver',v3_B1));keep_B2<-!is.na(m$age)&m$age>=18&v3_cc(m,c('family_id','fatty_liver',v3_B2))
 v3_table(rbind(v3_counts(all_cl,'All valid clinical ultrasounds before family restriction'),v3_counts(cl,'Clinical: >=2 observed generations'),v3_counts(m,'Microbiome: >=2 observed generations'),v3_counts(m[keep_B1,],'B1 complete cases'),v3_counts(m[keep_B2,],'B2 adult common complete cases')),'00_Cohort_counts_EN',out)
 definitions<-data.frame(set=c('C1','C2_alcohol','C2_BMI','B1','B2','P1','P2'),variables=c('younger age + sex + observed older count when variable','C1 + younger alcohol','C1 + younger alcohol + continuous adult BMI',paste(v3_B1,collapse=' + '),paste(v3_B2,collapse=' + '),'No regression covariates; generation-role node permutations',paste(v3_P2,collapse=' + ')))
 v3_table(definitions,'00_Model_definitions_EN',out)
 miss<-do.call(rbind,lapply(unique(c(v3_B2,'family_id')),function(v)data.frame(variable=v,N=nrow(m),missing=sum(is.na(m[[v]])),adult_missing=sum(is.na(m[[v]])&m$age>=18,na.rm=TRUE))))
 v3_table(miss,'00_Microbiome_covariate_missingness_EN',out)
 flow<-data.frame(stage=c('Clinical source records','Clinical valid ultrasound/generation','Clinical >=2 generation family cohort','Microbiome valid ultrasound/QC','Microbiome >=2 generation family cohort','B1 complete cases','B2 adult common complete cases'),N=c(nrow(raw),nrow(all_cl),nrow(cl),nrow(d$metadata),nrow(m),sum(keep_B1),sum(keep_B2)))
 v3_table(flow,'00_Sample_flow_EN',out)
 result<-list(clinical=cl,metadata=m,counts=counts,relative=relative,features=ft,candidates=rownames(counts)[candidate],flow=flow,input_fingerprints=fingerprints,protocol=v3_protocol)
 saveRDS(result,file.path(out,'private','expanded_analysis_data.rds'),compress='gzip')
 v3_log('Cohorts ready: clinical ',nrow(cl),'; microbiome ',nrow(m),'; B1 ',sum(keep_B1),'; B2 ',sum(keep_B2),'; candidate SGBs ',sum(candidate))
 result
}
v3_descriptive <- function(x,cfg,out) {
 rows<-list();fields<-c(age='Age, years',BMI18='BMI, kg/m2',WAIST18='Waist, cm',syst18='Systolic BP, mmHg',DKVAI18='Diastolic BP, mmHg',GLUK18='Glucose, mmol/L',HbA1c18='HbA1c, mmol/mol',TRIGLY18='Triglycerides, mmol/L',HDLKOL18='HDL cholesterol, mmol/L',ALAT18='ALT, U/L',ASAT18='AST, U/L',alkodpd18='Alcohol, servings/day',dietscore_g18='Diet score, 0-27')
 for(type in c('Clinical','Microbiome')){
  m<-if(type=='Clinical')x$clinical else x$metadata;wide<-data.frame(Characteristic=c(unname(fields),'Female, n/N (%)','Hepatic steatosis, n/N (%)','Age <18, n/N (%)'),check.names=FALSE)
  for(g in c('Overall','G0','G1','G2')){
   z<-if(g=='Overall')m else m[m$generation==g,,drop=FALSE];values<-vapply(names(fields),function(v){a<-z[[v]];a<-a[is.finite(a)];if(!length(a))return('Not available');q<-quantile(a,c(.25,.5,.75));sprintf('%.1f [%.1f, %.1f]; n=%d',q[2],q[1],q[3],length(a))},character(1))
   binary<-function(a){n<-sum(!is.na(a));sprintf('%d/%d (%.1f%%)',sum(a,na.rm=TRUE),n,100*mean(a,na.rm=TRUE))}
   values<-c(values,binary(z$sex=='Female'),binary(z$fatty_liver==1),binary(z$age<18));wide[[paste0(g,' (N=',nrow(z),')')]]<-values
  };v3_table(wide,paste0('Table1_',type,'_by_generation_EN'),out);rows[[type]]<-wide
 }
 z<-x$metadata;qc<-data.frame(sample_count=nrow(z),full_nonzero_SGBs=nrow(x$relative),candidate_union_SGBs=length(x$candidates),relative_sum_max_error=max(abs(colSums(x$relative)-1)),input_integer_counts=all(x$counts==floor(x$counts)),low_total_flags=sum(z$flag_low_total),low_feature_flags=sum(z$flag_low_features))
 v3_table(qc,'01_Microbiome_QC_EN',out);rows
}
v3_clinical <- function(x,cfg,out) {
 cfg$cg_alpha<-.05;cfg$cg_max_iterations<-100L;cfg$cg_working_correlation<-'independence'
 models<-desc<-flows<-list()
 for(g in v3_pairs){key<-paste(g,collapse='-');b<-fl_cg_make_comparison(v3_scope(x$clinical,g),g);d<-b$data
  core<-v3_cc(d,c('age','sex'));adult<-core&!is.na(d$age)&d$age>=18&v3_cc(d,c('alkodpd18','BMI18'))
  specs<-list(Unadjusted=list(core,character()),C1=list(core,c('age','sex','n_older_observed')),Adult_C1=list(adult,c('age','sex','n_older_observed')),C2_alcohol=list(adult,c('age','sex','n_older_observed','alkodpd18')),C2_BMI=list(adult,c('age','sex','n_older_observed','alkodpd18','BMI18')))
  for(nm in names(specs)){
   z<-droplevels(d[specs[[nm]][[1]],,drop=FALSE]);spec<-data.frame(model=nm,sample_set=if(nm%in%c('Unadjusted','C1'))'all_age_same_sample'else'adult_same_sample',covariates=paste(specs[[nm]][[2]],collapse='+'),primary=nm=='C1',references='Long2019;LiangZeger1986',study_adaptation='Primary label; not verified parent-child. Pair-specific observed generations; no third-generation requirement.')
   r<-fl_cg_fit(z,spec,key,cfg);models[[paste(key,nm)]]<-r$result
  }
  desc[[key]]<-fl_cg_describe(d,key);flows[[key]]<-b$flow
  saveRDS(d,file.path(out,'private',paste0('clinical_pairs_',key,'.rds')))
 }
 result<-list(models=ov_bind(models),descriptive=ov_bind(desc),flow=ov_bind(flows))
 v3_table(result$models,'Table2_Clinical_cross_generation_all_models_EN',out);v3_table(result$descriptive,'Clinical_observed_prevalence_EN',out);v3_table(result$flow,'Clinical_pair_sample_flow_EN',out)
 v3_log('Clinical models complete');result
}
v3_nuisance <- function(m,fields,age_spline=TRUE) {
 m<-droplevels(m);terms<-setdiff(fields,'age');knots<-list()
 if('age'%in%fields){
  if(age_spline){groups<-if('generation'%in%fields)levels(droplevels(factor(m$generation)))else 'Overall'
   for(g in groups){j<-if(g=='Overall')seq_len(nrow(m))else which(m$generation==g);df<-min(3L,length(unique(m$age[j]))-1L);if(df<1)next
    b<-if(df==1)matrix(m$age[j],ncol=1)else splines::ns(m$age[j],df=df)
    for(k in seq_len(ncol(b))){v<-paste0('age_',g,'_',k);m[[v]]<-0;m[j,v]<-b[,k];terms<-c(terms,v)}
   }
  }else terms<-c(terms,'age')
 }
 constant<-terms[vapply(m[,terms,drop=FALSE],function(v)length(unique(v))<2,logical(1))]
 Z<-model.matrix(reformulate(setdiff(terms,constant)),m);q<-qr(Z,tol=1e-9);keep<-sort(q$pivot[seq_len(q$rank)]);aliases<-setdiff(colnames(Z),colnames(Z)[keep]);Z<-Z[,keep,drop=FALSE]
 list(metadata=m,Z=Z,nuisance=setdiff(terms,constant),dropped=c(constant,aliases))
}
v3_micro_samples <- function(x) {
 m<-x$metadata;B1<-v3_cc(m,c('fatty_liver','family_id',v3_B1));B2<-!is.na(m$age)&m$age>=18&v3_cc(m,c('fatty_liver','family_id',v3_B2))
 list(B1=list(mask=B1,fields=v3_B1),Adult_B1_same_sample=list(mask=B2,fields=v3_B1),B2=list(mask=B2,fields=v3_B2))
}
v3_design <- function(x,mask,fields) {
 m<-droplevels(x$metadata[mask,,drop=FALSE]);m<-m[order(m$family_id,m$sample_id),,drop=FALSE]
 nz<-v3_nuisance(m,fields);m<-nz$metadata;Z<-nz$Z[,colnames(nz$Z)!='(Intercept)',drop=FALSE]
 centers<-colMeans(Z);scales<-apply(Z,2,sd);Z<-sweep(sweep(Z,2,centers,'-'),2,scales,'/');map<-data.frame(column=sprintf('nuis_%03d',seq_len(ncol(Z))),original=colnames(Z),center=centers,scale=scales);colnames(Z)<-map$column
 md<-data.frame(family_id=factor(m$family_id),fatty_liver=factor(m$fatty_liver,levels=0:1,labels=c('No','Yes')),Z,row.names=m$sample_id,check.names=FALSE)
 mm<-model.matrix(reformulate(c('fatty_liver',map$column)),md);stopifnot(qr(mm)$rank==ncol(mm))
 ids<-x$candidates;A<-t(x$relative[ids,match(m$sample_id,colnames(x$relative)),drop=FALSE]);pv<-colMeans(A>0)
 list(metadata=m,design=md,relative=A,mapping=map,bases='fatty_liver',terms='fatty_liverYes',nuisance=map$column,formula=paste('~ fatty_liver +',paste(map$column,collapse=' + '),'+ (1 | family_id)'),filter=data.frame(feature_id=ids,prevalence=pv,retained=pv>=.10),dropped=nz$dropped)
}
v3_alpha <- function(x,cfg,out) {
 result<-list();interactions<-list()
 for(nm in names(v3_micro_samples(x))){s<-v3_micro_samples(x)[[nm]];d<-v3_design(x,s$mask,s$fields);m<-d$metadata;md<-d$design;md$shannon<-m$shannon
  for(g in c('G0','G1','G2'))md[[paste0('case_',g)]]<-as.integer(m$generation==g)*m$fatty_liver
  case<-paste0('case_',c('G0','G1','G2'));forms<-list(Overall=reformulate(c('fatty_liver',d$nuisance),response='shannon'),Generation=reformulate(c(case,d$nuisance),response='shannon'))
  for(type in names(forms)){
   f<-forms[[type]];mm<-model.matrix(f,md);stopifnot(qr(mm)$rank==ncol(mm));fit<-geepack::geeglm(f,data=md,id=as.integer(factor(md$family_id)),family=gaussian,corstr='independence',std.err='san.se')
   stopifnot(fit$geese$error==0);b<-coef(fit);V<-vcov(fit);terms<-if(type=='Overall')'fatty_liverYes'else case
   for(term in terms){i<-match(term,names(b));se<-sqrt(V[i,i]);g<-if(type=='Overall')'Overall'else sub('case_','',term);j<-if(g=='Overall')rep(TRUE,nrow(m))else m$generation==g
    result[[paste(nm,type,g)]]<-data.frame(model=nm,generation=g,N=sum(j),cases=sum(m$fatty_liver[j]),families=length(unique(m$family_id[j])),beta=b[i],SE=se,CI_low=b[i]-qnorm(.975)*se,CI_high=b[i]+qnorm(.975)*se,P_raw=2*pnorm(-abs(b[i]/se)),status='ESTIMATED',effect='Shannon difference: steatosis minus no steatosis')
   }
   if(type=='Generation'){ix<-match(case,names(b));H<-rbind(c(1,-1,0),c(0,1,-1));delta<-H%*%b[ix];W<-H%*%V[ix,ix]%*%t(H);chi<-as.numeric(crossprod(delta,solve(W,delta)));interactions[[nm]]<-data.frame(model=nm,chi_square=chi,df=2,P_raw=pchisq(chi,2,lower.tail=FALSE))}
  };fl_write_csv(d$mapping,file.path(out,'audit',paste0('alpha_design_',nm,'.csv')))
 }
 result<-list(effects=ov_bind(result),interaction=ov_bind(interactions));v3_table(result$effects,'Table4_Shannon_associations_EN',out);v3_table(result$interaction,'Shannon_generation_heterogeneity_EN',out);v3_log('Shannon models complete');result
}
v3_endpoint <- function(m,fields) {
 if(!length(fields))return(matrix(1/sqrt(nrow(m)),nrow(m),1))
 z<-v3_nuisance(m,fields);qr.Q(qr(z$Z))
}
v3_pair_projection <- function(Y,Qa,Qb) {
 rowfit<-drop(Qa%*%crossprod(Qa,rowMeans(Y)));colfit<-drop(Qb%*%crossprod(Qb,colMeans(Y)))
 Y-outer(rowfit,rep(1,ncol(Y)))-outer(rep(1,nrow(Y)),colfit)+mean(Y)
}
v3_pair_test <- function(m,Y,same,fields,cfg,out,key) {
 a<-which(m$generation==strsplit(key,'-',fixed=TRUE)[[1]][1]);b<-which(m$generation==strsplit(key,'-',fixed=TRUE)[[1]][2])
 Qa<-v3_endpoint(m[a,,drop=FALSE],fields);Qb<-v3_endpoint(m[b,,drop=FALSE],fields)
 E<-v3_pair_projection(Y,Qa,Qb);X<-v3_pair_projection(1*same,Qa,Qb);ss<-sum(X^2);stopifnot(ss>1e-8)
 X<-X/sqrt(ss);num<-sum(E*X);df<-length(Y)-ncol(Qa)-ncol(Qb);rss<-sum(E^2)-num^2;observed<-num/sqrt(rss/df)
 stopifnot(max(abs(crossprod(Qa,rowMeans(E))))<1e-7,max(abs(crossprod(Qb,colMeans(E))))<1e-7)
 B<-cfg$v3_permutations;set.seed(cfg$v3_seed+match(key,c('G0-G1','G1-G2','G0-G2'))+if(length(fields))100L else 0L)
 blocks<-function(z){key<-if(length(fields))interaction(z$sex,floor(z$age/5),drop=TRUE)else factor(rep(1,nrow(z)));list(groups=unname(split(seq_len(nrow(z)),key)))}
 rows<-fl_msa_node_permutations(blocks(m[a,,drop=FALSE]),length(a),B);cols<-fl_msa_node_permutations(blocks(m[b,,drop=FALSE]),length(b),B)
 prefix<-file.path(out,'private',paste0('pair_',if(length(fields))'P2_'else'P1_',key))
 values<-list(E,X,Qa,Qb,rows,cols);suffix<-c('_E.bin','_X.bin','_QA.bin','_QB.bin','_ROW.bin','_COL.bin')
 for(i in seq_along(values)){con<-file(paste0(prefix,suffix[i]),'wb');tryCatch(writeBin(if(i<5)as.double(values[[i]])else as.integer(values[[i]]),con,endian='little'),finally=close(con))}
 status<-system2(cfg$python,c(shQuote(file.path(cfg$project_root,'scripts/two_generation_pair_kernel.py')),'--prefix',shQuote(prefix),'--na',length(a),'--nb',length(b),'--ka',ncol(Qa),'--kb',ncol(Qb),'--permutations',B,'--df',df),stdout=paste0(prefix,'.log'),stderr=paste0(prefix,'.log'))
 if(status!=0)stop('Pair kernel failed: ',prefix)
 null<-read.csv(paste0(prefix,'_t.csv'))$t_null;stopifnot(length(null)==B,all(is.finite(null)))
 # Independently recompute first five permuted projections in R.
 check<-vapply(1:5,function(i){e<-E[rows[i,],cols[i,],drop=FALSE];r<-v3_pair_projection(e,Qa,Qb);n<-sum(X*r);n/sqrt((sum(r^2)-n^2)/df)},numeric(1))
 stopifnot(max(abs(check-head(null,5)))<1e-7)
 p<-(1+sum(abs(null)>=abs(observed)-1e-10))/(B+1)
 fl_write_csv(data.frame(iteration=seq_len(B),t_null=null),paste0(prefix,'_null.csv'))
 unlink(paste0(prefix,suffix),recursive=FALSE)
 data.frame(comparison=key,model=if(length(fields))'P2'else'P1',N=nrow(m),families=length(unique(m$family_id)),older_N=length(a),younger_N=length(b),IF_pairs=sum(same),BF_pairs=sum(!same),IF_mean=mean(Y[same]),BF_mean=mean(Y[!same]),IF_minus_BF=num/sqrt(ss),t_statistic=observed,P_raw=p,permutations=B,p_resolution=1/(B+1),rank_older=ncol(Qa),rank_younger=ncol(Qb),R_python_first5_error=max(abs(check-head(null,5))),status='ESTIMATED')
}
v3_similarity <- function(x,cfg,out) {
 results<-distribution<-flow<-list()
 for(g in v3_pairs){key<-paste(g,collapse='-');baseline<-v3_scope(x$metadata,g)
  for(model in c('P1','P2')){
   fields<-if(model=='P1')character()else v3_P2
   m<-if(length(fields))v3_scope(baseline[v3_cc(baseline,fields),,drop=FALSE],g)else baseline
   A<-x$relative[,match(m$sample_id,colnames(x$relative)),drop=FALSE]
   D<-as.matrix(vegan::vegdist(t(A),method='bray'));a<-which(m$generation==g[1]);b<-which(m$generation==g[2]);Y<-D[a,b,drop=FALSE];same<-outer(as.character(m$family_id[a]),as.character(m$family_id[b]),'==')
   v3_log('Similarity ',key,' ',model,': ',length(a),' x ',length(b),' nodes')
   r<-v3_pair_test(m,Y,same,fields,cfg,out,key);results[[paste(key,model)]]<-r
   for(type in c('BF','IF')){v<-if(type=='IF')Y[same]else Y[!same];q<-quantile(v,c(0,.25,.5,.75,1));distribution[[paste(key,model,type)]]<-data.frame(comparison=key,model=model,type=type,pairs=length(v),minimum=q[1],Q1=q[2],median=q[3],Q3=q[4],maximum=q[5],mean=mean(v))}
   flow[[paste(key,model)]]<-data.frame(comparison=key,model=model,baseline_N=nrow(baseline),complete_pair_eligible_N=nrow(m),cases=sum(m$fatty_liver))
   saveRDS(list(metadata=m,same=same,distance=Y),file.path(out,'private',paste0('pair_source_',key,'_',model,'.rds')))
   rm(D,Y,A);gc(FALSE)
  }
 }
 result<-list(effects=ov_bind(results),distribution=ov_bind(distribution),flow=ov_bind(flow))
 v3_table(result$effects,'Table3_Family_microbiome_similarity_EN',out);v3_table(result$distribution,'Pair_distance_distributions_EN',out);v3_table(result$flow,'Pair_sample_flow_EN',out);result
}
v3_composition_legacy_check <- function(x,cfg,out) {
 results<-audit<-list()
 for(nm in names(v3_micro_samples(x))){
  sc<-v3_micro_samples(x)[[nm]];d<-v3_design(x,sc$mask,sc$fields);m<-d$metadata;md<-d$design;md$fatty_liver<-m$fatty_liver
  A<-x$relative[,match(m$sample_id,colnames(x$relative)),drop=FALSE];D<-vegan::vegdist(t(A),method='bray')
  # Remove only exactly redundant condition columns after adding family fixed effects.
  Z<-model.matrix(reformulate(c('family_id',d$nuisance)),md);q<-qr(Z,tol=1e-9);keep<-q$pivot[seq_len(q$rank)];Z<-Z[,sort(keep),drop=FALSE]
  kept_nuis<-d$nuisance[d$nuisance%in%colnames(Z)];condition<-c('family_id',kept_nuis)
  Z<-model.matrix(reformulate(condition),md);stopifnot(qr(Z)$rank==ncol(Z))
  if(sum(qr.resid(qr(Z),m$fatty_liver)^2)<1e-8)stop('No residual steatosis variation: ',nm)
  form<-as.formula(paste('D ~ fatty_liver + Condition(',paste(condition,collapse=' + '),')'))
  v3_log('Composition ',nm,': ',nrow(m),' people; family-conditioned dbRDA')
  fit<-vegan::dbrda(form,data=md,add='lingoes');stopifnot(fit$CCA$qrank==1L)
  setup<-fl_dv_fast_setup(fit,m,Z)
  for(scheme in c('family','family_generation')){
   block<-if(scheme=='family')as.character(m$family_id)else as.character(interaction(m$family_id,m$generation,drop=TRUE))
   P<-fl_dv_permutations(block,cfg$v3_permutations,cfg$v3_seed+match(nm,names(v3_micro_samples(x)))*10+match(scheme,c('family','family_generation')))
   null<-fl_dv_fast_null(setup,P,cfg,file.path(out,'private'),paste(nm,scheme,sep='_'))
   test<-vegan::permutest(fit,permutations=P[1:10,,drop=FALSE],model='reduced',parallel=1L);err<-max(abs(head(null,10)-as.numeric(test$F.perm)));stopifnot(err<1e-7)
   obs<-as.numeric(test$F.0);p<-(1+sum(null>=obs-sqrt(.Machine$double.eps)))/(length(null)+1)
   results[[paste(nm,scheme)]]<-data.frame(model=nm,scheme=scheme,N=nrow(m),cases=sum(m$fatty_liver),families=length(unique(m$family_id)),pseudo_F=obs,R2_total=fit$CCA$tot.chi/fit$tot.chi,R2_partial=fit$CCA$tot.chi/(fit$CCA$tot.chi+fit$CA$tot.chi),P_raw=p,permutations=length(null),status='ESTIMATED')
   sizes<-table(block);disc<-vapply(split(m$fatty_liver,block),function(v)length(unique(v))>1,logical(1))
   audit[[paste(nm,scheme)]]<-data.frame(model=nm,scheme=scheme,blocks=length(sizes),movable_people=sum(sizes[sizes>1]),discordant_blocks=sum(disc),vegan_kernel_error=err)
   tmp<-paste0(file.path(out,'private',paste0('work_kernel_',nm,'_',scheme)),c('_E.bin','_Q.bin','_P.bin','_F.csv','.log'));unlink(tmp,recursive=FALSE)
   rm(P);gc(FALSE)
  }
  saveRDS(list(model=fit,metadata=m,condition_columns=colnames(Z)),file.path(out,'private',paste0('composition_',nm,'.rds')));rm(fit,setup,D,A,Z);gc(FALSE)
 }
 result<-list(effects=ov_bind(results),diagnostics=ov_bind(audit));v3_table(result$effects,'Table5_Community_composition_EN',out);v3_table(result$diagnostics,'Composition_permutation_diagnostics_EN',out);result
}
v3_composition <- function(x,cfg,out) {
 stopifnot(as.character(packageVersion('vegan'))=='2.7.6')
 results<-audit<-list()
 for(nm in names(v3_micro_samples(x))){
  sc<-v3_micro_samples(x)[[nm]];d<-v3_design(x,sc$mask,sc$fields);m<-d$metadata;md<-d$design;md$fatty_liver<-m$fatty_liver
  A<-x$relative[,match(m$sample_id,colnames(x$relative)),drop=FALSE]
  Z<-model.matrix(reformulate(c('family_id',d$nuisance)),md);q<-qr(Z,tol=1e-9);keep<-q$pivot[seq_len(q$rank)];Z<-Z[,sort(keep),drop=FALSE];kept_nuis<-d$nuisance[d$nuisance%in%colnames(Z)];condition<-c('family_id',kept_nuis);Z<-model.matrix(reformulate(condition),md)
  stopifnot(qr(Z)$rank==ncol(Z));signature<-digest::digest(list(protocol=v3_protocol,relative=A,design=Z,outcome=m$fatty_liver,vegan='2.7.6'),algo='sha256')
  fitpath<-file.path(out,'private',paste0('composition_fit_checkpoint_',nm,'.rds'))
  if(file.exists(fitpath)){old<-readRDS(fitpath);stopifnot(identical(old$signature,signature));fit<-old$fit;setup<-old$setup;v3_log('Resume cached composition fit: ',nm)}else{
   v3_log('Composition ',nm,': fitting ',nrow(m),' people; exact vegan 2.7-6 protocol')
   D<-vegan::vegdist(t(A),method='bray');form<-as.formula(paste('D ~ fatty_liver + Condition(',paste(condition,collapse=' + '),')'));fit<-vegan::dbrda(form,data=md,add='lingoes');stopifnot(fit$CCA$qrank==1L)
   setup<-fl_dv_fast_setup(fit,m,Z);saveRDS(list(signature=signature,fit=fit,setup=setup),fitpath,compress='gzip');rm(D);gc(FALSE)
  }
  r<-m$fatty_liver-ave(m$fatty_liver,m$family_id);Qc<-setup$Q[,-1,drop=FALSE]
  compact<-fit;compact$Ybar<-setup$E;compact$pCCA$QR<-qr(Qc);compact$pCCA$rank<-ncol(Qc);compact$CCA$QR<-qr(cbind(Qc,r))
  for(scheme in c('family','family_generation')){
   cache<-file.path(out,'private',paste0('composition_result_',nm,'_',scheme,'.rds'))
   if(file.exists(cache)){old<-readRDS(cache);stopifnot(identical(old$signature,signature));results[[paste(nm,scheme)]]<-old$result;audit[[paste(nm,scheme)]]<-old$audit;v3_log('Resume composition permutation: ',nm,' ',scheme);next}
   block<-if(scheme=='family')as.character(m$family_id)else as.character(interaction(m$family_id,m$generation,drop=TRUE));P<-fl_dv_permutations(block,cfg$v3_permutations,cfg$v3_seed+match(nm,names(v3_micro_samples(x)))*10+match(scheme,c('family','family_generation')))
   prefix<-file.path(out,'private',paste0('vegan276_',nm,'_',scheme));wb<-function(v,s){con<-file(paste0(prefix,s),'wb');tryCatch(writeBin(v,con,endian='little'),finally=close(con))}
   wb(as.double(setup$E),'_E.bin');wb(as.double(Qc),'_Q.bin');wb(as.double(r),'_X.bin');wb(as.integer(P),'_P.bin')
   v3_log('Composition ',nm,' ',scheme,': ',nrow(P),' current-version permutations')
   status<-system2(cfg$python,c(shQuote(file.path(cfg$project_root,'scripts/partial_dbrda_vegan276_kernel.py')),'--prefix',shQuote(prefix),'--n',nrow(m),'--k',ncol(Qc),'--permutations',nrow(P),'--df',setup$residual_df),stdout=paste0(prefix,'.log'),stderr=paste0(prefix,'.log'));if(status!=0)stop('Current vegan kernel failed: ',prefix)
   null<-read.csv(paste0(prefix,'_F.csv'))$F_null;stopifnot(length(null)==nrow(P),all(is.finite(null)))
   # Family columns are invariant under these permutations. The compact vegan
   # object keeps the identical conditional response and restores the full df.
   v3_log('Composition ',nm,' ',scheme,': independent compact/full vegan checks')
   test<-vegan::permutest(compact,permutations=P[1:10,,drop=FALSE],model='reduced',parallel=1L);ratio<-setup$residual_df/test$df[2];reference<-as.numeric(test$F.perm)*ratio;err<-max(abs(head(null,10)-reference))
   fl_write_csv(data.frame(draw=1:10,kernel=head(null,10),compact_vegan=reference,absolute_error=abs(head(null,10)-reference)),file.path(out,'audit',paste0('composition_equivalence_',nm,'_',scheme,'.csv')))
   if(err>1e-7)stop('Current-version compact check failed: max error ',err)
   full<-vegan::permutest(fit,permutations=P[1,,drop=FALSE],model='reduced',parallel=1L);fullerr<-abs(as.numeric(full$F.perm)[1]-reference[1]);v3_log('Composition ',nm,' ',scheme,': compact error=',signif(err,3),'; original error=',signif(fullerr,3))
   if(fullerr>1e-7)stop('Compact/full original vegan check failed: ',fullerr)
   obs<-as.numeric(full$F.0);p<-(1+sum(null>=obs-sqrt(.Machine$double.eps)))/(length(null)+1)
   result<-data.frame(model=nm,scheme=scheme,N=nrow(m),cases=sum(m$fatty_liver),families=length(unique(m$family_id)),pseudo_F=obs,R2_total=fit$CCA$tot.chi/fit$tot.chi,R2_partial=fit$CCA$tot.chi/(fit$CCA$tot.chi+fit$CA$tot.chi),P_raw=p,permutations=length(null),status='ESTIMATED')
   sizes<-table(block);disc<-vapply(split(m$fatty_liver,block),function(v)length(unique(v))>1,logical(1));diag<-data.frame(model=nm,scheme=scheme,blocks=length(sizes),movable_people=sum(sizes[sizes>1]),discordant_blocks=sum(disc),vegan_kernel_error=err,full_original_vegan_error=fullerr,vegan_version='2.7.6',permutation_protocol='Reduced response and conditions co-permuted; exposure re-residualized per draw')
   results[[paste(nm,scheme)]]<-result;audit[[paste(nm,scheme)]]<-diag;saveRDS(list(signature=signature,result=result,audit=diag),cache)
   fl_write_csv(data.frame(iteration=seq_along(null),F_null=null),paste0(prefix,'_null.csv'));unlink(paste0(prefix,c('_E.bin','_Q.bin','_X.bin','_P.bin','_F.csv')),recursive=FALSE);rm(P);gc(FALSE)
  }
  rm(fit,setup,compact,A,Z);gc(FALSE)
 }
 result<-list(effects=ov_bind(results),diagnostics=ov_bind(audit));v3_table(result$effects,'Table5_Community_composition_EN',out);v3_table(result$diagnostics,'Composition_permutation_diagnostics_EN',out);result
}
v3_taxa <- function(x,cfg,out) {
 results<-diagnostics<-designs<-list();scenarios<-v3_micro_samples(x)
 for(nm in names(scenarios)){
  sc<-scenarios[[nm]];d<-v3_design(x,sc$mask,sc$fields);ids<-d$filter$feature_id[d$filter$retained];A<-d$relative[,ids,drop=FALSE]
  expected<-apply(A,2,get('LOG',asNamespace('Maaslin2')));stopifnot(identical(rownames(expected),as.character(d$metadata$sample_id)))
  replacement<-vapply(seq_len(ncol(A)),function(j)min(A[A[,j]>0,j])/2,numeric(1));v3_table(data.frame(model=nm,feature_id=ids,zero_replacement=replacement),paste0('taxa_zero_replacement_',nm),out)
  pieces<-ledger<-list();batches<-split(ids,ceiling(seq_along(ids)/20L))
  for(i in seq_along(batches)){
   features<-batches[[i]];signature<-digest::digest(list(protocol=v3_protocol,relative=A[,features,drop=FALSE],design=d$design,extract=body(ov_extract),package=as.character(packageVersion('Maaslin2'))),algo='sha256')
   folder<-file.path(cfg$project_root,'outputs/24_Two_generation_reanalysis/_checkpoints/private',nm,paste0(sprintf('%03d_',i),substr(signature,1,12)));dir.create(folder,recursive=TRUE,showWarnings=FALSE);checkpoint<-file.path(folder,'statistics.rds');start<-Sys.time()
   if(file.exists(checkpoint)){z<-readRDS(checkpoint);stopifnot(identical(z$signature,signature));piece<-z$statistics;reuse<-TRUE}
   else{
    v3_log('Taxa ',nm,' batch ',i,'/',length(batches),' (',length(features),' SGBs; N=',nrow(A),')')
    sink(file.path(folder,'console.log'));ans<-tryCatch(Maaslin2::Maaslin2(input_data=as.data.frame(A[,features,drop=FALSE]),input_metadata=d$design,output=file.path(folder,'official'),fixed_effects=c(d$bases,d$nuisance),random_effects='family_id',normalization='NONE',transform='LOG',standardize=FALSE,min_prevalence=0,min_abundance=0,min_variance=0,analysis_method='LM',cores=1,correction='BH',plot_heatmap=FALSE,plot_scatter=FALSE,save_models=TRUE),finally=sink())
    piece<-ov_extract(ans,features,d,nm,expected[,features,drop=FALSE]);saveRDS(list(signature=signature,statistics=piece),checkpoint);reuse<-FALSE;rm(ans);gc(FALSE)
   }
   pieces[[i]]<-piece;ledger[[i]]<-data.frame(model=nm,batch=i,total_batches=length(batches),seconds=as.numeric(difftime(Sys.time(),start,units='secs')),reused=reuse,completed_SGBs=sum(lengths(batches[seq_len(i)])))
   fl_write_csv(ov_bind(ledger),file.path(out,'audit',paste0('taxa_progress_',nm,'.csv')))
  }
  e<-ov_expand(ov_bind(lapply(pieces,`[[`,'effects')),d,x$candidates,nm);dg<-ov_bind(lapply(pieces,`[[`,'diagnostics'));f<-x$features[match(e$feature_id,x$features$feature_id),,drop=FALSE]
  t<-data.frame(model=nm,feature_id=e$feature_id,SGB_id=sub('^t__','',e$feature_id),species_label=ov_species_label(f$taxonomy),taxonomy=f$taxonomy,N=nrow(A),cases=sum(d$metadata$fatty_liver),families=length(unique(d$metadata$family_id)),prevalence=d$filter$prevalence[match(e$feature_id,d$filter$feature_id)],beta_log2=e$estimate,SE=e$SE,df=e$df,CI_low=e$CI_low,CI_high=e$CI_high,geometric_mean_ratio=2^e$estimate,ratio_CI_low=2^e$CI_low,ratio_CI_high=2^e$CI_high,P_raw=e$p_raw,analysis_valid=e$analysis_valid,status=e$status,singular=e$singular)
  t$nominal_raw_P_lt_005<-t$analysis_valid&!is.na(t$P_raw)&t$P_raw<.05
  results[[nm]]<-t;diagnostics[[nm]]<-dg;designs[[nm]]<-d
  v3_table(t,paste0('Table6_SGB_',nm,'_all_results_EN'),out);v3_table(t[t$nominal_raw_P_lt_005,,drop=FALSE],paste0('SGB_',nm,'_nominal_candidates_EN'),out);v3_table(d$filter,paste0('SGB_',nm,'_filter_EN'),out)
  fl_write_csv(d$mapping,file.path(out,'audit',paste0('taxa_design_',nm,'.csv')))
 }
 summary<-ov_bind(lapply(names(results),function(nm){t<-results[[nm]];data.frame(model=nm,N=t$N[1],cases=t$cases[1],families=t$families[1],candidate_SGBs=nrow(t),tested=sum(!grepl('^FILTERED',t$status)),valid=sum(t$analysis_valid),unreliable=sum(!t$analysis_valid&!grepl('^FILTERED',t$status)),nominal=sum(t$nominal_raw_P_lt_005),lower=sum(t$nominal_raw_P_lt_005&t$beta_log2<0,na.rm=TRUE),higher=sum(t$nominal_raw_P_lt_005&t$beta_log2>0,na.rm=TRUE))}))
 v3_table(summary,'SGB_model_summary_EN',out);v3_table(ov_bind(diagnostics),'SGB_model_diagnostics_EN',out)
 list(tables=results,diagnostics=diagnostics,designs=designs,summary=summary)
}
v3_stage_cache <- function(stage,compute,out) {
 path<-file.path(out,'private',paste0('stage_',stage,'.rds'))
 method<-get(paste0('v3_',stage));signature<-digest::digest(list(protocol=v3_protocol,method=body(method),B1=v3_B1,B2=v3_B2,P2=v3_P2),algo='sha256')
 if(file.exists(path)){z<-readRDS(path);if(identical(attr(z,'stage_signature'),signature)){v3_log('Resume completed stage: ',stage);return(z)}else v3_log('Stage source changed; rebuilding ',stage)}
 v3_log('START ',stage);z<-compute();attr(z,'stage_signature')<-signature;saveRDS(z,path,compress='gzip');v3_log('COMPLETE ',stage);z
}
v3_root <- function() {
 override<-Sys.getenv('FATTYLIVER_ROOT',unset='');if(nzchar(override))return(normalizePath(override,winslash='/',mustWork=TRUE))
 candidates<-c(getwd(),dirname(getwd()),'E:/YFS/workspace/FattyLiver_generation_analysis')
 for(p in candidates)if(file.exists(file.path(p,'R/reanalysis_20261024.R')))return(normalizePath(p,winslash='/',mustWork=TRUE))
 stop('Cannot locate fatty-liver project')
}
v3_main <- function(args=character(),stage='pipeline',root=v3_root()) {
 if(getRversion()<'4.5'){
  exe<-file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe')
  status<-system2(exe,c(shQuote(file.path(root,'scripts/24_two_generation_pipeline.R')),paste0('--stage=',stage),vapply(args,shQuote,character(1))))
  if(status!=0)stop('Two-generation reanalysis failed; see run log');return(invisible(NULL))
 }
 if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
 .libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),.libPaths()))
 for(f in c('config.R','R/common.R','R/entry_helpers.R','R/confounders.R','R/nature_export_approved.R','R/reanalysis_numerical_helpers.R','R/reanalysis_clinical_mi.R'))sys.source(file.path(root,f),envir=.GlobalEnv)
 cfg<-fl_config(root);cfg$output_module<-'24_Two_generation_reanalysis';cfg$run_id<-'two_generation_20261024';cfg$run_id_explicit<-FALSE;cfg$v3_seed<-20261024L;cfg$v3_permutations<-9999L
 stages<-sub('^--stage=','',grep('^--stage=',args,value=TRUE));if(length(stages))stage<-tail(stages,1)
 known<-args%in%c('--run','--check','--help','--self-test')|grepl('^--stage=|^--resume=',args)
 if(any(!known))stop('Unknown arguments: ',paste(args[!known],collapse=', '))
 if(!stage%in%c('pipeline','prepare','descriptive','clinical','clinical_mi','alpha','diversity','similarity','composition','taxa','visualization','tables'))stop('Unknown stage: ',stage)
 if('--help'%in%args){cat('24_two_generation_pipeline.R [--stage=pipeline|prepare|descriptive|clinical|clinical_mi|alpha|similarity|composition|taxa|visualization|tables] [--resume=RUN_ID] [--check]\n');return(invisible(NULL))}
 deps<-c('Maaslin2','lme4','lmerTest','geepack','detectseparation','vegan','permute','digest','ggplot2','ggrepel','patchwork','svglite','jsonlite','openxlsx','mice','miceadds','posterior')
 missing<-deps[!vapply(deps,requireNamespace,logical(1),quietly=TRUE)];if(length(missing))stop('Missing project packages: ',paste(missing,collapse=','))
 stopifnot(as.character(packageVersion('Maaslin2'))=='1.24.1')
 if('--check'%in%args){for(f in c('R/reanalysis_20261024.R','R/reanalysis_numerical_helpers.R','scripts/24_two_generation_pipeline.R'))parse(file.path(root,f));cat('PASS: syntax and all required packages; no models fitted\n');return(invisible(NULL))}
 if('--self-test'%in%args){source(file.path(root,'tests/two_generation_protocol_test.R'),local=new.env(parent=.GlobalEnv),encoding='UTF-8');return(invisible(NULL))}
 resume<-sub('^--resume=','',grep('^--resume=',args,value=TRUE))
 if(length(resume)){
  stopifnot(length(resume)==1,grepl('^[A-Za-z0-9_.-]+$',resume));module<-file.path(root,'outputs',cfg$output_module)
  possible<-file.path(module,c('_running','current'),resume);hit<-possible[dir.exists(possible)];if(length(hit)!=1)stop('Resume run not found or ambiguous')
  out<-hit;cfg$run_id<-resume;ctx<-list(cfg=cfg,module=module,run_dir=out)
 }else{ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir}
 for(p in c('tables','figures','reports','audit','private'))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
 writeLines(c('RUNNING',v3_protocol,stage),file.path(out,'run_status.txt'))
 logpath<-file.path(out,'audit','pipeline.log');logcon<-file(logpath,open='at',encoding='UTF-8');sink(logcon,split=TRUE);log_open<-TRUE;on.exit(if(log_open){sink();close(logcon)},add=TRUE)
 Sys.setenv(PYTHONPATH=paste(c(file.path(root,'.figure_qa_python'),Sys.getenv('PYTHONPATH')),collapse=.Platform$path.sep),OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4')
 v3_log('Protocol ',v3_protocol,'; run ',cfg$run_id,'; requested stage ',stage)
 writeLines(cfg$run_id,file.path(root,'outputs',cfg$output_module,'LATEST_RUN.txt'))
 x<-v3_stage_cache('prepare',function()v3_prepare(cfg,out),out)
 stopifnot(identical(unname(tools::md5sum(x$input_fingerprints$path)),x$input_fingerprints$md5))
 if(stage=='prepare')return(invisible(out))
 if(stage%in%c('pipeline','descriptive','visualization','tables'))desc<-v3_stage_cache('descriptive',function()v3_descriptive(x,cfg,out),out)
 if(stage=='descriptive')return(invisible(out))
 if(stage%in%c('pipeline','clinical','visualization','tables'))clinical<-v3_stage_cache('clinical',function()v3_clinical(x,cfg,out),out)
 if(stage=='clinical')return(invisible(out))
 if(stage%in%c('pipeline','clinical_mi','visualization','tables'))clinical_mi<-v3_stage_cache('clinical_mi',function()v3_clinical_mi(x,cfg,out),out)
 if(stage=='clinical_mi')return(invisible(out))
 if(stage%in%c('pipeline','alpha','diversity','visualization','tables'))alpha<-v3_stage_cache('alpha',function()v3_alpha(x,cfg,out),out)
 if(stage=='alpha')return(invisible(out))
 # Each following branch is independently runnable from saved expanded inputs.
 if(stage%in%c('pipeline','similarity','visualization','tables'))similarity<-v3_stage_cache('similarity',function()v3_similarity(x,cfg,out),out)
 if(stage=='similarity')return(invisible(out))
 if(stage%in%c('pipeline','composition','diversity','visualization','tables'))composition<-v3_stage_cache('composition',function()v3_composition(x,cfg,out),out)
 if(stage%in%c('composition','diversity'))return(invisible(out))
 if(stage%in%c('pipeline','taxa','visualization','tables'))taxa<-v3_stage_cache('taxa',function()v3_taxa(x,cfg,out),out)
 if(stage=='taxa')return(invisible(out))
 figures<-v3_visualize(x,clinical,alpha,similarity,composition,taxa,cfg,out)
 figures<-ov_bind(list(figures,v3_mi_visualize(clinical,clinical_mi,cfg,out)))
 v3_table(figures,'Figure_export_manifest_EN',out)
 v3_finish(x,clinical,alpha,similarity,composition,taxa,figures,cfg,out,clinical_mi)
 sink();close(logcon);log_open<-FALSE
 writeLines(c('COMPLETED',v3_protocol),file.path(out,'run_status.txt'));final<-if(basename(dirname(out))=='current')out else fl_output_complete_run(ctx)
 view<-v3_publish_view(final,cfg);v3_log('COMPLETED: ',final,'; viewing folder: ',view)
 invisible(list(authoritative=final,viewer=view))
}
v3_finish <- function(x,clinical,alpha,similarity,composition,taxa,figures,cfg,out,clinical_mi) {
 v3_log('Validating results and building English workbook')
 stopifnot(!anyDuplicated(x$metadata$sample_id),max(abs(colSums(x$relative)-1))<1e-10,all(rowSums(table(x$clinical$family_id,factor(x$clinical$generation,levels=c('G0','G1','G2')))>0)>=2),all(rowSums(table(x$metadata$family_id,x$metadata$generation)>0)>=2))
 protected_match<-match(as.character(x$metadata$participant_id),as.character(x$clinical$participant_id));stopifnot(!anyNA(protected_match))
 for(v in c('family_id','generation','sex','fatty_liver'))stopifnot(identical(as.character(x$metadata[[v]]),as.character(x$clinical[[v]][protected_match])))
 for(g in c('G0-G1','G1-G2','G0-G2')){
  z<-clinical$models[clinical$models$comparison==g,];stopifnot(length(unique(z$n[z$model%in%c('Unadjusted','C1')]))==1,length(unique(z$n[z$model%in%c('Adult_C1','C2_alcohol','C2_BMI')]))==1)
 }
 stopifnot(identical(taxa$designs$Adult_B1_same_sample$metadata$sample_id,taxa$designs$B2$metadata$sample_id),all(similarity$effects$R_python_first5_error<1e-7),all(composition$diagnostics$vegan_kernel_error<1e-7),all(composition$diagnostics$full_original_vegan_error<1e-7))
 for(stem in figures$figure)for(ext in c('pdf','png','svg','tiff'))stopifnot(file.exists(file.path(out,'figures',paste0(stem,'.',ext))))
 stopifnot(all(clinical_mi$models$valid_imputations==50L),all(is.finite(clinical_mi$diagnostics$Rhat)),all(clinical_mi$diagnostics$Rhat<1.1))
 methods<-c('# Analysis protocol and figure legends','',
 'Outcome: ultrasound-detected hepatic steatosis, RASVAMAKSA18 codes 2/3 versus code 1; codes 4/5 and missing excluded. MASLD is not classified or analysed.',
 'Scope: primary labels only (tutkno80perhe). Each cohort requires at least two observed eligible generations, not all three. Each G0-G1, G1-G2 or G0-G2 comparison uses only families with eligible members in those two roles. A third role is not required. Clinical models use the full clinical SAS source and do not require a stool sample. Microbiome analyses require valid clinical linkage, ultrasound and existing microbiome QC. Minors remain in all-age analyses. Baseline scope is defined before covariate completeness; pair comparisons additionally require both roles to remain observable in their own analysis.',
 'Primary labels do not verify direct parent-child or grandparent-grandchild relationships or current cohabitation. No perhe2/perhe3 labels are used.',
 'C1: younger age (linear, years), sex and observed older-member count when variable. One younger person per row; exposure is at least one observed older case. Unadjusted and C1 models use the same age/sex-complete sample. Family GEE, working independence, robust sandwich CI/P. C2: adults with jointly complete age, sex, alcohol and BMI; Adult C1, C2 alcohol, and C2 alcohol+BMI are fitted to the same people. No older mean age/female proportion, geography, smoking or diet is added to clinical models.',
 'B1: generation, age, sex, preparation batch, sequencing pool, sampling season and antibiotic use in the preceding three months. Age retains the existing within-generation natural splines (3 df per generation); season retains sine/cosine of stool month. Batch/pool are factors. Exactly redundant nuisance columns are removed by QR, never by P value. No survey region/urbanicity, matrix total, freezing delay, stool-ultrasound gap, bowel symptoms or DNA concentration is in this adopted adjustment set.',
 'B2: B1 plus current smoking/tobacco use (binary), alcohol servings/day, diet score 0-27, adult physical-activity index PAI18, log10 energy intake ENERCK18, and continuous adult BMI. Adult B1 and B2 use exactly the same complete sample. Unmatched/missing nutrient metadata are treated as missing, not as zero. B2 is a complete-case adult association. Old restricted-cohort MI files are not reused.',
 'Clinical C2 multiple-imputation supplement: missing adult alcohol/BMI are newly imputed for each expanded comparison, with 50 datasets and predictive mean matching (five donors). Family random-intercept two-level PMM is used for repeated younger members; ordinary PMM is used when all younger members have distinct primary labels. Predictors include age, sex, observed disease, older-case exposure, older-member count, log1p alcohol and BMI. Observed cells, disease, age, sex and labels are protected. GEE log-ORs and robust variances are pooled using Rubin rules. Diagnostics and all per-imputation fits are retained. This supplement does not impute microbial B2 covariates. The new participant-level RDS is private/adult_clinical_imputed_two_generation.rds; original data are preserved.',
 'Shannon: all nonzero SGBs, using counts divided by the full sample count total, natural-log entropy. Family-cluster GEE with B1/B2. Overall and generation-specific disease coefficients are estimated; generation heterogeneity uses a joint Wald test. These are adjusted mean Shannon differences.',
 'Community composition: Bray-Curtis on full-denominator relative abundance; partial dbRDA with Lingoes correction, conditioning on family and B1/B2 covariates. The estimand is a within-family conditional composition association, different from the marginal Shannon or taxon mixed model. Reduced-model permutations stay within family; a stricter family-by-generation scheme is also reported. 9,999 permutations follow vegan 2.7.6 semantics: reduced response and conditioning variables are co-permuted and the fixed exposure is re-residualized at every draw. Each optimized kernel is checked against ten compact equivalent vegan draws and one original full-design vegan draw. Total corrected R2 and partial R2 are both reported.',
 'P1: mean Bray-Curtis IF versus all BF cross-generation dyads among the same pair-eligible nodes, no regression covariates. Each two-person pair counted once; people can occur in many pairs. Unrestricted node permutations separately within older/younger generation roles preserve repeated-node identities; dyads are never treated as independent observations for P values.',
 'P2: both endpoints age (natural splines, 3 df), sex, preparation batch, sequencing pool, seasonal sine/cosine, recent antibiotics. No same-region/batch/pool indicators are added. Residual MRQAP with role-specific sex-by-five-year node permutation restrictions, a continuation of the prior exchangeability sensitivity. Exact additive endpoint projections avoid a large dyadic design matrix; first five null statistics independently verified in R. Two-sided studentized permutation P, 9,999 draws. P1 and P2 have separately reported sample counts. No dyad-independent CI is reported.',
 'Taxa: candidate union rebuilt as prevalence >=10% in any baseline generation; each model then retains candidates with overall prevalence >=10% in its actual sample. Filtered rows and nonconverged models remain visible. Official MaAsLin2 1.24.1, normalization NONE (already relative), transformation LOG, standardize FALSE, family random intercept. Zeros are replaced by half the smallest positive abundance for that SGB in that model sample, then log2; positive values are unchanged. No rarefaction, batch correction, retained-feature reclosure or outcome SD scaling. Beta is the adjusted log2 geometric-mean ratio of zero-replaced relative abundance; 2^beta is not an absolute microbial-load ratio.',
 'Taxon inference uses official MaAsLin2 coefficients and Satterthwaite inference checked against lmerTest. Model diagnostics include convergence, optimizer Hessian, singular random effects, residual skewness and family ICC. Nonconverged fits are not called valid. Singular but otherwise converged fits remain flagged. Official batch BH outputs are private and not reported as whole-study FDR.',
 'All reported P values are raw at the user request. Raw P<0.05 taxon candidates are exploratory, without multiplicity control. Representative species are selected by smallest valid raw P among species detected in at least ten people in each disease group, up to eight; all zeros and observed values are included in descriptive box plots.',
 'Figures: boxes median/IQR, whiskers 1.5 IQR; participant dots grey, box outlines black. Clinical and Shannon forest intervals are 95% robust CIs. Similarity bars are means or adjusted coefficients, not independent sample estimates. Composition bars show total corrected R2, without invented uncertainty intervals. Volcano points are valid SGB models, vertical zero and horizontal raw P=0.05; blue lower and orange higher. Binary reference blue #79CAFB, comparison orange #F7A24F. Generation G0 #F1B656, G1 #397FC7, G2 #040676. All figures exported PDF, SVG, PNG and LZW TIFF at 600 dpi.',
 'References: Long et al. Liver Int. 2019;39:740-747. doi:10.1111/liv.13956. Gacesa et al. Nature. 2022;604:732-739. doi:10.1038/s41586-022-04567-7. Ma et al. Nat Commun. 2024;15:3612. doi:10.1038/s41467-024-47859-4. Valles-Colomer et al. Nat Microbiol. 2022;7:87-96. doi:10.1038/s41564-021-01021-8. Dekker et al. Psychometrika. 2007;72:563-581. doi:10.1007/s11336-007-9016-1. Legendre & Anderson. Ecol Monogr. 1999;69:1-24. Mallick et al. PLoS Comput Biol. 2021;17:e1009442. doi:10.1371/journal.pcbi.1009442. White et al. Stat Med. 2011;30:377-399. doi:10.1002/sim.4067.',
 'Study adaptations (not verbatim published protocols): all-age models with minors, primary-label exposure and observation count, within-generation age basis, pool effects, role-specific node permutations, rectangular additive-endpoint MRQAP, family-conditioned restricted community permutations. This revised protocol was specified after earlier results were known; it is not a preregistration.')
 writeLines(methods,file.path(out,'reports/Methods_and_figure_legends_EN.md'),useBytes=TRUE)
 legends<-c('# Figure legends','',
  paste0('Figure 1. Study populations and generation structure. The clinical cohort includes ',nrow(x$clinical),' participants in ',length(unique(x$clinical$family_id)),' primary families. The microbiome cohort includes ',nrow(x$metadata),' participants in ',length(unique(x$metadata$family_id)),' primary families. At least two observed eligible generations are required; a third generation is not required. Minors are retained. Boxes show medians and interquartile ranges, with 1.5-IQR whiskers; grey points are participants. Bars show observed ultrasound-steatosis prevalence. No regression adjustment is applied.'),'',
  'Figure 2. Cross-generation association of ultrasound hepatic steatosis. Panel a shows observed younger-member prevalence according to whether at least one observed older member in the same primary family has steatosis; labels show cases/participants. Panel b shows unadjusted and C1-adjusted odds ratios with robust 95% confidence intervals, fitted to the same age/sex-complete sample in each comparison. C1 includes younger age, sex and older-member observation count when variable. All models use family-cluster GEE. Primary labels do not verify direct parent-child relationships.','',
  'Figure 3. Within-family and between-family gut microbiome similarity. Panel a shows mean Bray-Curtis distances for all IF and BF cross-generation pairs among the same eligible nodes (P1). Panel b shows adjusted IF-minus-BF distance (P2). Lower distance indicates more similar composition. P2 includes both endpoints age, sex, preparation batch, sequencing pool, sampling season and recent antibiotics. Node permutations preserve repeated participation in pairs; dyads are not independent sample sizes. No independent-dyad confidence intervals are shown.','',
  paste0('Figure 4. Hepatic steatosis and Shannon diversity. Observed distributions and B1-adjusted generation-specific mean differences among ',taxa$summary$N[taxa$summary$model=='B1'],' participants (',taxa$summary$cases[taxa$summary$model=='B1'],' cases). Differences are steatosis minus no steatosis; horizontal intervals are robust 95% confidence intervals. B1 includes age, sex, generation, preparation batch, sequencing pool, sampling season and recent antibiotics; GEE accounts for primary-family clustering. Age is encoded using generation-specific splines and season using sine/cosine terms.'),'',
  'Figure 5. Conditional association of steatosis with overall microbial composition. Bars show the percentage of total Lingoes-corrected Bray-Curtis variation attributable to steatosis in partial dbRDA, conditioning on primary family and B1/B2 covariates. Adult B1 and B2 use identical complete cases. This is a within-family conditional association. Permutation P values and partial R2 are provided in Table5_Community_composition_EN; no confidence interval is inferred from pair counts.','',
  'Figure 6. Overall species-level abundance associations: all-age B1 (a) and adult B2 (b). Panels use separate axis ranges. Each point is a valid species-level genome bin (SGB) mixed model. The horizontal axis shows the adjusted log2 ratio for zero-replaced relative abundance (steatosis/no steatosis); the vertical axis is -log10(raw P). Official MaAsLin2 LOG models include a primary-family random intercept and the indicated B1/B2 adjustment set. For each SGB, zeros are replaced by half its minimum positive value in that model sample; positive values are unchanged. Blue/orange indicate lower/higher abundance at raw P<0.05; grey indicates P>=0.05. Raw-P candidates are exploratory and do not represent multiplicity-controlled discoveries.','',
  'Figure 7 (B1 and B2). Relative-abundance distributions for representative SGBs. Up to eight valid SGB models with the smallest raw P values are shown, requiring detection in at least ten participants in each disease group. Distributions include zeros and all observations. Boxes show medians and interquartile ranges, whiskers extend to 1.5 IQR, and grey points represent participants. Panels have separate abundance scales. These distributions are unadjusted; adjusted effects and model diagnostics are reported in the corresponding numerical tables. SGB IDs are retained where species labels are unresolved or ambiguous.','',
  'Supplementary Figure 1. Adult clinical C2 extensions. Adult C1, C1 plus alcohol, and C1 plus alcohol and continuous BMI are fitted to the same complete-case sample within each comparison. Points and intervals are family-GEE odds ratios and robust 95% confidence intervals. Low-event diagnostic flags remain in the full clinical table.','',
  'Supplementary Figure 2. Adult Shannon associations before and after B2 adjustment. Adult B1 and B2 use exactly the same complete-case participants. B2 adds current smoking/tobacco use, alcohol, diet score, physical activity, total energy and continuous BMI to B1. Points are adjusted mean differences with robust 95% confidence intervals.','',
  'Supplementary Figure 3. Adult clinical C2 BMI results in complete cases and after multiple imputation. Both methods adjust for younger age, sex, observed older-member count when variable, alcohol and BMI, with family-cluster GEE. Fifty new imputed datasets target missing adult alcohol/BMI within each expanded comparison, using family-aware predictive mean matching and Rubin pooling. In this run only alcohol is missing; BMI is observed. Complete-case and imputation samples differ in size and case counts; differences must not be interpreted solely as a covariate-adjustment effect. Original data are preserved.','',
  'All reported P values are raw. IF, within family; BF, between families; B1/B2/C1/C2, adjustment sets defined in Adjustment_sets_and_literature_EN. The outcome is ultrasound steatosis; MASLD is not analysed.')
 writeLines(legends,file.path(out,'reports','Figure_legends_EN.md'),useBytes=TRUE)
 summary<-c('# 新调整集及两代可用范围：结果总结','',paste0('临床：',nrow(x$clinical),'人、',length(unique(x$clinical$family_id)),'个主家系；菌群：',nrow(x$metadata),'人、',length(unique(x$metadata$family_id)),'个主家系。'),'','## C1：临床跨代关联',capture.output(print(clinical$models[clinical$models$model=='C1',c('comparison','n','cases','OR','CI_low','CI_high','p_value','status')],row.names=FALSE)),'','## C2：成人共同样本',capture.output(print(clinical$models[clinical$models$model%in%c('Adult_C1','C2_alcohol','C2_BMI'),c('comparison','model','n','cases','OR','CI_low','CI_high','p_value','status')],row.names=FALSE)),'','## 家系内外菌群相似性',capture.output(print(similarity$effects[,c('comparison','model','N','IF_mean','BF_mean','IF_minus_BF','P_raw')],row.names=FALSE)),'','## Shannon',capture.output(print(alpha$effects[,c('model','generation','N','cases','beta','CI_low','CI_high','P_raw')],row.names=FALSE)),'','## 整体组成',capture.output(print(composition$effects,row.names=FALSE)),'','## SGB关联',capture.output(print(taxa$summary,row.names=FALSE)),'','本轮均为超声脂肪肝，不分析MASLD。仅使用主家系；家系标签不证明直接亲子或同住。所有P为原始P；菌种P<0.05仅作为探索性候选。成人结果来自新的共同完整病例，未沿用旧插补数据。','原始SAS/RDS及既有图表保留。新数据保存在本轮private/expanded_analysis_data.rds；查看文件夹不包含逐人数据。')
 writeLines(enc2utf8(summary),file.path(out,'reports/Results_summary_zh.md'),useBytes=TRUE)
 mi_summary<-c('','## C2：新范围成人多重插补（补充分析）',capture.output(print(clinical_mi$models,row.names=FALSE)),'',capture.output(print(clinical_mi$flow,row.names=FALSE)),'','临床C2按饮酒/BMI缺失设定重新插补50组；本轮实际仅饮酒存在缺失，BMI均已观察，无需插补。原始数据保留；菌群B2仍为共同完整病例。完整病例与插补结果的样本范围不同，不把估计差异全部归因于增加协变量。','菌种模型报告相对丰度的对数尺度关联；较多菌的负系数不等于绝对菌量普遍降低。未经多重检验校正的候选不作为确证性发现。')
 summary<-gsub('成人结果来自新的共同完整病例，未沿用旧插补数据。','成人完整病例结果与新范围临床C2插补补充结果分开报告，未沿用旧插补数据。',summary,fixed=TRUE)
 writeLines(enc2utf8(c(summary,mi_summary)),file.path(out,'reports/Results_summary_zh.md'),useBytes=TRUE)
 # Compact manuscript-facing tables accompany complete numerical results.
 pretty_p<-function(p)ifelse(is.na(p),'',ifelse(p<.0001,format(p,scientific=TRUE,digits=3),sprintf('%.4f',p)))
 pretty_or<-function(z)ifelse(is.finite(z$OR),sprintf('%.2f (%.2f-%.2f)',z$OR,z$CI_low,z$CI_high),paste0('NE: ',z$status))
 z<-clinical$models[clinical$models$model%in%c('Unadjusted','C1'),];display<-data.frame(Comparison=z$comparison,Model=ifelse(z$model=='C1','C1 adjusted','Unadjusted'),Participants=z$n,Cases=z$cases,`OR (95% CI)`=pretty_or(z),`P value`=pretty_p(z$p_value),check.names=FALSE)
 v3_table(display,'Main_Table2_Clinical_associations_EN',out)
 z<-clinical$models[clinical$models$model%in%c('Adult_C1','C2_alcohol','C2_BMI'),];display<-data.frame(Comparison=z$comparison,Model=z$model,Participants=z$n,Cases=z$cases,`OR (95% CI)`=pretty_or(z),`P value`=pretty_p(z$p_value),Status=z$status,check.names=FALSE)
 v3_table(display,'Supplementary_Table1_Adult_C2_EN',out)
 z<-alpha$effects;display<-data.frame(Model=z$model,Generation=z$generation,Participants=z$N,Cases=z$cases,`Shannon difference (95% CI)`=sprintf('%.3f (%.3f to %.3f)',z$beta,z$CI_low,z$CI_high),`P value`=pretty_p(z$P_raw),check.names=FALSE)
 v3_table(display,'Main_Table4_Shannon_effects_EN',out)
 adjustment_table<-data.frame(
  Analysis=c('Descriptive','C1: all ages','C2: adults','B1: individual microbiome','B2: adult microbiome','P1: raw similarity','P2: adjusted similarity','Clinical C2 missing-data supplement'),
  Variables=c('No regression covariates','Younger age; sex; observed older-member count if variable; family GEE','C1 + alcohol; then + continuous BMI; same adult sample','Age; sex; pooled generation; preparation batch; sequencing pool; stool season; recent antibiotics; family dependence','B1 + smoking; alcohol; diet score; physical activity; total energy; continuous BMI; common adult complete cases','No regression covariates; role-specific node permutations','Both endpoint age; sex; preparation batch; sequencing pool; season; recent antibiotics; residual node permutations','Impute missing adult alcohol/BMI only; observed age/sex/exposure/outcome fixed; family-aware PMM; 50 datasets; pooled GEE'),
  Literature=c('Long 2019 Table 1; Ma 2024 Fig 2b','Long 2019 Methods 2.5 and Table 2','Long 2019 Table 2','Gacesa 2022 Methods; Ma 2024 Fig 2a','Ma 2024 Fig 2a','Gacesa 2022 Extended Data Fig 4; Valles-Colomer 2022','Gacesa 2022; Ma 2024; Dekker 2007','White 2011, doi:10.1002/sim.4067'),
  Adaptation=c('Observed eligible cohort','Primary labels; all ages; observation-count adjustment','Current measured alcohol and BMI','Generation-specific age spline; sequencing pool; family structure','Current measured lifestyle indices','Cross-generation roles; repeated-node inference','Endpoint formula and restricted rectangular MRQAP','Two-level PMM for repeated family members; ordinary PMM for singleton families'),check.names=FALSE)
 v3_table(adjustment_table,'Adjustment_sets_and_literature_EN',out)
 file.copy(file.path(cfg$project_root,'analysis_plan/current_two_generation_20261024.md'),file.path(out,'reports','Adopted_analysis_protocol_EN.md'),overwrite=TRUE)
 md_table<-function(d){esc<-function(a)gsub('|','/',as.character(a),fixed=TRUE);c(paste0('| ',paste(esc(names(d)),collapse=' | '),' |'),paste0('| ',paste(rep('---',ncol(d)),collapse=' | '),' |'),vapply(seq_len(nrow(d)),function(i)paste0('| ',paste(esc(unlist(d[i,,drop=FALSE])),collapse=' | '),' |'),character(1)))}
 z<-clinical$models[clinical$models$model=='C1',];c1_tab<-data.frame(比较=z$comparison,N=z$n,病例=z$cases,`OR（95%CI）`=pretty_or(z),原始P=pretty_p(z$p_value),check.names=FALSE)
 z<-clinical$models[clinical$models$model=='C2_BMI',];cc_tab<-data.frame(比较=z$comparison,N=z$n,病例=z$cases,`OR（95%CI）`=pretty_or(z),原始P=pretty_p(z$p_value),状态=z$status,check.names=FALSE)
 z<-clinical_mi$models[clinical_mi$models$model=='C2_BMI_MI',];mi_tab<-data.frame(比较=z$comparison,N=z$N,病例=z$cases,`OR（95%CI）`=pretty_or(z),原始P=pretty_p(z$P_raw),状态=z$status,check.names=FALSE)
 z<-similarity$effects;sim_tab<-data.frame(比较=z$comparison,模型=z$model,N=z$N,IF均值=sprintf('%.4f',z$IF_mean),BF均值=sprintf('%.4f',z$BF_mean),`IF−BF`=sprintf('%.4f',z$IF_minus_BF),原始P=pretty_p(z$P_raw),check.names=FALSE)
 z<-alpha$effects[alpha$effects$model=='B1',];alpha_tab<-data.frame(代际=z$generation,N=z$N,病例=z$cases,`Shannon差值（95%CI）`=sprintf('%.3f (%.3f to %.3f)',z$beta,z$CI_low,z$CI_high),原始P=pretty_p(z$P_raw),check.names=FALSE)
 z<-composition$effects;comp_tab<-data.frame(模型=z$model,置换=z$scheme,N=z$N,`总变异解释度（%）`=sprintf('%.4f',100*z$R2_total),原始P=pretty_p(z$P_raw),check.names=FALSE)
 zh<-c('# 两代可用范围与新调整集：结果总结','',paste0('临床：',nrow(x$clinical),'人、',length(unique(x$clinical$family_id)),'个主家系。菌群：',nrow(x$metadata),'人、',length(unique(x$metadata$family_id)),'个主家系；B1完整样本',taxa$summary$N[taxa$summary$model=='B1'],'人；成人B1/B2共同完整样本',taxa$summary$N[taxa$summary$model=='B2'],'人。'),'','所有结果均为超声脂肪肝，不分析MASLD；全龄分析保留未成年人。家系只要求相应两代可用，不要求第三代参加。临床分析不要求有菌群样本。主家系标签不证明直接亲子关系或同住。','',
 '## 1. 临床跨代：C1全龄主模型','',md_table(c1_tab),'','C1调整年轻成员年龄、性别；年长成员观察人数变化时加入人数，GEE处理家系聚类。粗模型与C1在相同样本中比较。','',
 '## 2. 成人饮酒＋BMI扩展','',md_table(cc_tab),'','上表为共同完整病例。下表为新范围50组多重插补的补充结果：','',md_table(mi_tab),'','本轮实际仅饮酒存在缺失，BMI均已观察。插补不改变年龄、性别、脂肪肝、暴露、家系或其他已观察值。插补与完整病例样本量不同，估计差异不能全部归因于增加协变量。稀疏事件警告保留。','',
 '## 3. 家系菌群相似性','',md_table(sim_tab),'','P1为无回归调整的原始比较；P2调整双方年龄、性别、建库批次、测序池、季节及近期抗生素。P2的IF−BF为调整系数，IF/BF均值仍为原始均值。节点置换处理同一人重复配对；配对数不能视为独立样本量。9999次置换的最小P值分辨率为0.0001。','',
 '## 4. Shannon多样性','',md_table(alpha_tab),'',paste0('B1代际差异联合检验原始P=',pretty_p(alpha$interaction$P_raw[alpha$interaction$model=='B1']),'。差值为脂肪肝减去无脂肪肝。'),'','B1调整年龄、性别、代际、建库批次、测序池、季节及近期抗生素，GEE处理家系聚类。成人B2在相同完整病例中增加吸烟、饮酒、饮食、活动、总能量和BMI；详细效应见英文表格。','',
 '## 5. 整体菌群组成','',md_table(comp_tab),'','这是控制家系及B1/B2后的家系内条件关联；分别报告家系内及家系内同代置换。解释度为占校正后总组成变异的比例，不能解释为临床预测准确度。','',
 '## 6. 总体菌种关联','',md_table(taxa$summary),'','官方MaAsLin2对数相对丰度混合模型，主家系随机截距；B1/B2按新调整集重新计算。nominal/lower/higher为原始P<0.05的探索性候选总数及方向，不是FDR控制后的发现。较多负系数不等于绝对菌量普遍降低。代表性原始丰度分布与完整调整效应、诊断分开报告。','',
 '## 本轮改变','', '- 从三代均参加改为相应两代可用，并从完整临床SAS重新构建临床样本。','- 所有现行编号入口改为C1/C2/B1/B2/P1/P2；不沿用旧范围的统计结果或插补数据。','- 新范围临床成人饮酒缺失重新插补50组；菌群B2仍使用共同完整病例。','- 输出英文工作簿、图注、调整集来源表；12张图均提供PDF、PNG、TIFF和SVG。','- 原始SAS/RDS及既有结果保留。逐人新数据在权威结果private/expanded_analysis_data.rds；新插补数据在private/adult_clinical_imputed_two_generation.rds；查看目录不含逐人数据。')
 writeLines(enc2utf8(zh),file.path(out,'reports','Results_summary_zh.md'),useBytes=TRUE)
 files<-list.files(file.path(out,'tables'),pattern='[.]csv$',full.names=TRUE);wb<-openxlsx::createWorkbook();index<-data.frame(sheet=character(),source=character());used<-character()
 for(i in seq_along(files)){
  sheet<-substr(gsub('[^A-Za-z0-9_]','_',tools::file_path_sans_ext(basename(files[i]))),1,27);sheet<-paste0(sprintf('%02d_',i),sheet);sheet<-substr(sheet,1,31)
  openxlsx::addWorksheet(wb,sheet);d<-read.csv(files[i],fileEncoding='UTF-8-BOM',check.names=FALSE,stringsAsFactors=FALSE);openxlsx::writeData(wb,sheet,d,headerStyle=openxlsx::createStyle(fontName='Times New Roman',fontSize=11,textDecoration='bold',border=c('top','bottom'),wrapText=TRUE),withFilter=FALSE)
  openxlsx::freezePane(wb,sheet,firstRow=TRUE);openxlsx::setColWidths(wb,sheet,cols=seq_len(ncol(d)),widths=pmin(42,pmax(12,nchar(names(d))+2)))
  pcols<-which(grepl('^P_raw$|^p_value$|^p_raw$',names(d)));if(length(pcols)&&nrow(d))openxlsx::addStyle(wb,sheet,openxlsx::createStyle(numFmt='0.00E+00'),rows=2:(nrow(d)+1),cols=pcols,gridExpand=TRUE)
  index<-rbind(index,data.frame(sheet=sheet,source=basename(files[i])))
 }
 openxlsx::addWorksheet(wb,'Contents');openxlsx::writeData(wb,'Contents',index);openxlsx::setColWidths(wb,'Contents',1:2,c(33,75))
 priority<-match(c('Table1_Clinical_by_generation_EN.csv','Table1_Microbiome_by_generation_EN.csv','Main_Table2_Clinical_associations_EN.csv','Supplementary_Table1_Adult_C2_EN.csv','Supplementary_Clinical_C2_multiple_imputation_EN.csv','Table3_Family_microbiome_similarity_EN.csv','Main_Table4_Shannon_effects_EN.csv','Table5_Community_composition_EN.csv','SGB_model_summary_EN.csv','Adjustment_sets_and_literature_EN.csv'),index$source);priority<-priority[!is.na(priority)]
 openxlsx::worksheetOrder(wb)<-c(nrow(index)+1L,priority,setdiff(seq_len(nrow(index)),priority))
 openxlsx::saveWorkbook(wb,file.path(out,'tables','Hepatic_steatosis_two_generation_tables_EN.xlsx'),overwrite=TRUE)
 # Source, protocol and session provenance are retained outside the public viewer.
 for(f in c('R/reanalysis_20261024.R','R/reanalysis_numerical_helpers.R','R/reanalysis_clinical_mi.R','scripts/24_two_generation_pipeline.R','scripts/two_generation_pair_kernel.py','scripts/partial_dbrda_vegan276_kernel.py'))file.copy(file.path(cfg$project_root,f),file.path(out,'audit',basename(f)),overwrite=TRUE)
 plot_helpers<-c('v3_style','v3_visualize','v3_mi_visualize','nf_theme','nf_skill_root','nf_alignment','nf_export','nf_qa_command','ov_export')
 plot_source<-file.path(out,'audit','plot_and_export_source.R');writeLines(unlist(lapply(plot_helpers,function(n)c(paste0(n,' <- '),deparse(get(n),width.cutoff=120L),''))),plot_source)
 invisible(parse(plot_source));nf_qa_command(cfg,'validate_figure.py',c(shQuote(plot_source),'--json'),file.path(out,'audit','plot_source_preflight.json'))
 capture.output(sessionInfo(),file=file.path(out,'audit','sessionInfo.txt'))
 stopifnot(identical(unname(tools::md5sum(x$input_fingerprints$path)),x$input_fingerprints$md5))
 writeLines(c('PASS: protected inputs unchanged; pair-specific two-generation scope; clinical and adult same-sample checks; QR design checks; permutation-kernel equivalence; all four figure formats.','No MASLD. Raw P exploratory reporting. Diagnostic flags remain in full tables.'),file.path(out,'audit','final_validation.txt'))
 invisible(TRUE)
}
v3_publish_view <- function(source,cfg) {
 if(file.exists(file.path(cfg$project_root,'config/results_viewer.json'))){
  gallery<-new.env(parent=baseenv());sys.source(file.path(cfg$project_root,'R/result_gallery.R'),envir=gallery)
  return(gallery$fl_publish_organized_results(cfg))
 }
 base<-file.path(cfg$project_root,'Results_表格与图片','Reanalysis_two_generation_20261024');target<-file.path(base,basename(source));dir.create(target,recursive=TRUE,showWarnings=FALSE)
 files<-unlist(lapply(c('tables','figures','reports'),function(p)list.files(file.path(source,p),recursive=TRUE,full.names=TRUE)),use.names=FALSE)
 rel<-substring(gsub('\\','/',files,fixed=TRUE),nchar(gsub('\\','/',source,fixed=TRUE))+2L)
 for(i in seq_along(files)){dest<-file.path(target,rel[i]);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(files[i],dest,overwrite=TRUE),identical(unname(tools::md5sum(files[i])),unname(tools::md5sum(dest))))}
 fl_write_csv(data.frame(file=rel,source_md5=unname(tools::md5sum(files))),file.path(target,'file_manifest.csv'))
 html<-c('<!doctype html><meta charset="utf-8"><style>body{font:15px Arial;margin:30px;line-height:1.6}img{width:720px;max-width:95%}li{margin:5px}a{color:#255a80}</style><h1>Hepatic steatosis: two-generation reanalysis</h1><p>C1 / C2 / B1 / B2 / P1 / P2. Ultrasound steatosis only.</p><p><a href="tables/Hepatic_steatosis_two_generation_tables_EN.xlsx">English tables workbook</a> | <a href="reports/Results_summary_zh.md">Results summary</a> | <a href="reports/Methods_and_figure_legends_EN.md">Methods and legends</a></p>')
 for(folder in c('figures','tables','reports')){r<-rel[startsWith(rel,paste0(folder,'/'))];html<-c(html,paste0('<h2>',folder,'</h2><ul>'));for(f in r){url<-utils::URLencode(f,reserved=FALSE);html<-c(html,paste0('<li><a href="',url,'">',basename(f),'</a></li>'));if(grepl('[.]png$',f))html<-c(html,paste0('<img src="',url,'">'))};html<-c(html,'</ul>')}
 writeLines(html,file.path(target,'index.html'),useBytes=TRUE)
 writeLines(c('<!doctype html><meta charset="utf-8"><h1>Two-generation reanalysis</h1>',paste0('<p><a href="',basename(source),'/index.html">Open current results: ',basename(source),'</a></p>')),file.path(base,'index.html'),useBytes=TRUE)
 writeLines(basename(source),file.path(base,'CURRENT.txt'));target
}
v3_style <- function() nf_theme()+ggplot2::theme(axis.text=ggplot2::element_text(size=6.5,colour='black'),axis.title=ggplot2::element_text(size=7),legend.text=ggplot2::element_text(size=6.5),plot.title=ggplot2::element_text(size=7,margin=ggplot2::margin(b=5)),strip.text=ggplot2::element_text(size=6.5),strip.background=ggplot2::element_blank(),plot.margin=ggplot2::margin(6,7,6,7))
v3_visualize <- function(x,clinical,alpha,similarity,composition,taxa,cfg,out,only=NULL) {
 v3_log('Rendering new figures from new result tables only')
 records<-list();binary<-setNames(unname(cfg$binary_colors),c('No steatosis','Steatosis'));familypal<-setNames(unname(cfg$binary_colors),c('BF','IF'))
 export<-function(p,name,width=183,height=115,multi=FALSE,panels=1L){
  if(!is.null(only)&&!name%in%only)return(invisible(NULL))
  if(multi){r<-nf_export(list(name=name,plot=p,width=width,height=height,panels=panels,export_tiff=TRUE),out,cfg)}else r<-ov_export(p,name,out,cfg,width,height)
  records[[name]]<<-data.frame(figure=name,width_mm=width,height_mm=height,formats='PDF;PNG;SVG;TIFF',dpi=600);invisible(NULL)
 }
 patch<-function(plots,ncol=2){plots<-lapply(seq_along(plots),function(i)plots[[i]]+v3_style()+ggplot2::labs(tag=letters[i]));patchwork::wrap_plots(plots,ncol=ncol)}
 # Flow is deliberately a simple, independent diagram: two different eligible cohorts.
 flow<-data.frame(x=c(0,0,0,2,2,2),y=rep(3:1,2),label=c(paste0('Clinical records\nn = ',x$flow$N[1]),paste0('Valid ultrasound and generation\nn = ',x$flow$N[2]),paste0('At least two generations\nn = ',nrow(x$clinical),'; ',length(unique(x$clinical$family_id)),' families'),paste0('Microbiome + valid ultrasound\nn = ',x$flow$N[4]),paste0('At least two generations\nn = ',nrow(x$metadata),'; ',length(unique(x$metadata$family_id)),' families'),paste0('B1: n = ',x$flow$N[6],'\nAdult common B2: n = ',x$flow$N[7])))
 p<-ggplot2::ggplot(flow,ggplot2::aes(x,y))+ggplot2::geom_rect(ggplot2::aes(xmin=x-.85,xmax=x+.85,ymin=y-.28,ymax=y+.28),fill='white',colour='black',linewidth=.3)+ggplot2::geom_text(ggplot2::aes(label=label),family='Arial',size=6.5/(72.27/25.4),lineheight=1.15)+ggplot2::geom_segment(data=data.frame(x=rep(c(0,2),each=2),y=rep(c(2.7,1.7),2)),ggplot2::aes(x=x,xend=x,y=y,yend=y-.4),arrow=grid::arrow(length=grid::unit(1.5,'mm'),type='closed'),linewidth=.3,inherit.aes=FALSE)+ggplot2::scale_x_continuous(limits=c(-1,3))+ggplot2::scale_y_continuous(limits=c(.5,3.5))+ggplot2::theme_void(base_family='Arial')
 export(p,'Figure1a_Cohort_flow',height=95)
 m<-x$clinical
 p1<-ggplot2::ggplot(m,ggplot2::aes(generation,age,fill=generation))+ggplot2::geom_boxplot(colour='black',linewidth=.3,width=.55,outlier.shape=NA)+ggplot2::geom_point(position=ggplot2::position_jitter(width=.12,height=0,seed=20261024),colour='#888888',size=.25,alpha=.20,show.legend=FALSE)+ggplot2::scale_fill_manual(values=cfg$generation_colors)+ggplot2::labs(x=NULL,y='Age (years)',title='Age distribution')+ggplot2::guides(fill='none')
 prev<-v3_counts(m,'Clinical');prev<-prev[prev$group%in%c('G0','G1','G2'),];prev$group<-factor(prev$group,levels=c('G0','G1','G2'))
 p2<-ggplot2::ggplot(prev,ggplot2::aes(group,prevalence_percent,fill=group))+ggplot2::geom_col(width=.55,colour='black',linewidth=.3)+ggplot2::geom_text(ggplot2::aes(label=paste0(cases,'/',N)),vjust=-.7,size=6/(72.27/25.4),family='Arial')+ggplot2::scale_fill_manual(values=cfg$generation_colors)+ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.15)))+ggplot2::labs(x=NULL,y='Hepatic steatosis (%)',title='Observed prevalence')+ggplot2::guides(fill='none')
 export(patch(list(p1,p2)),'Figure1b_Age_and_prevalence',height=100,multi=TRUE,panels=2)
 # Clinical prevalence uses exactly the C1 analysis sample, not all descriptive rows.
 observed<-ov_bind(lapply(v3_pairs,function(g){z<-readRDS(file.path(out,'private',paste0('clinical_pairs_',paste(g,collapse='-'),'.rds')));z<-z[v3_cc(z,c('age','sex')),];do.call(rbind,lapply(0:1,function(e){a<-z[z$older_any_case==e,];data.frame(comparison=paste(g,collapse='-'),exposure=if(e)'Older case'else'No older case',N=nrow(a),cases=sum(a$fatty_liver),percent=100*mean(a$fatty_liver))}))}))
 observed$comparison<-factor(observed$comparison,levels=c('G0-G1','G1-G2','G0-G2'));observed$exposure<-factor(observed$exposure,levels=c('No older case','Older case'));cpal<-setNames(unname(cfg$binary_colors),levels(observed$exposure));v3_table(observed,'Figure2_prevalence_source_EN',out)
 p1<-ggplot2::ggplot(observed,ggplot2::aes(comparison,percent,fill=exposure))+ggplot2::geom_col(position=ggplot2::position_dodge(.75),width=.65,colour='black',linewidth=.25)+ggplot2::geom_text(ggplot2::aes(label=sprintf('%.1f%%\n%d/%d',percent,cases,N)),position=ggplot2::position_dodge(.75),vjust=-.35,size=5.5/(72.27/25.4),family='Arial')+ggplot2::scale_fill_manual(values=cpal)+ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.25)))+ggplot2::labs(x=NULL,y='Younger steatosis (%)',fill=NULL,title=NULL)
 r<-clinical$models[clinical$models$model%in%c('Unadjusted','C1')&is.finite(clinical$models$OR),];r$model<-factor(r$model,levels=c('Unadjusted','C1'));r$comparison<-factor(r$comparison,levels=rev(c('G0-G1','G1-G2','G0-G2')))
 pd<-ggplot2::position_dodge(width=.45)
 p2<-ggplot2::ggplot(r,ggplot2::aes(OR,comparison,fill=model,shape=model))+ggplot2::geom_vline(xintercept=1,linetype=2,colour='#888888',linewidth=.3)+ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low,xmax=CI_high),orientation='y',width=0,position=pd,colour='black',linewidth=.35)+ggplot2::geom_point(position=pd,size=1.9,colour='black',stroke=.35)+ggplot2::scale_x_log10()+ggplot2::scale_fill_manual(values=setNames(unname(cfg$binary_colors),levels(r$model)),labels=c('Unadjusted','C1 adjusted'))+ggplot2::scale_shape_manual(values=c(21,23),labels=c('Unadjusted','C1 adjusted'))+ggplot2::labs(x='Odds ratio (log scale)',y=NULL,fill=NULL,shape=NULL,title=NULL)
 export(patch(list(p1,p2)),'Figure2_Clinical_cross_generation',height=110,multi=TRUE,panels=2)
 adult<-clinical$models[clinical$models$model%in%c('Adult_C1','C2_alcohol','C2_BMI')&is.finite(clinical$models$OR),];adult$model<-factor(adult$model,levels=c('Adult_C1','C2_alcohol','C2_BMI'),labels=c('Adult C1','C2: alcohol','C2: alcohol + BMI'))
 adult$comparison<-factor(adult$comparison,levels=c('G0-G1','G1-G2','G0-G2'))
 adult$model<-factor(as.character(adult$model),levels=rev(c('Adult C1','C2: alcohol','C2: alcohol + BMI')))
 p<-ggplot2::ggplot(adult,ggplot2::aes(OR,model))+ggplot2::geom_vline(xintercept=1,linetype=2,colour='#888888',linewidth=.3)+ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low,xmax=CI_high),orientation='y',width=0,linewidth=.35)+ggplot2::geom_point(size=1.6,colour='#555555')+ggplot2::facet_wrap(~comparison,nrow=1)+ggplot2::scale_x_log10()+ggplot2::labs(x='Odds ratio (log scale)',y=NULL)+v3_style()
 export(p,'Supplement1_Clinical_adult_C2',height=95)
 if(is.null(similarity))return(ov_bind(records))
 s<-similarity$distribution[similarity$distribution$model=='P1',];s$type<-factor(s$type,levels=c('BF','IF'));s$comparison<-factor(s$comparison,levels=rev(c('G0-G1','G1-G2','G0-G2')))
 connectors<-merge(s[s$type=='BF',c('comparison','mean')],s[s$type=='IF',c('comparison','mean')],by='comparison',suffixes=c('_BF','_IF'))
 p1<-ggplot2::ggplot(s,ggplot2::aes(mean,comparison,fill=type,shape=type))+ggplot2::geom_segment(data=connectors,ggplot2::aes(x=mean_IF,xend=mean_BF,y=comparison,yend=comparison),inherit.aes=FALSE,colour='#B5B5B5',linewidth=.85)+ggplot2::geom_point(size=2.1,colour='black',stroke=.3)+ggplot2::geom_text(ggplot2::aes(label=sprintf('%.4f',mean)),nudge_y=-.23,size=5.7/(72.27/25.4),family='Arial',show.legend=FALSE)+ggplot2::scale_fill_manual(values=familypal)+ggplot2::scale_shape_manual(values=c(BF=21,IF=24))+ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(.10,.12)))+ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=.55))+ggplot2::labs(x='Mean Bray-Curtis distance',y=NULL,fill=NULL,shape=NULL,title=NULL)
 s<-similarity$effects[similarity$effects$model=='P2',];s$comparison<-factor(s$comparison,levels=rev(c('G0-G1','G1-G2','G0-G2')))
 p2<-ggplot2::ggplot(s,ggplot2::aes(comparison,IF_minus_BF))+ggplot2::geom_hline(yintercept=0,linetype=2,colour='#888888',linewidth=.3)+ggplot2::geom_col(fill=unname(cfg$binary_colors['Comparison']),width=.48,colour='black',linewidth=.25)+ggplot2::geom_text(ggplot2::aes(label=sprintf('%.4f',IF_minus_BF)),hjust=1.15,size=5.7/(72.27/25.4),family='Arial')+ggplot2::scale_y_continuous(limits=c(min(s$IF_minus_BF)-.015,.004),expand=ggplot2::expansion(mult=0))+ggplot2::scale_x_discrete(expand=ggplot2::expansion(add=.55))+ggplot2::coord_flip()+ggplot2::labs(x=NULL,y='Adjusted IF - BF distance',title=NULL)
 export(patch(list(p1,p2)),'Figure3_Family_microbiome_similarity',height=85,multi=TRUE,panels=2)
 m<-x$metadata[v3_micro_samples(x)$B1$mask,];m$steatosis<-factor(m$fatty_liver,levels=0:1,labels=c('No steatosis','Steatosis'))
 p1<-ggplot2::ggplot(m,ggplot2::aes(generation,shannon,fill=steatosis))+ggplot2::geom_boxplot(position=ggplot2::position_dodge(.75),width=.65,outlier.shape=NA,colour='black',linewidth=.25)+ggplot2::geom_point(position=ggplot2::position_jitterdodge(jitter.width=.18,dodge.width=.75,seed=20261024),colour='#888888',size=.22,alpha=.22,show.legend=FALSE)+ggplot2::scale_fill_manual(values=binary)+ggplot2::labs(x=NULL,y='Shannon index',fill=NULL,title=NULL)
 a<-alpha$effects[alpha$effects$model=='B1'&alpha$effects$generation!='Overall',];a$generation<-factor(a$generation,levels=rev(c('G0','G1','G2')))
 p2<-ggplot2::ggplot(a,ggplot2::aes(beta,generation,colour=generation))+ggplot2::geom_vline(xintercept=0,linetype=2,colour='#888888',linewidth=.3)+ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low,xmax=CI_high),orientation='y',width=0,linewidth=.4)+ggplot2::geom_point(size=1.9,shape=18)+ggplot2::scale_colour_manual(values=cfg$generation_colors,breaks=c('G0','G1','G2'))+ggplot2::guides(colour=ggplot2::guide_legend(nrow=1))+ggplot2::labs(x='Adjusted Shannon difference',y=NULL,colour=NULL,title=NULL)
 export(patch(list(p1,p2)),'Figure4_Shannon_diversity',height=110,multi=TRUE,panels=2)
 a<-alpha$effects[alpha$effects$model%in%c('Adult_B1_same_sample','B2'),];a$model<-factor(a$model,levels=c('Adult_B1_same_sample','B2'),labels=c('Adult B1','Adult B2'))
 a$generation<-factor(a$generation,levels=rev(c('Overall','G0','G1','G2')))
 p<-ggplot2::ggplot(a,ggplot2::aes(beta,generation,colour=model,shape=model))+ggplot2::geom_vline(xintercept=0,linetype=2,colour='#888888',linewidth=.3)+ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low,xmax=CI_high),orientation='y',width=0,position=pd,colour='black',linewidth=.3)+ggplot2::geom_point(position=pd,size=1.8)+ggplot2::scale_colour_manual(values=setNames(unname(cfg$binary_colors),levels(a$model)))+ggplot2::labs(x='Adjusted Shannon difference',y=NULL,colour=NULL,shape=NULL)+v3_style()
 export(p,'Supplement2_Shannon_adult_B2',height=110)
 if(is.null(taxa))return(ov_bind(records))
 if(!is.null(composition)){
 c<-composition$effects[composition$effects$scheme=='family',];c$model<-factor(c$model,levels=c('B1','Adult_B1_same_sample','B2'),labels=c('All-age B1','Adult B1','Adult B2'));c$R2_percent<-100*c$R2_total
 p<-ggplot2::ggplot(c,ggplot2::aes(model,R2_percent))+ggplot2::geom_col(width=.55,fill='#9A9A9A',colour='black',linewidth=.3)+ggplot2::geom_text(ggplot2::aes(label=sprintf('%.3f%%',R2_percent)),vjust=-.6,size=6/(72.27/25.4),family='Arial')+ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.18)))+ggplot2::labs(x=NULL,y='Total composition variance explained (%)')+v3_style()
 export(p,'Figure5_Community_composition',height=95)
 }
 volcanoes<-list()
 for(nm in c('B1','B2')){
  t<-taxa$tables[[nm]];z<-t[t$analysis_valid,];z$y<--log10(pmax(z$P_raw,.Machine$double.xmin));z$direction<-factor(ifelse(z$P_raw>=.05,'P >= 0.05',ifelse(z$beta_log2<0,'Lower in steatosis','Higher in steatosis')),levels=c('Lower in steatosis','P >= 0.05','Higher in steatosis'))
  labels<-ov_bind(lapply(c('Lower in steatosis','Higher in steatosis'),function(k){j<-which(z$direction==k);z[head(j[order(z$P_raw[j])],4),,drop=FALSE]}))
  palette<-c('Lower in steatosis'=unname(cfg$binary_colors['Reference']),'P >= 0.05'='#B8B8B8','Higher in steatosis'=unname(cfg$binary_colors['Comparison']))
  extent<-max(abs(z$beta_log2))*2.1
  labels<-ov_bind(lapply(split(labels,labels$direction,drop=TRUE),function(a){a<-a[order(a$y,decreasing=TRUE),];side<-ifelse(a$beta_log2>0,1,-1);a$label_y<-seq(1,.55,length.out=nrow(a))*max(z$y);a$label_x<-side*.98*extent;a$line_end_x<-side*.58*extent;a$text_hjust<-ifelse(side>0,1,0);a}))
  p<-ggplot2::ggplot(z,ggplot2::aes(beta_log2,y))+ggplot2::geom_hline(yintercept=-log10(.05),linetype=2,colour='#888888',linewidth=.3)+ggplot2::geom_vline(xintercept=0,linetype=2,colour='#999999',linewidth=.3)+ggplot2::geom_point(ggplot2::aes(colour=direction),size=1.2,alpha=.85)+ggplot2::geom_segment(data=labels,ggplot2::aes(xend=line_end_x,yend=label_y),linewidth=.2,colour='#555555')+ggplot2::geom_text(data=labels,ggplot2::aes(x=label_x,y=label_y,label=SGB_id,hjust=text_hjust),family='Arial',size=5.7/(72.27/25.4))+ggplot2::scale_colour_manual(values=palette)+ggplot2::scale_x_continuous(limits=c(-extent,extent),expand=ggplot2::expansion(mult=0))+ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.02,.15)))+ggplot2::labs(x='Adjusted log2 fold change (steatosis / no steatosis)',y='-log10(raw P)',colour=NULL,title=paste0(nm,': ',unique(t$N),' participants; ',unique(t$cases),' cases; ',nrow(z),' valid SGB models'))+v3_style()
  volcanoes[[nm]]<-p+ggplot2::labs(title=paste0(if(nm=='B1')'All-age B1'else'Adult B2',' (n = ',unique(t$N),')\n',unique(t$cases),' cases; ',nrow(z),' SGBs'),x='Adjusted log2 fold change\n(steatosis / no steatosis)')
  v3_table(z,paste0('Figure6_',nm,'_volcano_source_EN'),out)
  d<-taxa$designs[[nm]];no<-d$metadata$fatty_liver==0;yes<-d$metadata$fatty_liver==1
  detected_no<-colSums(d$relative[no,,drop=FALSE]>0);detected_yes<-colSums(d$relative[yes,,drop=FALSE]>0)
  supported<-t$analysis_valid&detected_no[t$feature_id]>=10&detected_yes[t$feature_id]>=10
  pick<-head(order(ifelse(supported,t$P_raw,Inf)),8);pick<-pick[is.finite(t$P_raw[pick])&supported[pick]]
  long<-ov_bind(lapply(t$feature_id[pick],function(id)data.frame(SGB=sub('^t__','',id),group=factor(d$metadata$fatty_liver,levels=0:1,labels=c('No steatosis','Steatosis')),relative_percent=100*d$relative[,id])))
  long$SGB<-factor(long$SGB,levels=t$SGB_id[pick]);p<-ggplot2::ggplot(long,ggplot2::aes(group,relative_percent,fill=group))+ggplot2::geom_boxplot(width=.6,colour='black',linewidth=.25,outlier.shape=NA)+ggplot2::geom_point(position=ggplot2::position_jitter(width=.12,height=0,seed=20261024),colour='#888888',size=.20,alpha=.22,show.legend=FALSE)+ggplot2::facet_wrap(~SGB,ncol=2,scales='free_y')+ggplot2::scale_fill_manual(values=binary)+ggplot2::labs(x=NULL,y='Relative abundance (%)',fill=NULL)+v3_style()+ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),panel.spacing=grid::unit(5,'mm'))
  export(p,paste0('Figure7_',nm,'_Representative_abundance'),height=200);v3_table(t[pick,],paste0('Figure7_',nm,'_display_species_EN'),out)
 }
 combined<-patch(volcanoes)+patchwork::plot_layout(guides='collect')+patchwork::plot_annotation(theme=ggplot2::theme(legend.position='bottom'))
 export(combined,'Figure6_Overall_and_adult_SGB_volcano',height=110,multi=TRUE,panels=2)
 record<-ov_bind(records)
 manifest<-file.path(out,'tables','Figure_export_manifest_EN.csv')
 if(!is.null(only)&&file.exists(manifest)){
  previous<-read.csv(manifest,check.names=FALSE);previous<-previous[!previous$figure%in%c(record$figure,'Figure6_B1_SGB_volcano','Figure6_B2_SGB_volcano'),,drop=FALSE];record<-ov_bind(list(previous,record))
 }
 v3_table(record,'Figure_export_manifest_EN',out);record
}
