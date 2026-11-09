# Inert extensions. Classification and source descriptions never change fitted
# coefficients, pathway selection, input normalization or original BH families.
pfs_primary <- function(r) r$effects[r$effects$model=='B1'&r$effects$generation=='Overall',,drop=FALSE]
pfs_classify <- function(r,cfg,out) {
  path<-file.path(cfg$project_root,'data/reference/pathway_redesign_20261102/map_metacyc-pwy_lineage.tsv')
  map<-read.delim(path,header=FALSE,sep='\t',quote='',comment.char='',stringsAsFactors=FALSE,
    col.names=c('feature_id','lineage'))
  stopifnot(ncol(map)==2,all(nzchar(map$feature_id)),all(nzchar(map$lineage)))
  map<-unique(map);tokens<-strsplit(map$lineage,'|',fixed=TRUE)
  map$level_1<-vapply(tokens,`[`,character(1),1)
  map$level_2<-vapply(tokens,function(v)if(length(v)>=2L)v[2]else v[1],character(1))
  map<-map[order(map$feature_id,map$lineage),,drop=FALSE]
  ann<-r$annotation
  rows<-lapply(ann$feature_id,function(id){
    z<-map[map$feature_id==id,,drop=FALSE]
    data.frame(feature_id=id,mapped=nrow(z)>0,lineages=nrow(z),
      all_lineages=paste(z$lineage,collapse='; '),
      all_level_1=paste(sort(unique(z$level_1)),collapse='; '),
      all_level_2=paste(sort(unique(z$level_2)),collapse='; '),
      display_parent=if(nrow(z))sort(unique(z$level_2))[1]else'Unmapped',
      multiple_parents=nrow(z)>1,stringsAsFactors=FALSE)
  })
  annotation<-merge(ann,ov_bind(rows),by='feature_id',sort=FALSE,all.x=TRUE)
  long<-merge(ann,map,by='feature_id',sort=FALSE,all.x=TRUE)
  long$mapping_status<-ifelse(is.na(long$lineage),'Unmapped','Historical MetaCyc lineage')
  primary<-pfs_primary(r);sig<-primary[tf_bool(primary$FDR_lt_005),,drop=FALSE]
  sig<-sig[order(sig$P_BH,sig$P_raw,sig$feature_id),,drop=FALSE]
  sig$display_parent<-annotation$display_parent[match(sig$feature_id,annotation$feature_id)]
  selected<-head(sig[!duplicated(sig$display_parent),,drop=FALSE],10L)
  selected$display_rank<-seq_len(nrow(selected));stopifnot(nrow(selected)>0L)
  selected<-merge(selected,annotation[,c('feature_id','all_lineages','all_level_2','mapped','multiple_parents')],
    by='feature_id',sort=FALSE,all.x=TRUE)
  selected<-selected[order(selected$display_rank),,drop=FALSE]
  selected$selection_rule<-'First BH discovery per lexical level-2 display parent, ordered by q, P and ID; maximum ten parents'
  e<-merge(long[,c('feature_id','level_1','level_2','lineage','mapping_status')],primary,
    by='feature_id',sort=FALSE,all.y=TRUE)
  e$level_1[is.na(e$level_1)]<-'Unmapped';e$level_2[is.na(e$level_2)]<-'Unmapped'
  catrows<-lapply(split(e,interaction(e$level_1,e$level_2,drop=TRUE)),function(z){
    z<-z[!duplicated(z$feature_id),,drop=FALSE]
    data.frame(level_1=z$level_1[1],level_2=z$level_2[1],tested=nrow(z),
      discoveries=sum(tf_bool(z$FDR_lt_005)),lower=sum(tf_bool(z$FDR_lt_005)&z$beta_log2<0),
      higher=sum(tf_bool(z$FDR_lt_005)&z$beta_log2>0))
  })
  provenance<-data.frame(resource='Eric Franzosa public MetaCyc lineage map',
    URL='https://forum.biobakery.org/uploads/short-url/hDx8TKiUH5eeweMXw4TY94zVYUQ.tsv',
    context_URL='https://forum.biobakery.org/t/metacyc-hierarchy-to-invetigate-identify-specific-pathways/1830',
    accessed='2026-11-02',version='Historical file shared in 2021; not current MetaCyc release or mandatory HUMAnN hierarchy',
    md5=unname(tools::md5sum(path)),
    SHA256='B69F027CA898EA1F7906B3DBFB47923E813348208BECE17AE32AB112BF0DC32D',
    display_rule='Lexically first level-2 parent for layout only; full multi-parent map retained')
  ade_table(annotation,'MetaCyc_annotations_all_pathways_EN',out)
  ade_table(long,'MetaCyc_all_parent_lineages_EN',out)
  ade_table(ov_bind(catrows),'Functional_category_descriptive_counts_EN',out)
  ade_table(provenance,'MetaCyc_mapping_provenance_EN',out)
  ade_table(selected,'Main_figure_representatives_EN',out)
  r$protected<-ov_bind(list(r$protected,ade_fingerprint(path)))
  r$classification<-list(annotation=annotation,lineages=long,counts=ov_bind(catrows),provenance=provenance)
  r$representatives<-selected;r
}

pfs_sources <- function(r,cfg,out) {
  x<-readRDS(file.path(out,'private/functional_input.rds'))
  z<-readRDS(x$source_files[1]);rn<-rownames(z);cn<-colnames(z)
  rd<-as.data.frame(SummarizedExperiment::rowData(z))
  source_rows<-which(grepl('|',rn,fixed=TRUE))
  labels<-rn[source_rows];taxon<-sub('^[^|]*\\|','',labels)
  id<-sub(':.*','',sub('\\|.*','',labels))
  stopifnot(all(id%in%x$features$feature_id),!anyDuplicated(labels))
  genus<-sub('^g__','',sub('\\..*$','',taxon))
  genus[!grepl('^g__',taxon)]<-'Unclassified'
  genus[grepl('unclassified|unidentified|unknown|^$',genus,ignore.case=TRUE)]<-'Unclassified'
  species<-if('Species'%in%names(rd))as.character(rd$Species[source_rows])else sub('.*s__','',taxon)
  species[is.na(species)]<-''
  source_annotation<-data.frame(source_row=source_rows,feature_id=id,source_label=labels,
    taxon_label=taxon,genus=gsub('_',' ',genus,fixed=TRUE),species=species,stringsAsFactors=FALSE)
  ids<-r$fits$B1$sample_ids;sample_indices<-match(ids,cn);mi<-match(ids,x$metadata$sample_id)
  stopifnot(!anyNA(sample_indices),!anyNA(mi),identical(ids,as.character(x$metadata$sample_id[mi])))
  m<-x$metadata[mi,,drop=FALSE];denom<-x$denominator[mi]
  cd<-as.data.frame(SummarizedExperiment::colData(z))
  v<-intersect(c('tutknom18','TUTKNOM18'),names(cd))
  if(length(v))stopifnot(all(as.character(cd[[v[1]]][sample_indices])==as.character(m$participant_id)))
  stopifnot(all(denom>0),nrow(m)==length(ids),sum(m$fatty_liver)==unique(r$fits$B1$effects$cases))
  desc<-list(dataset='assay001',n_samples=length(cn),n_source_rows=length(rn),
    zero_based_row_indices=source_rows-1L,zero_based_sample_indices=sample_indices-1L)
  manifest<-file.path(out,'private/source_import_indices.json')
  jsonlite::write_json(desc,manifest,auto_unbox=TRUE,pretty=TRUE)
  chunkdir<-file.path(out,'private/source_chunks')
  ade_log('Reading all ',length(source_rows),' actual HUMAnN stratified rows in bounded chunks')
  status<-system2(cfg$python,c(shQuote(file.path(cfg$project_root,'tools/read_humann_h5_chunks.py')),
    shQuote(x$source_files[2]),shQuote(manifest),shQuote(chunkdir)))
  if(status!=0L)stop('Read-only source import failed')
  chunks<-jsonlite::fromJSON(file.path(chunkdir,'chunks.json'))$chunks
  groups<-list(All=seq_len(nrow(m)),No_steatosis=which(m$fatty_liver==0),Steatosis=which(m$fatty_liver==1))
  sums<-matrix(0,nrow(x$relative),nrow(m),dimnames=list(rownames(x$relative),ids))
  summaries<-checks<-list()
  for(j in seq_len(nrow(chunks))){
    cc<-chunks[j,];ix<-seq.int(cc$start_one_based,length.out=cc$rows)
    con<-file(file.path(chunkdir,cc$file),'rb')
    a<-tryCatch(matrix(readBin(con,'double',n=cc$rows*nrow(m),size=8,endian='little'),nrow=cc$rows),finally=close(con))
    stopifnot(ncol(a)==nrow(m),all(is.finite(a)),all(a>=0))
    a<-sweep(a,2,denom,'/');rownames(a)<-as.character(ix)
    sa<-rowsum(a,id[ix],reorder=FALSE)
    ii<-match(rownames(sa),rownames(sums));sums[ii,]<-sums[ii,,drop=FALSE]+sa
    for(g in names(groups)){
      k<-groups[[g]]
      summaries[[length(summaries)+1L]]<-data.frame(source_index=ix,group=g,N=length(k),
        mean_RA_percent=rowMeans(a[,k,drop=FALSE])*100,
        detection_percent=rowMeans(a[,k,drop=FALSE]>0)*100)
    }
    # An independently written row-wise reduction checks the chunk summary.
    checkix<-unique(c(1L,nrow(a)));ii<-groups$Steatosis
    ref<-vapply(checkix,function(k)sum(a[k,ii])/length(ii)*100,numeric(1))
    measured<-tail(summaries,1)[[1]]$mean_RA_percent[checkix]
    checks[[j]]<-data.frame(chunk=j,max_mean_error=max(abs(ref-measured)),finite_nonnegative=TRUE)
    rm(a,sa);gc(FALSE)
  }
  s<-ov_bind(summaries);k<-s$source_index
  for(v in names(source_annotation))s[[v]]<-source_annotation[[v]][k]
  s<-s[,c(names(source_annotation),'group','N','mean_RA_percent','detection_percent')]
  stopifnot(!anyDuplicated(s[,c('source_row','group')]))
  # Sum the per-species means: this is exactly the mean of the summed observed
  # strata, not a reconstruction of HUMAnN's community pathway estimate.
  genus<-stats::aggregate(mean_RA_percent~feature_id+genus+group,s,FUN=sum)
  totals<-stats::aggregate(mean_RA_percent~feature_id+group,s,FUN=sum)
  names(totals)[3]<-'summed_stratified_mean_RA_percent'
  genus<-merge(genus,totals,by=c('feature_id','group'),sort=FALSE,all.x=TRUE)
  genus$source_share_percent<-ifelse(genus$summed_stratified_mean_RA_percent>0,
    100*genus$mean_RA_percent/genus$summed_stratified_mean_RA_percent,NA_real_)
  s<-merge(s,totals,by=c('feature_id','group'),sort=FALSE,all.x=TRUE)
  s$source_share_percent<-ifelse(s$summed_stratified_mean_RA_percent>0,
    100*s$mean_RA_percent/s$summed_stratified_mean_RA_percent,NA_real_)
  ids_path<-unique(id);community<-x$relative[,match(ids,colnames(x$relative)),drop=FALSE]
  audit<-ov_bind(lapply(ids_path,function(path){
    a<-sums[path,];b<-community[path,];den<-pmax(b,.Machine$double.eps)
    data.frame(feature_id=path,N=length(a),stratified_carriers=sum(a>0),community_carriers=sum(b>0),
      source_rows=sum(id==path),mean_sum_strata_percent=mean(a)*100,mean_community_RA_percent=mean(b)*100,
      mean_strata_to_community=if(mean(b)>0)mean(a)/mean(b)else NA_real_,
      max_absolute_nonadditivity=max(abs(a-b)),
      maximum_ratio_in_community_carriers=if(any(b>0))max((a/den)[b>0])else NA_real_)
  }))
  # Check that genus aggregation retained all supplied strata, including zeros.
  expected<-ov_bind(lapply(names(groups),function(g)data.frame(feature_id=rownames(sums),group=g,
    expected_mean_RA_percent=rowMeans(sums[,groups[[g]],drop=FALSE])*100)))
  tt<-merge(totals,expected,by=c('feature_id','group'))
  stopifnot(max(abs(tt$summed_stratified_mean_RA_percent-tt$expected_mean_RA_percent))<1e-9)
  reps<-r$representatives$feature_id
  pool<-genus[genus$group=='All'&genus$feature_id%in%reps&genus$genus!='Unclassified',,drop=FALSE]
  scores<-stats::aggregate(source_share_percent~genus,pool,FUN=sum)
  scores<-scores[order(-scores$source_share_percent,scores$genus),,drop=FALSE]
  top<-head(scores$genus,8L)
  genus$display_genus<-ifelse(genus$genus=='Unclassified','Unclassified',ifelse(genus$genus%in%top,genus$genus,'Other named genera'))
  display<-stats::aggregate(mean_RA_percent~feature_id+group+display_genus,genus,FUN=sum)
  display<-merge(display,totals,by=c('feature_id','group'),all.x=TRUE,sort=FALSE)
  display$source_share_percent<-ifelse(display$summed_stratified_mean_RA_percent>0,
    100*display$mean_RA_percent/display$summed_stratified_mean_RA_percent,NA_real_)
  display<-display[display$feature_id%in%reps,,drop=FALSE]
  add<-expand.grid(feature_id=reps,group=names(groups),display_genus=c(top,'Other named genera','Unclassified'),stringsAsFactors=FALSE)
  display<-merge(add,display,by=c('feature_id','group','display_genus'),all.x=TRUE,sort=FALSE)
  # A missing taxon is a true zero only when that pathway/group has observed
  # strata. An entirely unstratified pathway/group has an undefined share.
  at<-match(paste(display$feature_id,display$group),paste(totals$feature_id,totals$group))
  missing<-is.na(display$source_share_percent)
  has_total<-!is.na(at)&totals$summed_stratified_mean_RA_percent[at]>0
  display$source_share_percent[missing&has_total]<-0
  display$descriptive_only<-TRUE
  ade_table(source_annotation,'HUMAnN_all_source_annotations_EN',out)
  ade_table(s,'HUMAnN_all_species_pathway_sources_EN',out)
  ade_table(genus,'HUMAnN_all_genus_pathway_sources_EN',out)
  ade_table(audit,'Community_vs_summed_strata_audit_EN',out)
  ade_table(display,'Main_figure_descriptive_source_shares_EN',out)
  ade_table(data.frame(rank=seq_along(top),genus=top,rule='Pooled source shares across main representative pathways; outcome-blind genus ranking'),
    'Displayed_source_genera_EN',out)
  fl_write_csv(ov_bind(checks),file.path(out,'audit/source_chunk_independent_mean_checks.csv'))
  r$sources<-list(genus=genus,totals=totals,audit=audit,display=display,top_genera=top,
    source_rows=length(source_rows),taxon_strata=length(unique(taxon)),samples=nrow(m),
    groups=data.frame(group=names(groups),N=lengths(groups)),
    max_aggregation_error=max(abs(tt$summed_stratified_mean_RA_percent-tt$expected_mean_RA_percent)))
  ade_log('Source summaries complete; all ',length(source_rows),' rows retained; no community-additivity assumption')
  r
}

pfs_taxon_context <- function(r,cfg,out) {
  src<-tf_source(cfg$project_root)
  baseline<-readRDS(src$files[1]);ann<-ade_labels(baseline$features)
  primary<-tf_baseline(src)$B1
  named<-ann[!is.na(ann$scientific_name),,drop=FALSE]
  named$key<-tolower(named$scientific_name)
  sources<-tf_read(file.path(out,'tables/HUMAnN_all_source_annotations_EN.csv'))
  tax<-unique(sources[,c('taxon_label','genus','species')])
  tax$key<-tolower(gsub('_',' ',sub('^s__','',tax$species),fixed=TRUE))
  rows<-lapply(seq_len(nrow(tax)),function(i){
    q<-tax[i,];hits<-named[named$key==q$key,,drop=FALSE]
    valid_name<-grepl('^(candidatus )?[a-z][a-z]+ [a-z][a-z-]+$',q$key)&&
      !grepl('\\b(bacterium|archaeon|unclassified|uncultured|unidentified|unknown|sp|spp)\\b',q$key,perl=TRUE)
    unique_hit<-valid_name&&nrow(hits)==1L
    id<-if(unique_hit)hits$feature_id else NA_character_
    k<-match(id,primary$feature_id)
    data.frame(taxon_label=q$taxon_label,humann_genus=q$genus,humann_species=q$species,
      exact_named_SGB_candidates=paste(hits$feature_id,collapse='; '),
      mapping_status=if(!valid_name)'HUMAnN name unresolved'else if(nrow(hits)==0L)'No exact named SGB'else if(nrow(hits)>1L)'Multiple SGBs; no forced match'else'Unique exact binomial match',
      feature_id=id,SGB_beta_log2=if(unique_hit)primary$beta_log2[k]else NA_real_,
      SGB_P_BH=if(unique_hit)primary$P_BH[k]else NA_real_,
      SGB_primary_discovery=unique_hit&&!is.na(k)&&tf_bool(primary$FDR_lt_005[k]),stringsAsFactors=FALSE)
  })
  crosswalk<-ov_bind(rows)
  links<-merge(sources,crosswalk,by='taxon_label',all.x=TRUE,sort=FALSE,
    suffixes=c('_pathway','_SGB'))
  ade_table(crosswalk,'HUMAnN_species_to_existing_SGB_exact_crosswalk_EN',out)
  ade_table(links,'HUMAnN_pathway_source_existing_SGB_context_EN',out)
  counts<-as.data.frame(table(mapping_status=crosswalk$mapping_status),stringsAsFactors=FALSE)
  names(counts)[2]<-'source_taxon_labels';ade_table(counts,'HUMAnN_to_SGB_mapping_counts_EN',out)
  r$source_taxonomy_context<-list(crosswalk=crosswalk,counts=counts,
    unique_matches=sum(crosswalk$mapping_status=='Unique exact binomial match'),
    matched_primary_SGBs=length(unique(crosswalk$feature_id[crosswalk$SGB_primary_discovery])))
  rm(baseline);gc(FALSE);r
}

pfs_composition <- function(r,cfg,out) {
  stopifnot(as.character(packageVersion('vegan'))=='2.7.6')
  x<-readRDS(file.path(out,'private/functional_input.rds'))
  results<-audit<-list();scopes<-v3_micro_samples(x)
  for(nm in names(scopes)){
    sc<-scopes[[nm]];d<-v3_design(x,sc$mask,sc$fields)
    m<-d$metadata;md<-d$design;md$fatty_liver<-m$fatty_liver
    retained<-r$fits[[nm]]$filter$feature_id[r$fits[[nm]]$filter$retained]
    stopifnot(identical(as.character(m$sample_id),r$fits[[nm]]$sample_ids))
    A<-x$relative[retained,match(m$sample_id,colnames(x$relative)),drop=FALSE]
    # Do not re-close the retained pathways: the supplied 620-row denominator
    # remains identical to that used by the feature models and source tables.
    Z<-model.matrix(reformulate(c('family_id',d$nuisance)),md)
    q<-qr(Z,tol=1e-9);keep<-q$pivot[seq_len(q$rank)];Z<-Z[,sort(keep),drop=FALSE]
    kept_nuis<-d$nuisance[d$nuisance%in%colnames(Z)];condition<-c('family_id',kept_nuis)
    Z<-model.matrix(reformulate(condition),md);stopifnot(qr(Z)$rank==ncol(Z))
    signature<-digest::digest(list(protocol=pfs_protocol,relative=A,design=Z,outcome=m$fatty_liver,
      seed=cfg$v3_seed,permutations=cfg$v3_permutations,vegan='2.7.6'),algo='sha256')
    fitpath<-file.path(out,'private',paste0('composition_fit_checkpoint_',nm,'.rds'))
    if(file.exists(fitpath)){
      old<-readRDS(fitpath);stopifnot(identical(old$signature,signature));fit<-old$fit;setup<-old$setup
    }else{
      ade_log('Functional composition ',nm,': ',nrow(m),' people, ',nrow(A),' eligible pathways')
      D<-vegan::vegdist(t(A),method='bray')
      form<-as.formula(paste('D ~ fatty_liver + Condition(',paste(condition,collapse=' + '),')'))
      fit<-vegan::dbrda(form,data=md,add='lingoes');stopifnot(fit$CCA$qrank==1L)
      setup<-fl_dv_fast_setup(fit,m,Z)
      saveRDS(list(signature=signature,fit=fit,setup=setup),fitpath,compress='gzip');rm(D);gc(FALSE)
    }
    residual_case<-m$fatty_liver-ave(m$fatty_liver,m$family_id);Qc<-setup$Q[,-1,drop=FALSE]
    compact<-fit;compact$Ybar<-setup$E;compact$pCCA$QR<-qr(Qc)
    compact$pCCA$rank<-ncol(Qc);compact$CCA$QR<-qr(cbind(Qc,residual_case))
    for(scheme in c('family','family_generation')){
      key<-paste(nm,scheme);cache<-file.path(out,'private',paste0('composition_result_',nm,'_',scheme,'.rds'))
      if(file.exists(cache)){
        old<-readRDS(cache);stopifnot(identical(old$signature,signature))
        results[[key]]<-old$result;audit[[key]]<-old$audit;next
      }
      block<-if(scheme=='family')as.character(m$family_id)else as.character(interaction(m$family_id,m$generation,drop=TRUE))
      P<-fl_dv_permutations(block,cfg$v3_permutations,
        cfg$v3_seed+match(nm,names(scopes))*10L+match(scheme,c('family','family_generation')))
      prefix<-file.path(out,'private',paste0('pathway_vegan276_',nm,'_',scheme))
      wb<-function(v,s){con<-file(paste0(prefix,s),'wb');tryCatch(writeBin(v,con,endian='little'),finally=close(con))}
      wb(as.double(setup$E),'_E.bin');wb(as.double(Qc),'_Q.bin')
      wb(as.double(residual_case),'_X.bin');wb(as.integer(P),'_P.bin')
      ade_log('Functional composition ',nm,' ',scheme,': ',nrow(P),' restricted draws')
      status<-system2(cfg$python,c(shQuote(file.path(cfg$project_root,'scripts/partial_dbrda_vegan276_kernel.py')),
        '--prefix',shQuote(prefix),'--n',nrow(m),'--k',ncol(Qc),'--permutations',nrow(P),'--df',setup$residual_df),
        stdout=paste0(prefix,'.log'),stderr=paste0(prefix,'.log'))
      if(status!=0L)stop('Audited vegan kernel failed: ',prefix)
      null<-read.csv(paste0(prefix,'_F.csv'))$F_null
      stopifnot(length(null)==nrow(P),all(is.finite(null)))
      ade_log('Independent compact and original vegan checks: ',nm,' ',scheme)
      test<-vegan::permutest(compact,permutations=P[1:10,,drop=FALSE],model='reduced',parallel=1L)
      ratio<-setup$residual_df/test$df[2];reference<-as.numeric(test$F.perm)*ratio
      err<-max(abs(head(null,10)-reference))
      fl_write_csv(data.frame(draw=1:10,kernel=head(null,10),compact_vegan=reference,
        absolute_error=abs(head(null,10)-reference)),
        file.path(out,'audit',paste0('composition_equivalence_',nm,'_',scheme,'.csv')))
      if(err>1e-7)stop('Current-version compact check failed: ',err)
      full<-vegan::permutest(fit,permutations=P[1,,drop=FALSE],model='reduced',parallel=1L)
      fullerr<-abs(as.numeric(full$F.perm)[1]-reference[1])
      if(fullerr>1e-7)stop('Original vegan check failed: ',fullerr)
      obs<-as.numeric(full$F.0)
      p<-(1+sum(null>=obs-sqrt(.Machine$double.eps)))/(length(null)+1)
      result<-data.frame(model=nm,scheme=scheme,N=nrow(m),cases=sum(m$fatty_liver),
        families=length(unique(m$family_id)),eligible_pathways=nrow(A),pseudo_F=obs,
        R2_total=fit$CCA$tot.chi/fit$tot.chi,
        R2_partial=fit$CCA$tot.chi/(fit$CCA$tot.chi+fit$CA$tot.chi),P_raw=p,
        permutations=length(null),status='ESTIMATED',
        estimand='Conditional on primary family and model covariates; Lingoes-corrected Bray-Curtis')
      sizes<-table(block);disc<-vapply(split(m$fatty_liver,block),function(v)length(unique(v))>1,logical(1))
      diag<-data.frame(model=nm,scheme=scheme,blocks=length(sizes),movable_people=sum(sizes[sizes>1]),
        discordant_blocks=sum(disc),vegan_kernel_error=err,full_original_vegan_error=fullerr,
        vegan_version='2.7.6',permutation_protocol='Reduced response and conditions co-permuted; exposure re-residualized per draw')
      results[[key]]<-result;audit[[key]]<-diag
      saveRDS(list(signature=signature,result=result,audit=diag),cache)
      fl_write_csv(data.frame(iteration=seq_along(null),F_null=null),paste0(prefix,'_null.csv'))
      # These four exact disposable files are under this run's private folder.
      unlink(paste0(prefix,c('_E.bin','_Q.bin','_X.bin','_P.bin','_F.csv')),recursive=FALSE)
      rm(P);gc(FALSE)
    }
    rm(fit,setup,compact,A,Z);gc(FALSE)
  }
  result<-list(effects=ov_bind(results),diagnostics=ov_bind(audit))
  # Prespecified primary functional test: all ages / within-family permutation.
  # The remaining five tests are explicit sensitivity checks, with a separate
  # descriptive six-test BH column; original permutation P values are retained.
  result$effects$role<-ifelse(result$effects$model=='B1'&result$effects$scheme=='family','Primary','Sensitivity')
  result$effects$P_BH_six_tests<-p.adjust(result$effects$P_raw,'BH')
  ade_table(result$effects,'Functional_community_composition_EN',out)
  ade_table(result$diagnostics,'Functional_community_permutation_diagnostics_EN',out)
  result
}

pfs_details <- function(r,cfg,out,figures_only=FALSE) {
  needed<-c('classification','composition','sources','taxon_context')
  staged<-file.path(out,'private/source_standalone.rds')
  if(!figures_only&&file.exists(staged)&&!all(c('classification','sources')%in%r$completed_stages)){
    s<-readRDS(staged);p<-pfs_primary(r);sigcols<-c('feature_id','beta_log2','P_raw','P_BH')
    sig<-digest::digest(p[order(p$feature_id),sigcols],algo='sha256')
    stopifnot(identical(s$protocol,pfs_protocol),identical(s$primary_signature,sig))
    for(k in c('classification','representatives','sources'))r[[k]]<-s[[k]]
    map<-file.path(cfg$project_root,'data/reference/pathway_redesign_20261102/map_metacyc-pwy_lineage.tsv')
    r$protected<-unique(ov_bind(list(r$protected,ade_fingerprint(map))))
    r$completed_stages<-unique(c(r$completed_stages,'classification','sources'));pfs_save(r,out)
    ade_log('Recovered actual-source summaries; exact fresh-primary coefficient/P/q signature verified')
  }
  if(figures_only&&!all(needed%in%r$completed_stages))stop('Figures-only requires completed classification, community and source stages; use --resume')
  if(!'classification'%in%r$completed_stages){
    r<-pfs_classify(r,cfg,out);r$completed_stages<-c(r$completed_stages,'classification');pfs_save(r,out)
  }
  if(!'composition'%in%r$completed_stages){
    r$composition<-pfs_composition(r,cfg,out);r$completed_stages<-c(r$completed_stages,'composition');pfs_save(r,out)
  }
  if(!'sources'%in%r$completed_stages){
    r<-pfs_sources(r,cfg,out);r$completed_stages<-c(r$completed_stages,'sources');pfs_save(r,out)
  }
  if(!'taxon_context'%in%r$completed_stages){
    r<-pfs_taxon_context(r,cfg,out);r$completed_stages<-c(r$completed_stages,'taxon_context');pfs_save(r,out)
  }
  r
}
