# 17: prespecified SGB association models; original input and MI remain immutable.
fl_sgb_config <- function(root) {
  cfg <- fl_config(root)
  cfg$sgb_cores <- max(1L, min(4L, as.integer(Sys.getenv("FATTYLIVER_SGB_CORES", "2"))))
  cfg$sgb_seed <- 20261022L
  cfg$sgb_display_n <- 12L
  cfg
}

fl_sgb_prepare <- function(cfg) {
  d <- fl_load_processed_microbiome(cfg$project_root)
  cv <- fl_cv_load(cfg, d$metadata); m <- cv$metadata
  stopifnot(nrow(m)==1275L, length(unique(m$family_id))==319L,
    all(table(m$family_id,m$generation)>0L), nrow(d$matrices$relative)==695L)
  fields <- c("age","sex","generation","family_id","fatty_liver","prep_batch","seq_pool",
    "survey_center","log_total","abx_recent","season_sin","season_cos","freeze_delay_log","stool_gap_log")
  mask <- complete.cases(m[,fields]); sample_flow <- data.frame(
    stage=c("Eligible three-generation cohort","Core complete cases","Missing core covariates"),
    n=c(nrow(m),sum(mask),sum(!mask)))
  m <- droplevels(m[mask,,drop=FALSE]); m <- m[order(m$family_id,m$sample_id),,drop=FALSE]
  m$generation <- factor(m$generation,levels=c("G0","G1","G2")); m$family_id <- factor(m$family_id)
  for(v in c("sex","prep_batch","seq_pool","survey_center"))m[[v]] <- factor(m[[v]])
  rel <- t(d$matrices$relative[,match(m$sample_id,colnames(d$matrices$relative)),drop=FALSE])
  lg <- t(d$matrices$log10[,match(m$sample_id,colnames(d$matrices$log10)),drop=FALSE])
  stopifnot(identical(rownames(rel),m$sample_id),all(is.finite(rel)),all(rel>=0),all(rel<=1),
    ncol(rel)==695L,all(is.finite(lg)))
  age_names <- character(); knots <- list()
  for(g in levels(m$generation)){
    j <- which(m$generation==g); B <- splines::ns(m$age[j],df=3)
    for(k in 1:3){nm<-paste0("age_",g,"_",k);m[[nm]]<-0;m[j,nm]<-B[,k];age_names<-c(age_names,nm)}
    knots[[g]] <- data.frame(generation=g,knots=paste(attr(B,"knots"),collapse=";"),
      boundary=paste(attr(B,"Boundary.knots"),collapse=";"))
  }
  nuisance <- c("generation",age_names,"sex","prep_batch","seq_pool","survey_center","log_total",
    "abx_recent","season_sin","season_cos","freeze_delay_log","stool_gap_log")
  constant <- nuisance[vapply(m[nuisance],function(x)length(unique(x))<2L,logical(1))]
  Z <- model.matrix(reformulate(setdiff(nuisance,constant)),m)
  qrZ <- qr(Z); keep <- sort(qrZ$pivot[seq_len(qrZ$rank)])
  aliases <- setdiff(colnames(Z),colnames(Z)[keep]); Z<-Z[,keep,drop=FALSE]
  Z <- Z[,colnames(Z)!="(Intercept)",drop=FALSE]
  cov_names <- sprintf("cov_%03d",seq_len(ncol(Z)))
  mapping <- data.frame(model_column=cov_names,original_column=colnames(Z))
  md <- data.frame(family_id=m$family_id,fatty_liver=m$fatty_liver,row.names=m$sample_id)
  for(k in seq_along(cov_names))md[[cov_names[k]]] <- Z[,k]
  for(g in levels(m$generation))md[[paste0("case_",g)]] <- as.numeric(m$generation==g)*m$fatty_liver
  form <- function(case)as.formula(paste("~",paste(c(case,cov_names,"(1 | family_id)"),collapse=" + ")),env=globalenv())
  for(cs in list("fatty_liver",paste0("case_",levels(m$generation)))){
    X<-model.matrix(lme4::nobars(form(cs)),md)
    if(qr(X)$rank!=ncol(X))stop("SGB fixed-effect design is not full rank")
  }
  counts <- do.call(rbind,lapply(levels(m$generation),function(g){i<-m$generation==g
    data.frame(generation=g,n=sum(i),cases=sum(m$fatty_liver[i]),minors=sum(m$age[i]<18),
      families=length(unique(m$family_id[i])))}))
  support <- do.call(rbind,lapply(c("Overall",levels(m$generation)),function(g){
    j<-if(g=="Overall")rep(TRUE,nrow(m))else m$generation==g
    a<-rel[j,,drop=FALSE];y<-m$fatty_liver[j]
    data.frame(feature_id=colnames(rel),generation=g,n=nrow(a),cases=sum(y),
      detected=colSums(a>0),detected_cases=colSums(a[y==1,,drop=FALSE]>0),
      detected_controls=colSums(a[y==0,,drop=FALSE]>0),
      prevalence=colMeans(a>0),mean_relative=colMeans(a))}))
  pv<-readLines(file.path(cfg$project_root,"data/processed/CURRENT"),n=1)
  sources<-c(file.path(cfg$project_root,"data/processed",pv,"microbiome_processed.rds"),cv$source_files)
  list(metadata=m,design=md,relative=rel,log10=lg,features=d$features,
    formula_overall=form("fatty_liver"),formula_generation=form(paste0("case_",levels(m$generation))),
    mapping=mapping,aliases=aliases,constant=constant,knots=do.call(rbind,knots),
    sample_counts=counts,sample_flow=sample_flow,support=support,
    source_manifest=data.frame(path=sources,md5=unname(tools::md5sum(sources))))
}

fl_sgb_status <- function(fit, raw_error=NA_character_) {
  if(!inherits(fit,"merMod"))return("NOT_ESTIMABLE")
  b<-lme4::fixef(fit);V<-tryCatch(as.matrix(vcov(fit)),error=function(e)NULL)
  if(is.null(V)||any(!is.finite(b))||any(!is.finite(V)))return("INVALID_VARIANCE")
  conv<-fit@optinfo$conv
  msg<-conv$lme4$messages
  msg<-msg[!grepl("boundary.*singular",msg,ignore.case=TRUE)]
  if(any(unlist(conv$opt)!=0)||length(msg))return("NONCONVERGED")
  if(lme4::isSingular(fit,tol=1e-4))return("ESTIMATED_SINGULAR")
  "ESTIMATED"
}

fl_sgb_extract <- function(fit_object,which_model,scope,predictors) {
  ans<-list()
  for(component in c("abundance","prevalence")){
    fd<-fit_object[[paste0("fit_data_",component)]];rr<-fd$results
    for(feature in unique(rr$feature)){
      fit<-fd$fits[[feature]];state<-fl_sgb_status(fit)
      for(term in predictors){
        q<-rr[rr$feature==feature & rr$name==term,,drop=FALSE]
        if(nrow(q)!=1L)q<-data.frame(coef=NA_real_,stderr=NA_real_,pval_individual=NA_real_,error="Missing coefficient")
        b<-q$coef[1];se<-q$stderr[1];df<-Inf
        if(component=="abundance"&&inherits(fit,"merMod")){
          cs<-tryCatch(coef(summary(fit)),error=function(e)NULL)
          if(!is.null(cs)&&term%in%rownames(cs)&&"df"%in%colnames(cs))df<-cs[term,"df"]
        }
        if(!is.finite(df)||df<=0)df<-Inf
        critical<-if(is.finite(df))qt(.975,df)else qnorm(.975)
        ok<-startsWith(state,"ESTIMATED")&&is.finite(b)&&is.finite(se)&&se>0
        ans[[length(ans)+1L]]<-data.frame(feature_id=feature,scope=scope,component=component,
          term=term,estimate=if(ok)b else NA_real_,SE=if(ok)se else NA_real_,df=df,
          CI_low=if(ok)b-critical*se else NA_real_,CI_high=if(ok)b+critical*se else NA_real_,
          p_raw=if(ok)q$pval_individual[1]else NA_real_,OR=if(ok&&component=="prevalence")exp(b)else NA_real_,
          OR_CI_low=if(ok&&component=="prevalence")exp(b-critical*se)else NA_real_,
          OR_CI_high=if(ok&&component=="prevalence")exp(b+critical*se)else NA_real_,
          fold_change=if(ok&&component=="abundance")2^b else NA_real_,
          status=if(ok)state else if(startsWith(state,"ESTIMATED"))"INVALID_COEFFICIENT"else state,
          package_warning=if("error"%in%names(q))as.character(q$error[1])else NA_character_,
          method=which_model,stringsAsFactors=FALSE)
      }
    }
  }
  do.call(rbind,ans)
}

fl_sgb_joint <- function(a,b) {
  if(is.finite(a)&&is.finite(b))return(pbeta(min(a,b),1,2))
  if(is.finite(a))return(a)
  if(is.finite(b))return(b)
  NA_real_
}

fl_sgb_heterogeneity <- function(fit_object) {
  H<-rbind(c(1,-1,0),c(0,1,-1));terms<-paste0("case_",c("G0","G1","G2"));out<-list()
  for(component in c("abundance","prevalence")){
    fits<-fit_object[[paste0("fit_data_",component)]]$fits
    for(feature in names(fits)){
      fit<-fits[[feature]];state<-fl_sgb_status(fit);stat<-p<-df2<-NA_real_
      if(startsWith(state,"ESTIMATED")){
        test<-tryCatch({
          b<-lme4::fixef(fit);if(!all(terms%in%names(b)))stop("Missing generation effect")
          if(component=="abundance"){
            L<-matrix(0,2,length(b),dimnames=list(NULL,names(b)));L[,terms]<-H
            q<-lmerTest::contestMD(fit,L,rhs=c(0,0),ddf="Satterthwaite")
            c(stat=q[["F value"]][1],p=q[["Pr(>F)"]][1],df2=q[["DenDF"]][1])
          }else{
            V<-as.matrix(vcov(fit))[terms,terms,drop=FALSE];delta<-H%*%b[terms];W<-H%*%V%*%t(H)
            chi<-as.numeric(t(delta)%*%solve(W,delta));c(stat=chi,p=pchisq(chi,2,lower.tail=FALSE),df2=NA_real_)
          }
        },error=function(e)NULL)
        if(is.null(test))state<-"NOT_ESTIMABLE_CONTRAST"else{stat<-test["stat"];p<-test["p"];df2<-test["df2"]}
      }
      out[[length(out)+1L]]<-data.frame(feature_id=feature,component=component,statistic=stat,df1=2,df2=df2,
        p_interaction_raw=p,status=state,test=if(component=="abundance")"Satterthwaite 2-df F"else"Augmented-GLMM 2-df Wald chi-square")
    }
  }
  do.call(rbind,out)
}

# The official median routine only requests one fitted feature at a time.
# A disk-backed S3 list keeps its unmodified algorithm and global reference,
# without retaining hundreds of large merMod objects in memory.
fl_sgb_disk_fits <- function(index) {
  stopifnot(!is.null(names(index)),!anyDuplicated(names(index)))
  structure(as.list(index),class="fl_sgb_disk_fits",cache=new.env(parent=emptyenv()))
}

`[[.fl_sgb_disk_fits` <- function(x,i,...) {
  path<-unclass(x)[[i]]
  if(is.null(path))return(NULL)
  cache<-attr(x,"cache")
  if(!identical(cache$path,path)){
    if(length(ls(cache,all.names=TRUE)))rm(list=ls(cache,all.names=TRUE),envir=cache)
    gc(FALSE)
    cache$fits<-readRDS(path)$fit$fit_data_abundance$fits
    cache$path<-path
  }
  feature<-if(is.character(i))i else names(x)[i]
  cache$fits[[feature]]
}

fl_sgb_fit <- function(x,cfg,out) {
  cores<-cfg$sgb_cores
  old_pb<-pbapply::pboptions(type="txt",use_lb=TRUE)
  on.exit(pbapply::pboptions(old_pb),add=TRUE)
  call_fit<-function(label,formula){
    dest<-file.path(out,"private",label);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    batches<-split(seq_len(ncol(x$relative)),ceiling(seq_len(ncol(x$relative))/20L))
    pieces<-list();ledger<-list()
    for(i in seq_along(batches)){
      ids<-batches[[i]];bd<-file.path(dest,sprintf("batch_%03d",i));dir.create(bd,showWarnings=FALSE)
      checkpoint<-file.path(bd,"complete_fit.rds")
      signature<-digest::digest(list(x$relative[,ids,drop=FALSE],x$design,deparse(formula),
        packageVersion("maaslin3")),algo="sha256")
      cat(format(Sys.time())," MaAsLin 3 ",label," batch ",i,"/",length(batches),
        "; features ",min(ids),"-",max(ids),"; cores=",cores,"\n",sep="");flush.console()
      started<-Sys.time()
      if(file.exists(checkpoint)){
        cached<-readRDS(checkpoint)
        if(!identical(cached$signature,signature))stop("Checkpoint does not match inputs/model")
        ans<-cached$fit;rm(cached)
      }else{
        set.seed(cfg$sgb_seed+i)
        ans<-maaslin3::maaslin3(input_data=as.data.frame(x$relative[,ids,drop=FALSE]),input_metadata=x$design,
          output=bd,formula=formula,normalization="NONE",transform="LOG",standardize=FALSE,
          min_prevalence=0,min_abundance=0,min_variance=0,zero_threshold=0,
          median_comparison_abundance=FALSE,median_comparison_prevalence=FALSE,
          warn_prevalence=FALSE,small_random_effects=FALSE,augment=TRUE,
          plot_summary_plot=FALSE,plot_associations=FALSE,save_models=TRUE,cores=cores,verbosity="INFO")
        # Raw package files retain full batch output; the checkpoint contains
        # only fitted objects/results needed for exact contrasts and median test.
        ans<-lapply(ans[c("fit_data_abundance","fit_data_prevalence")],function(fd)fd[c("results","fits")])
        saveRDS(list(signature=signature,fit=ans),checkpoint,compress=FALSE)
      }
      # Extract exactly the same statistics before releasing this batch's fits.
      # Only small result frames and a checkpoint index survive across batches.
      if(label=="overall"){
        rr<-ans$fit_data_abundance$results
        rr<-rr[rr$name=="fatty_liver",,drop=FALSE]
        index<-setNames(rep(checkpoint,length(ans$fit_data_abundance$fits)),names(ans$fit_data_abundance$fits))
        pieces[[i]]<-list(effects=fl_sgb_extract(ans,"MaAsLin3","Overall","fatty_liver"),
          median_paras=rr,fit_index=index)
      }else{
        pieces[[i]]<-list(effects=fl_sgb_extract(ans,"MaAsLin3","Generation",paste0("case_",c("G0","G1","G2"))),
          heterogeneity=fl_sgb_heterogeneity(ans))
      }
      ledger[[i]]<-data.frame(stage=label,batch=i,first_feature=min(ids),last_feature=max(ids),
        elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),status="COMPLETE")
      fl_write_csv(do.call(rbind,ledger),file.path(out,"audit",paste0(label,"_progress.csv")))
      rm(ans);gc(FALSE)
    }
    fields<-names(pieces[[1]])
    setNames(lapply(fields,function(field){
      values<-lapply(pieces,`[[`,field)
      if(field=="fit_index")do.call(c,values)else do.call(rbind,values)
    }),fields)
  }
  overall<-call_fit("overall",x$formula_overall)
  main<-overall$effects
  fl_write_csv(main,file.path(out,"tables/SGB_overall_components.csv"))
  # Official paper-tag median-comparison routine, restricted to the prespecified
  # disease coefficient. Reuses fitted models; does not refit or re-normalize.
  cat("MaAsLin 3 median-reference abundance sensitivity\n");flush.console()
  rr_med<-overall$median_paras
  # The public package return has renamed individual P columns after joint tests.
  # Restore the internal routine's pre-joint schema, excluding stale joint P/Q.
  rr_med$pval<-rr_med$pval_individual
  rr_med<-rr_med[,setdiff(names(rr_med),c("pval_individual","qval_individual","pval_joint","qval_joint")),drop=FALSE]
  median_fun<-getFromNamespace("run_median_comparison","maaslin3")
  required_median_args<-c("paras","fits","ordereds","groups","metadata","random_effects_formula",
    "median_comparison_threshold","subtract_median","model")
  if(!identical(names(formals(median_fun)),required_median_args))stop("Unverified MaAsLin3 median API")
  set.seed(cfg$sgb_seed+1L)
  disk_fits<-fl_sgb_disk_fits(overall$fit_index)
  med<-median_fun(paras=rr_med,fits=disk_fits,ordereds=character(),groups=character(),
    metadata=x$design,random_effects_formula=~(1|family_id),median_comparison_threshold=0,
    subtract_median=TRUE,model="linear")
  med$median_test_status<-ifelse(is.finite(med$null_hypothesis)&is.finite(med$pval),"ESTIMATED","NOT_ESTIMABLE")
  # The official routine passes pre-existing error rows through unchanged.
  # Their original coefficient/P must not masquerade as a median-reference test.
  med$coef[med$median_test_status!="ESTIMATED"]<-NA_real_
  med$pval[med$median_test_status!="ESTIMATED"]<-NA_real_
  # These are median-reference estimates/tests, not directly measured absolute abundance.
  fl_write_csv(med[,setdiff(names(med),c("qval","qval_joint")),drop=FALSE],
    file.path(out,"tables/SGB_median_reference_sensitivity.csv"))
  rm(overall,disk_fits);gc(FALSE)
  generation<-call_fit("generation",x$formula_generation)
  effects<-generation$effects
  effects$generation<-sub("case_","",effects$term)
  hetero<-generation$heterogeneity
  fl_write_csv(effects,file.path(out,"tables/SGB_generation_effects.csv"))
  fl_write_csv(hetero,file.path(out,"tables/SGB_generation_heterogeneity.csv"))
  rm(generation);gc(FALSE)
  # Prespecified MaAsLin 2 sensitivity: same sample, features and nuisance space.
  # It models all observations with the already fixed log10 pseudocount table.
  if(!requireNamespace("Maaslin2",quietly=TRUE))stop("Maaslin2 sensitivity dependency is absent")
  cat("MaAsLin 2 abundance sensitivity\n");flush.console()
  md<-x$design[,c("fatty_liver",x$mapping$model_column,"family_id"),drop=FALSE]
  m2_pieces<-list();m2_ledger<-list()
  batches<-split(seq_len(ncol(x$log10)),ceiling(seq_len(ncol(x$log10))/20L))
  for(i in seq_along(batches)){
    ids<-batches[[i]];bd<-file.path(out,"private/maaslin2",sprintf("batch_%03d",i))
    dir.create(bd,recursive=TRUE,showWarnings=FALSE)
    checkpoint<-file.path(bd,"complete_results.rds")
    signature<-digest::digest(list(x$log10[,ids,drop=FALSE],md,x$mapping$model_column,
      packageVersion("Maaslin2"),"LM_NONE_NONE_unstandardized_min_abundance_-Inf"),algo="sha256")
    cat(format(Sys.time())," MaAsLin 2 batch ",i,"/",length(batches),"; features ",min(ids),"-",max(ids),"\n",sep="");flush.console()
    started<-Sys.time()
    if(file.exists(checkpoint)){
      cached<-readRDS(checkpoint)
      if(!identical(cached$signature,signature))stop("MaAsLin2 checkpoint does not match inputs/model")
      m2_pieces[[i]]<-cached$results;rm(cached)
    }else{
      m2<-Maaslin2::Maaslin2(input_data=as.data.frame(x$log10[,ids,drop=FALSE]),input_metadata=md,
        output=bd,fixed_effects=c("fatty_liver",x$mapping$model_column),
        random_effects="family_id",normalization="NONE",transform="NONE",standardize=FALSE,
        min_prevalence=0,min_abundance=-Inf,min_variance=0,analysis_method="LM",cores=cores,
        plot_heatmap=FALSE,plot_scatter=FALSE,save_models=TRUE)
      rr_batch<-m2$results;m2_pieces[[i]]<-rr_batch[rr_batch$metadata=="fatty_liver",,drop=FALSE]
      # The official package's saved model files remain in this batch folder.
      saveRDS(list(signature=signature,results=m2_pieces[[i]]),checkpoint)
      rm(m2,rr_batch);gc(FALSE)
    }
    m2_ledger[[i]]<-data.frame(stage="maaslin2",batch=i,first_feature=min(ids),last_feature=max(ids),
      elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),status="COMPLETE")
    fl_write_csv(do.call(rbind,m2_ledger),file.path(out,"audit/maaslin2_progress.csv"))
  }
  rr<-do.call(rbind,m2_pieces)
  fl_write_csv(rr[,setdiff(names(rr),"qval"),drop=FALSE],file.path(out,"tables/SGB_MaAsLin2_sensitivity.csv"))
  list(overall=main,generation=effects,heterogeneity=hetero,sensitivity=rr,median_sensitivity=med)
}

fl_sgb_results <- function(x,r,out) {
  f<-x$features;features<-f$feature_id;wide<-f
  # Leave ambiguous taxonomy in the full table; no arbitrary species assignment.
  for(component in c("prevalence","abundance")){
    z<-r$overall[r$overall$component==component,,drop=FALSE];j<-match(features,z$feature_id)
    for(v in c("estimate","SE","CI_low","CI_high","OR","OR_CI_low","OR_CI_high","fold_change","p_raw","status","package_warning"))
      wide[[paste(component,v,sep="_")]]<-z[[v]][j]
  }
  wide$p_joint_raw<-mapply(fl_sgb_joint,wide$prevalence_p_raw,wide$abundance_p_raw)
  wide$joint_components<-as.integer(is.finite(wide$prevalence_p_raw))+as.integer(is.finite(wide$abundance_p_raw))
  for(g in c("Overall","G0","G1","G2")){
    a<-x$support[x$support$generation==g,];j<-match(features,a$feature_id)
    for(v in c("n","cases","detected","detected_cases","detected_controls","prevalence"))wide[[paste(g,v,sep="_")]]<-a[[v]][j]
  }
  for(component in c("prevalence","abundance")){
    a<-r$heterogeneity[r$heterogeneity$component==component,];j<-match(features,a$feature_id)
    wide[[paste0(component,"_interaction_p_raw")]]<-a$p_interaction_raw[j]
    wide[[paste0(component,"_interaction_status")]]<-a$status[j]
  }
  sm<-r$sensitivity;j<-match(features,sm$feature)
  wide$MaAsLin2_estimate_log10<-sm$coef[j];wide$MaAsLin2_p_raw<-sm$pval[j]
  ms<-r$median_sensitivity;j<-match(features,ms$feature)
  wide$median_reference_effect<-ms$coef[j];wide$median_reference_null<-ms$null_hypothesis[j]
  wide$median_reference_p_raw<-ms$pval[j]
  wide$median_reference_status<-ms$median_test_status[j]
  wide$abundance_direction_concordant<-ifelse(is.finite(wide$abundance_estimate)&is.finite(wide$MaAsLin2_estimate_log10),
    sign(wide$abundance_estimate)==sign(wide$MaAsLin2_estimate_log10),NA)
  # Plot selection is recorded before fitting: 12 smallest joint raw P among
  # valid two-component models with >=10 detected cases and controls.
  eligible<-wide$joint_components==2L & wide$Overall_detected_cases>=10 & wide$Overall_detected_controls>=10
  i<-which(eligible);i<-i[order(wide$p_joint_raw[i],wide$feature_id[i])];i<-head(i,12L)
  wide$shown_in_figure<-seq_len(nrow(wide))%in%i
  fl_write_csv(wide,file.path(out,"tables/SGB_complete_results_EN.csv"))
  fl_write_csv(x$support,file.path(out,"tables/SGB_detection_support.csv"))
  fl_write_csv(wide[i,,drop=FALSE],file.path(out,"tables/SGB_figure_selection.csv"))
  list(complete=wide,selected=wide$feature_id[i])
}

fl_sgb_figure <- function(r,selected,out,cfg) {
  if(!length(selected))stop("No data-supported two-component SGB fits available for the planned figure")
  cols<-c(G0="#F1B656",G1="#397FC7",G2="#040676")
  make_plot<-function(scope,component){
    z<-r[[if(scope=="Overall")"overall"else"generation"]]
    z<-z[z$feature_id%in%selected & z$component==component,,drop=FALSE]
    z$Feature<-factor(sub("^t__","",z$feature_id),levels=rev(sub("^t__","",selected)))
    if(scope=="Overall")z$Group<-factor("Overall")else z$Group<-factor(z$generation,levels=names(cols))
    z$valid<-startsWith(z$status,"ESTIMATED")&is.finite(z$estimate)&is.finite(z$CI_low)&is.finite(z$CI_high)
    pd<-ggplot2::position_dodge(width=if(scope=="Overall")0 else .65)
    if(component=="prevalence"){
      z$effect<-z$OR;z$lo<-z$OR_CI_low;z$hi<-z$OR_CI_high;null<-1
    }else{z$effect<-z$estimate;z$lo<-z$CI_low;z$hi<-z$CI_high;null<-0}
    z$valid<-z$valid&is.finite(z$effect)&is.finite(z$lo)&is.finite(z$hi)
    if(component=="prevalence")z$valid<-z$valid&z$lo>0&z$effect>0&z$hi>0
    a<-ggplot2::ggplot(z,ggplot2::aes(effect,Feature,colour=Group,group=Group))+
      ggplot2::geom_vline(xintercept=null,linetype=2,linewidth=.25,colour="#777777")+
      ggplot2::geom_errorbar(data=z[z$valid,],ggplot2::aes(xmin=lo,xmax=hi),orientation="y",width=0,
        position=pd,linewidth=.30)+
      ggplot2::geom_point(data=z[z$valid,],position=pd,size=1.55,shape=18)+
      ggplot2::scale_colour_manual(values=if(scope=="Overall")c(Overall="#39877F")else cols,
        guide=if(scope=="Overall")"none"else"legend")+
      ggplot2::labs(x=if(component=="prevalence")"Adjusted detection OR (log scale)"else"Adjusted nonzero abundance difference (log2)",
        y=NULL,title=paste(if(scope=="Overall")"Overall"else"By generation",if(component=="prevalence")"detection"else"abundance"))
    if(component=="prevalence")a<-a+ggplot2::scale_x_log10()
    a
  }
  p<-nf_patch(list(make_plot("Overall","prevalence"),make_plot("Overall","abundance"),
    make_plot("Generation","prevalence"),make_plot("Generation","abundance")),"AB\nCD",collect_guides=TRUE)
  spec<-list(name="Fig06_SGB_associations",plot=p,width=183,height=165,panels=4,export_tiff=TRUE)
  nf_export(spec,out,cfg)
  saveRDS(p,file.path(out,"private/figure_plot.rds"))
  invisible(p)
}

fl_sgb_report <- function(x,r,summary,out,cfg) {
  z<-summary$complete
  count<-function(v)sum(is.finite(v)&v<.05)
  stats<-data.frame(metric=c("Participants","Steatosis cases","Primary families","SGBs",
    "Two-component overall tests","Joint raw P < 0.05","Prevalence raw P < 0.05",
    "Abundance raw P < 0.05","Prevalence interaction raw P < 0.05","Abundance interaction raw P < 0.05",
    "Median-reference raw P < 0.05"),
    value=c(nrow(x$metadata),sum(x$metadata$fatty_liver),length(unique(x$metadata$family_id)),nrow(z),
      sum(z$joint_components==2),count(z$p_joint_raw),count(z$prevalence_p_raw),count(z$abundance_p_raw),
      count(z$prevalence_interaction_p_raw),count(z$abundance_interaction_p_raw),count(z$median_reference_p_raw)))
  fl_write_csv(stats,file.path(out,"tables/SGB_analysis_summary.csv"))
  fl_write_csv(as.data.frame(table(component=r$overall$component,status=r$overall$status)),
    file.path(out,"tables/SGB_model_status_summary.csv"))
  notes<-c("# SGB associations with ultrasound-detected hepatic steatosis", "",
    paste0("Eligible cohort: 319 primary families with members in all three generations; analysis sample: ",
      nrow(x$metadata)," participants and ",sum(x$metadata$fatty_liver)," cases. Minors retained. 695 prespecified SGBs."),
    "Primary MaAsLin 3 analysis models detection and nonzero log2 relative abundance. The estimand is association on the relative scale, so the official documented median_comparison_abundance=FALSE setting is used. The original full-feature denominator is preserved (normalization=NONE); no repeat TSS and no pseudocount for the nonzero model.",
    "Family random intercepts are retained. Fixed nuisance adjustment: generation, generation-specific 3-d.f. age splines, sex, survey centre, library batch, sequencing pool, log10 input-matrix total, sampling season sine/cosine, log1p recorded freezing interval, log1p absolute stool-ultrasound interval, recent antibiotics. QR removes only aliases in the nuisance column space. Input-matrix total is not verified raw read depth.",
    "Detection uses the published augmented logistic mixed model. Augmentation is retained exactly as implemented by the paper-tag package; these are conditional model odds ratios, not family-sandwich GEE estimates. The complete-case mean family size is below four: the paper-tag package warns that prevalence mixed-model coefficients may be inflated and P values too small in this setting. This warning is retained and detection results require particular caution, rather than being treated as confirmed discoveries. Singular/convergence diagnostics are also retained.",
    "Generation-specific case coefficients are algebraically equivalent to disease-by-generation interactions with generation main effects included in the nuisance space. Heterogeneity uses a 2-d.f. Satterthwaite F test for abundance and a 2-d.f. Wald test for the augmented detection model. The two components are reported separately. Generation and chronological age are not causally separable by these data alone.",
    "Overall two-component P uses Beta(1,2) CDF(min(component P)), as in MaAsLin 3. If only one component is estimable its P is retained and joint_components=1 explicitly records the fallback. All public P values are raw; package-internal BH values are not used for discovery or figure selection.",
    "95% abundance confidence intervals use the fitted Satterthwaite t distribution; detection intervals use the fitted normal approximation on the log-odds scale. Singular fits remain flagged; failed fits have missing estimates, never zero. Per-feature nonzero support and each generation's case support are provided.",
    "Two prespecified sensitivities address distinct issues: (1) the original MaAsLin 3 paper-tag median-comparison routine, reusing the exact same abundance fits with the disease coefficient tested against its cohort-wide median; (2) MaAsLin 2 on the previously saved log10(relative + fixed pseudocount) matrix, with normalization and transform disabled. These estimands differ from primary nonzero relative abundance. Neither is independent validation; median-reference inference is not an experimental absolute microbial count.",
    "Display rule fixed before fitting: the 12 smallest overall joint raw P among valid two-component models with at least 10 detected cases and 10 detected controls. All 695 SGBs remain in the complete table. Labels retain SGB IDs rather than assigning an arbitrary name when taxonomy is ambiguous. Raw P < 0.05 is exploratory, not a multiplicity-controlled discovery.",
    "No adult MI-SGB analysis is silently run. Existing adult MI was designed for clinical outcomes and needs compatibility review before use for SGB outcomes/interactions. No phenotype or microbial zero is imputed here. No family-specific taxa follow-up, transmission, mediation, prediction or new pathway model is run in this first-stage module.",
    "References: Nickols et al., Nature Methods (2026), doi:10.1038/s41592-025-02923-9; Mallick et al., PLOS Computational Biology (2021), doi:10.1371/journal.pcbi.1009442; Nearing et al., Nature Communications (2022), doi:10.1038/s41467-022-28034-z.")
  writeLines(notes,file.path(out,"reports/Methods_and_figure_legend_EN.md"),useBytes=TRUE)
  zh<-c("# 17 脂肪肝相关SGB及三代差异", "",paste0("主要样本",nrow(x$metadata),"人、",sum(x$metadata$fatty_liver),"例、319个主家系；695个预定SGB，保留未成年人。"),
    "主要模型区分检出关联与非零相对丰度关联，并处理主家系随机截距及既定人口学、技术和采样协变量。",
    "主分析针对相对丰度零原假设；另用MaAsLin 3论文原始中位数比较步骤检查组成参照敏感性。两者的估计目标及P不能混称绝对丰度发现。",
    paste0("原始P<0.05：总体单菌检验",count(z$p_joint_raw),"个；检出分量",count(z$prevalence_p_raw),"个；丰度分量",count(z$abundance_p_raw),"个。均为探索性计数，不是多重检验控制后的发现。"),
    paste0("代际差异原始P<0.05：检出分量",count(z$prevalence_interaction_p_raw),"个；丰度分量",count(z$abundance_interaction_p_raw),"个。G2仅17例，需结合区间、支持数与拟合状态。"),
    "输出完整英文结果表、单菌逐代支持表、代际效应与交互表、MaAsLin2及中位数参照敏感性，以及一张四面板主图。图片沿用既定三代配色、183mm宽、Arial和四种格式。",
    "本模块不运行候选菌家系跟进、不复用临床MI硬套单菌、不运行菌株传播或中介。原始RDS及既有结果不变。")
  writeLines(enc2utf8(zh),file.path(out,"reports/results_zh.md"),useBytes=TRUE)
  esc<-function(v){v<-gsub("&","&amp;",as.character(v),fixed=TRUE);gsub("<","&lt;",v,fixed=TRUE)}
  tbl<-function(dat){paste0("<table><tr>",paste0("<th>",esc(names(dat)),"</th>",collapse=""),"</tr>",
    paste(vapply(seq_len(nrow(dat)),function(i)paste0("<tr>",paste0("<td>",esc(vapply(dat[i,,drop=FALSE],function(a)if(is.numeric(a))format(signif(a,4))else as.character(a),character(1))),"</td>",collapse=""),"</tr>"),character(1)),collapse=""),"</table>")}
  display<-z[match(summary$selected,z$feature_id),c("official_id","Overall_detected_cases","Overall_detected_controls",
    "prevalence_OR","prevalence_p_raw","abundance_estimate","abundance_p_raw","p_joint_raw",
    "prevalence_interaction_p_raw","abundance_interaction_p_raw"),drop=FALSE]
  page<-c("<!doctype html><meta charset='utf-8'><title>SGB associations</title>",
    "<style>body{font:15px/1.65 Arial,sans-serif;color:#333;max-width:1400px;margin:30px auto;padding:0 20px}img{max-width:100%}table{border-collapse:collapse;font-size:12px}th,td{padding:7px;border-bottom:1px solid #ccc;text-align:left}.scroll{overflow:auto}</style>",
    "<h1>SGB associations with hepatic steatosis</h1><p>Exploratory raw P values; fixed three-generation primary-family cohort.</p>",tbl(stats),
    "<p><a href='tables/SGB_complete_results_EN.csv'>Complete English results</a> · <a href='reports/Methods_and_figure_legend_EN.md'>Methods and legend</a></p>",
    "<p><a href='figures/Fig06_SGB_associations.pdf'>PDF</a> · <a href='figures/Fig06_SGB_associations.png'>PNG</a> · <a href='figures/Fig06_SGB_associations.tiff'>TIFF</a> · <a href='figures/Fig06_SGB_associations.svg'>SVG</a></p>",
    "<img src='figures/Fig06_SGB_associations.png'><h2>Prespecified display subset</h2><div class='scroll'>",tbl(display),"</div>")
  writeLines(enc2utf8(page),file.path(out,"index.html"),useBytes=TRUE)
}

fl_sgb_run <- function(cfg) {
  resume<-Sys.getenv("FATTYLIVER_SGB_RESUME_RUN",unset="")
  if(nzchar(resume)){
    if(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$",resume))stop("Unsafe resume run ID")
    module<-file.path(cfg$project_root,"outputs",cfg$output_module)
    out<-file.path(module,"_running",resume)
    if(!dir.exists(out))stop("Resume requires an existing interrupted _running directory")
    marker<-file.path(out,"run_status.txt")
    if(!file.exists(marker))stop("Resume status marker is missing")
    previous_status<-readLines(marker,n=1L,warn=FALSE)
    if(!previous_status%in%c("RUNNING","CANCELLED_BY_USER"))
      stop("Only an interrupted or explicitly cancelled run can be resumed")
    # Supplying the resume run ID explicitly at a later invocation is the user's
    # restart action; merely editing this entry never starts or un-cancels a run.
    if(previous_status=="CANCELLED_BY_USER")file.copy(marker,
      file.path(out,"audit",paste0("cancellation_before_resume_",format(Sys.time(),"%Y%m%d_%H%M%S"),".txt")),overwrite=FALSE)
    cfg$run_id<-resume;ctx<-list(cfg=cfg,module=module,run_dir=out)
  }else{ctx<-fl_output_start_run(cfg);out<-ctx$run_dir}
  for(p in c("tables","figures","reports","private","audit"))dir.create(file.path(out,p),showWarnings=FALSE)
  writeLines("RUNNING",file.path(out,"run_status.txt"))
  tryCatch({
    x<-fl_sgb_prepare(cfg)
    oldmanifest<-file.path(out,"audit/input_manifest.csv")
    if(nzchar(resume)&&file.exists(oldmanifest)){
      old<-read.csv(oldmanifest,stringsAsFactors=FALSE,check.names=FALSE)
      if(!identical(as.character(old$path),as.character(x$source_manifest$path))||
         !identical(as.character(old$md5),as.character(x$source_manifest$md5)))
        stop("Resume input manifest does not match original inputs")
    }
    fl_write_csv(x$source_manifest,file.path(out,"audit/input_manifest.csv"))
    fl_write_csv(x$mapping,file.path(out,"audit/nuisance_design_mapping.csv"))
    fl_write_csv(data.frame(alias=x$aliases),file.path(out,"audit/aliased_nuisance_columns.csv"))
    fl_write_csv(x$knots,file.path(out,"audit/age_spline_knots.csv"))
    fl_write_csv(x$sample_flow,file.path(out,"tables/sample_flow.csv"))
    fl_write_csv(x$sample_counts,file.path(out,"tables/sample_counts_by_generation.csv"))
    saveRDS(x$metadata,file.path(out,"private/model_metadata.rds"))
    saveRDS(x$design,file.path(out,"private/model_design.rds"))
    for(f in c("scripts/17_sgb_association.R","R/sgb_association.R","config/sgb_analysis_20261022.md"))
      file.copy(file.path(cfg$project_root,f),file.path(out,"audit",basename(f)))
    r<-fl_sgb_fit(x,cfg,out);summary<-fl_sgb_results(x,r,out)
    oldpath<-Sys.getenv("PYTHONPATH",unset="");on.exit(Sys.setenv(PYTHONPATH=oldpath),add=TRUE)
    Sys.setenv(PYTHONPATH=paste(c(file.path(cfg$project_root,".figure_qa_python"),if(nzchar(oldpath))oldpath),collapse=.Platform$path.sep))
    oldtheme<-ggplot2::theme_set(nf_theme());on.exit(ggplot2::theme_set(oldtheme),add=TRUE)
    bundle<-file.path(out,"audit/figure_source_bundle.R")
    writeLines(c(readLines(file.path(cfg$project_root,"R/nature_export_approved.R"),encoding="UTF-8"),
      readLines(file.path(cfg$project_root,"R/sgb_association.R"),encoding="UTF-8")),bundle,useBytes=TRUE)
    nf_qa_command(cfg,"validate_figure.py",c(shQuote(bundle),"--json"),file.path(out,"audit/source_preflight.json"))
    fl_sgb_figure(r,summary$selected,out,cfg);fl_sgb_report(x,r,summary,out,cfg)
    stopifnot(identical(unname(tools::md5sum(x$source_manifest$path)),x$source_manifest$md5))
    saveRDS(r,file.path(out,"private/derived_results.rds"))
    capture.output(sessionInfo(),file=file.path(out,"audit/sessionInfo.txt"))
    writeLines("COMPLETED",file.path(out,"run_status.txt"));final<-fl_output_complete_run(ctx)
    key<-file.path(cfg$project_root,"Results_表格与图片/Results_重点表格/tables/English/SGB_associations")
    dir.create(key,recursive=TRUE,showWarnings=FALSE)
    for(f in list.files(file.path(final,"tables"),full.names=TRUE))file.copy(f,key,overwrite=TRUE)
    cat("Completed SGB associations: ",final,"\n",sep="");invisible(final)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(e)})
}
