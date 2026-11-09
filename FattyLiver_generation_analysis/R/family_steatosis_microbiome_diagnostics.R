# Diagnostic visualizations only; does not refit or change inferential models.
fp_diagnostics <- function(r,cfg,out) {
  audit<-file.path(out,'audit');cache<-file.path(cfg$project_root,'outputs',cfg$output_module,'_checkpoints/private')
  all_dg<-ov_bind(lapply(list.files(audit,pattern='_diagnostics.csv$',full.names=TRUE),tf_read))
  summary<-ov_bind(lapply(split(all_dg,list(all_dg$comparison,all_dg$model),drop=TRUE),function(z)data.frame(
    comparison=z$comparison[1],model=z$model[1],fits=nrow(z),official_matches=sum(z$official_beta_P_match),
    response_max_error=max(z$response_max_error,na.rm=TRUE),singular=sum(z$status=='ESTIMATED_SINGULAR'),
    unreliable=sum(!z$status%in%c('ESTIMATED','ESTIMATED_SINGULAR')),median_residual_skewness=median(z$residual_skewness,na.rm=TRUE),
    residual_absolute_skewness_gt_2=sum(abs(z$residual_skewness)>2,na.rm=TRUE))))
  fl_write_csv(summary,file.path(audit,'Model_diagnostic_summary.csv'))
  chosen<-list();stats<-list();page<-0L
  pdf<-file.path(audit,'B1_residual_assumption_diagnostics.pdf');grDevices::cairo_pdf(pdf,width=8,height=4,family='Arial')
  on.exit(grDevices::dev.off(),add=TRUE)
  for(key in c('G0-G1','G1-G2','G0-G2')){
    dg<-all_dg[all_dg$comparison==key&all_dg$model=='B1'&all_dg$status%in%c('ESTIMATED','ESTIMATED_SINGULAR'),]
    ids<-unique(c(sort(dg$feature_id)[1],dg$feature_id[which.max(abs(dg$residual_skewness))]))
    md<-readRDS(file.path(out,'private',paste0('Linked_younger_',key,'.rds')))
    md<-md[md$fatty_liver==0&v3_cc(md,c(v3_P2,'n_older_observed')),,drop=FALSE]
    cps<-list.files(file.path(cache,paste0(key,'_B1')),pattern='^batch_[0-9]+\\.rds$',recursive=TRUE,full.names=TRUE)
    for(id in ids){fit<-NULL
      for(cp in cps){z<-readRDS(cp)$value;if(!id%in%z$effects$feature_id)next
        i<-as.integer(sub('batch_([0-9]+)\\.rds','\\1',basename(cp)))
        p<-file.path(dirname(cp),sprintf('official_%03d/fit/fits/models.rds',i));ff<-readRDS(p)[[id]]
        if(!is.null(ff)&&identical(sort(rownames(model.frame(ff))),sort(as.character(md$sample_id)))){fit<-ff;break}
      }
      if(is.null(fit))stop('Cannot locate the matching official diagnostic fit: ',key,' ',id)
      rv<-as.numeric(residuals(fit));fv<-as.numeric(fitted(fit));std<-rv/sd(rv);page<-page+1L
      draw<-function(){par(mfrow=c(1,2),mar=c(4,4,3,1),family='Arial',cex=.85)
        qqnorm(std,main=paste(key,sub('^t__','',id),'residual Q-Q'),pch=16,col='#888888',cex=.4);qqline(std,col='black')
        plot(fv,rv,pch=16,col=grDevices::adjustcolor('#777777',.45),cex=.4,xlab='Fitted LOG abundance',ylab='Residual',main='Residual versus fitted');abline(h=0,lty=2)
      }
      draw()
      png<-file.path(audit,sprintf('Residual_diagnostics_%s_%s.png',key,sub('^t__','',id)))
      grDevices::png(png,width=1200,height=600,res=150,type='cairo',bg='white');tryCatch(draw(),finally=grDevices::dev.off())
      sw<-shapiro.test(rv)
      stats[[paste(key,id)]]<-data.frame(comparison=key,model='B1',feature_id=id,N=length(rv),
        selection=if(id==sort(dg$feature_id)[1])'First feature ID, independent of significance'else'Maximum absolute residual skewness',
        Shapiro_W=unname(sw$statistic),Shapiro_P=sw$p.value,diagnostic_page=page,
        interpretation='Diagnostic only; a small normality P does not by itself change the prespecified published LOG model')
    }
  }
  fl_write_csv(ov_bind(stats),file.path(audit,'Residual_diagnostic_checks.csv'))
  writeLines(c('# Numerical and assumption review',
    'Official package coefficients/P values and fitted transformed responses are checked for every valid model. Optimizer and singularity states are retained in complete tables.',
    'Residual Q-Q and residual-versus-fitted plots cover a deterministic first SGB and the maximum absolute residual-skewness SGB in each all-age comparison. These are diagnostics, not additional discovery tests. Zero-rich LOG outcomes can show nonnormal residuals or unequal spread. The adopted published model is retained; diagnostics are not proof of universal Gaussian adequacy or causal identification.',
    'All adult paired-sample comparisons use exactly the same participant IDs. Primary-family labels remain grouping variables rather than verified pedigrees.'),file.path(audit,'Numerical_and_assumption_review.md'))
  invisible(summary)
}

fp_scope_table <- function(out) {
  pop<-readRDS(file.path(out,'private/analysis_population.rds'))
  z<-ov_bind(lapply(v3_pairs,function(g){key<-paste(g,collapse='-')
    md<-readRDS(file.path(out,'private',paste0('Linked_younger_',key,'.rds')))
    md<-md[md$fatty_liver==0&v3_cc(md,c(v3_P2,'n_older_observed')),,drop=FALSE]
    older<-pop$clinical[pop$clinical$generation==g[1],,drop=FALSE]
    older_stool_families<-unique(older$family_id[older$participant_id%in%pop$metadata$participant_id])
    data.frame(comparison=key,unaffected_younger_B1_N=nrow(md),with_at_least_one_older_stool_N=sum(md$family_id%in%older_stool_families),
      no_observed_older_stool_N=sum(!md$family_id%in%older_stool_families),
      added_beyond_prior_two_microbiome_generation_baseline_N=sum(!md$sample_id%in%pop$prior_microbiome_ids))
  }))
  v3_table(z,'00_Clinical_versus_microbiome_linkage_scope_EN',out);z
}

fp_validate <- function(r,cfg,out) {
  t<-r$results;checks<-list();add<-function(label,ok,detail=''){
    checks[[length(checks)+1L]]<<-data.frame(check=label,status=if(ok)'PASS'else'FAIL',detail=detail)
    if(!ok)stop('Module 27 validation failed: ',label)
  }
  add('One result per feature/comparison/model',!anyDuplicated(t[,c('feature_id','comparison','model')]))
  fp_scope_table(out)
  for(model in unique(t$model)){
    k<-which(t$model==model&t$eligible);good<-tf_bool(t$analysis_valid[k])&is.finite(t$P_raw[k]);p<-rep(1,length(k));p[good]<-t$P_raw[k][good]
    q<-p.adjust(p,'BH');add(paste(model,'joint BH independently recomputed'),max(abs(q[good]-t$P_BH[k[good]]))<1e-12)
    add(paste(model,'joint eligible test count'),all(t$FDR_tests[t$model==model]==length(k)))
  }
  valid<-tf_bool(t$analysis_valid)
  add('All valid effect intervals and P values finite',all(is.finite(as.matrix(t[valid,c('beta_log2','SE','CI_low','CI_high','P_raw','P_BH')]))))
  add('Confidence intervals contain estimates',all(t$CI_low[valid]<=t$beta_log2[valid]&t$CI_high[valid]>=t$beta_log2[valid]))
  add('Abundance ratio transform',max(abs(log2(t$abundance_ratio[valid])-t$beta_log2[valid]))<1e-10)
  for(key in unique(t$comparison)){
    a<-readRDS(file.path(out,'private',paste0(key,'_Adult_B1_same_sample_design.rds')))$metadata
    b<-readRDS(file.path(out,'private',paste0(key,'_B2_design.rds')))$metadata
    add(paste(key,'adult B1/B2 identical participant IDs'),identical(a$sample_id,b$sample_id))
    for(model in c('B1','Adult_B1_same_sample','B2')){
      d<-readRDS(file.path(out,'private',paste0(key,'_',model,'_design.rds')))
      add(paste(key,model,'unaffected younger, unique individual rows'),all(d$metadata$fatty_liver==0)&&!anyDuplicated(d$metadata$sample_id))
      add(paste(key,model,'older exposure definition'),all(d$metadata$older_any_case==as.integer(d$metadata$n_older_cases>0))&&all(d$metadata$n_older_observed>=1))
      z<-t[t$comparison==key&t$model==model,]
      add(paste(key,model,'result N matches design'),all(z$N==nrow(d$metadata))&&all(z$exposed_N==sum(d$metadata$older_any_case)))
      add(paste(key,model,'family random intercept eligibility'),identical(d$random,any(table(d$metadata$family_id)>1)))
    }
  }
  files<-tf_read(file.path(out,'audit/input_fingerprints.csv'))
  add('Original source fingerprints unchanged',identical(as.character(files$md5),unname(tools::md5sum(files$path))))
  scripts<-file.path(cfg$project_root,c('scripts/27_family_steatosis_microbiome.R','R/family_steatosis_microbiome.R',
    'R/family_steatosis_microbiome_figures.R','R/family_steatosis_microbiome_diagnostics.R','config.R'))
  fl_write_csv(data.frame(path=scripts,md5=unname(tools::md5sum(scripts))),file.path(out,'audit/analysis_source_fingerprints.csv'))
  fl_write_csv(ov_bind(checks),file.path(out,'audit/Validation_checks.csv'))
  # Validate the actual source together with its imported export/theme helpers.
  bundle<-file.path(out,'audit/Plotting_source_bundle.R')
  writeLines(unlist(lapply(file.path(cfg$project_root,c('R/family_steatosis_microbiome_figures.R','R/taxa_followup_figures.R',
    'R/nature_export_approved.R')),readLines,warn=FALSE)),bundle,useBytes=TRUE)
  log<-file.path(out,'audit/Plotting_source_static_preflight.txt')
  status<-system2(cfg$python,c(shQuote(file.path(nf_skill_root(),'scripts/validate_figure.py')),shQuote(bundle),'--backend','r'),stdout=log,stderr=log)
  if(status!=0)stop('Combined plotting-source preflight requires review; see ',log)
  invisible(ov_bind(checks))
}
