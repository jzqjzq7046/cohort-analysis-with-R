.fc_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.fc_root,'R/family_candidate_followup.R'),encoding='UTF-8')
fc_redraw<-function(root=.fc_root) {
  if(getRversion()<'4.5') {
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts/50b_redraw_family_candidates.R'))))
    if(status!=0)stop('Project-R display revision failed');return(invisible(NULL))
  }
  sys.source(file.path(root,'R/family_four_stage.R'),envir=.GlobalEnv);cfg<-f4_setup(root)
  parent<-fc_adopted(root,'50_Family_candidate_characterization')
  r<-readRDS(file.path(parent,'private/analysis_result.rds'));stopifnot(identical(r$protocol,fc_protocol))
  protected<-data.frame(path=file.path(parent,'figures',paste0('Figure50_Candidate_model_evidence.',c('pdf','png','tiff','svg'))))
  protected$md5<-unname(tools::md5sum(protected$path));stopifnot(!anyNA(protected$md5))
  cfg$output_module<-'50_Family_candidate_characterization';cfg$run_id<-'family_candidates_display_20261104_v2';cfg$run_id_explicit<-FALSE
  ctx<-fl_output_start_run(cfg);out<-ctx$run_dir
  tryCatch({
    for(d in c('tables','reports','audit','private','figures')) {
      files<-list.files(file.path(parent,d),recursive=TRUE,full.names=TRUE,all.files=TRUE)
      for(file in files) {
        relative<-substring(file,nchar(parent)+2L);destination<-file.path(out,relative)
        dir.create(dirname(destination),recursive=TRUE,showWarnings=FALSE)
        stopifnot(file.copy(file,destination,overwrite=TRUE))
      }
    }
    fc_family_data(fc_adopted(root,'45_Family_de_novo_discovery'),r$candidates$feature_id,out)
    fc_validate(r,out);fc_reports(r,out);v3_table(fc_refs(),'References_and_adaptations_EN',out)
    fc_figures(r,cfg,out,redraw=c('family','host'))
    stopifnot(identical(unname(tools::md5sum(file.path(out,'figures',basename(protected$path)))),protected$md5),
      identical(unname(tools::md5sum(file.path(out,'private/analysis_result.rds'))),unname(tools::md5sum(file.path(parent,'private/analysis_result.rds')))))
    fl_write_csv(protected,file.path(out,'audit/unchanged_model_evidence_source.csv'))
    writeLines(c('Display-only revision: white heatmap midpoint; family paired estimation plot.',
      'All saved disease models, correlations, P/q, people, feature selection and model-evidence figure bytes are unchanged.',
      'Additional bootstrap intervals describe equal-family mean raw differences; they do not replace adjusted discovery inference.',
      'Display-input correction: the previous stage-24 plotting matrix lacked 42 discovery samples and displayed 268 families. The revised family plot uses the exact protected discovery matrix, all 731 participants and all 289 primary families, with no missing sample matches. Saved discovery and host-correlation results are unchanged.'),
      file.path(out,'reports/Display_revision_EN.txt'))
    ade_protect(r$protected);writeLines(c('COMPLETED',fc_protocol),file.path(out,'run_status.txt'))
    target<-fl_output_complete_run(ctx);writeLines(basename(target),file.path(ctx$module,'LATEST_RUN.txt'))
    cat('COMPLETED:',target,'\n');invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
if(sys.nframe()==0L)fc_redraw()else if(isTRUE(getOption('fatty_liver.autorun',TRUE)))fc_redraw()
