.fc_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.fc_root,'R/family_candidate_followup.R'),encoding='UTF-8')
fc_revision_v3<-function(root=.fc_root) {
  if(getRversion()<'4.5') {
    status<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts/50c_revise_model_evidence_restore_boxes.R'))))
    if(status!=0)stop('Project-R figure revision failed');return(invisible(NULL))
  }
  sys.source(file.path(root,'R/family_four_stage.R'),envir=.GlobalEnv);cfg<-f4_setup(root)
  parent<-fc_adopted(root,'50_Family_candidate_characterization')
  r<-readRDS(file.path(parent,'private/analysis_result.rds'));stopifnot(identical(r$protocol,fc_protocol))
  frozen<-c(file.path(parent,'private/analysis_result.rds'),file.path(parent,'tables',
    c('All_selected_candidates_EN.csv','Saved_model_evidence_EN.csv','Host_tests_EN.csv','Host_display_EN.csv')),
    file.path(parent,'figures',paste0('Figure50_Candidate_host_phenotypes_v2.',c('pdf','png','tiff','svg'))))
  frozen_hashes<-unname(tools::md5sum(frozen));stopifnot(!anyNA(frozen_hashes))
  cfg$output_module<-'50_Family_candidate_characterization';cfg$run_id<-'family_candidates_display_20261104_v3';cfg$run_id_explicit<-FALSE
  ctx<-fl_output_start_run(cfg);out<-ctx$run_dir
  tryCatch({
    for(d in c('tables','reports','audit','private','figures')) {
      for(file in list.files(file.path(parent,d),recursive=TRUE,full.names=TRUE,all.files=TRUE)) {
        relative<-substring(file,nchar(parent)+2L);destination<-file.path(out,relative)
        dir.create(dirname(destination),recursive=TRUE,showWarnings=FALSE)
        stopifnot(file.copy(file,destination,overwrite=TRUE))
      }
    }
    # Retire replaced viewing payloads privately; the complete former run remains in history.
    previous<-file.path(out,'private/previous_display');dir.create(previous,recursive=TRUE,showWarnings=FALSE)
    files<-list.files(file.path(out,'figures'),full.names=TRUE)
    replaced<-files[grepl('^Figure50_(Candidate_model_evidence|Family_abundance_profiles)',basename(files))]
    for(file in replaced)stopifnot(file.rename(file,file.path(previous,basename(file))))
    for(folder in c('tables','reports')) {
      files<-list.files(file.path(out,folder),full.names=TRUE)
      old<-files[grepl('^(Family_paired_difference_estimates_EN[.]|Display_revision_EN[.])',basename(files))]
      for(file in old)stopifnot(file.rename(file,file.path(previous,basename(file))))
    }
    manifest<-tf_read(file.path(out,'tables/Figure_export_manifest_EN.csv'))
    manifest<-manifest[!grepl('^Figure50_(Candidate_model_evidence|Family_abundance_profiles)',manifest$figure),]
    v3_table(manifest,'Figure_export_manifest_EN',out)
    fc_family_data(fc_adopted(root,'45_Family_de_novo_discovery'),r$candidates$feature_id,out)
    fc_validate(r,out);fc_reports(r,out);v3_table(fc_refs(),'References_and_adaptations_EN',out)
    fc_figures(r,cfg,out,redraw=c('evidence','family'))
    relative<-substring(frozen,nchar(parent)+2L)
    new_hashes<-unname(tools::md5sum(file.path(out,relative)))
    stopifnot(identical(new_hashes,frozen_hashes))
    fl_write_csv(data.frame(file=relative,unchanged=TRUE,md5=frozen_hashes),file.path(out,'audit/unchanged_inference_and_white_heatmap.csv'))
    writeLines(c('Corrected figure mapping: model-evidence image is the requested Figure 2; family-abundance image is Figure 3.',
      'Model evidence redesigned as effect-value heatmaps with original BH significance annotations, after Chen et al. 2024 Figure 6A.',
      'Left LOG abundance ratios and right within-family LDM percentage-point effects retain distinct colour scales and units; no pooled effect-size interpretation.',
      'Family abundance restored to connected state means with coloured boxplots and grey points; all 731 discovery participants and 289 families are retained.',
      'White-centered host heatmap, saved analysis object, host associations and candidate/model-evidence numerical tables are byte-identical to the previous completed run.',
      'No model refit, permutation, new bootstrap, P/q recalculation or reselection occurred.'),file.path(out,'reports/Display_revision_v3_EN.txt'))
    ade_protect(r$protected);writeLines(c('COMPLETED',fc_protocol),file.path(out,'run_status.txt'))
    target<-fl_output_complete_run(ctx);writeLines(basename(target),file.path(ctx$module,'LATEST_RUN.txt'))
    cat('COMPLETED:',target,'\n');invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
if(sys.nframe()==0L)fc_revision_v3()else if(isTRUE(getOption('fatty_liver.autorun',TRUE)))fc_revision_v3()
