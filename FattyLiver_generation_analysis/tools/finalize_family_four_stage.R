# Reconcile already-saved analyses after rendered QA. No model fit/permutation.
root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(root,'R/family_four_stage.R'),encoding='UTF-8');cfg<-f4_setup(root)
keys<-c(discovery='family_discovery',functionality='family_function');runs<-list();res<-list();qa<-list()
for(k in names(keys)) {
  mod<-file.path(root,'outputs',unname(cfg$output_modules[keys[[k]]]))
  p<-file.path(mod,'current',readLines(file.path(mod,'LATEST_RUN.txt'),1))
  stopifnot(readLines(file.path(p,'run_status.txt'),1)=='COMPLETED')
  r<-readRDS(file.path(p,'private/analysis_result.rds'));stopifnot(identical(r$protocol,f4_protocol))
  f4_validate(r,p);ade_protect(r$protected)
  stems<-sub('\\.pdf$','',list.files(file.path(p,'figures'),pattern='\\.pdf$'))
  stopifnot(length(stems)>0)
  for(stem in stems) {
    fs<-file.path(p,'figures',paste0(stem,c('.pdf','.png','.svg','.tiff')))
    stopifnot(all(file.exists(fs)),all(file.info(fs)$size>100))
    a<-jsonlite::fromJSON(file.path(p,'audit',paste0(stem,'.effective-text-audit.json')))
    b<-jsonlite::fromJSON(file.path(p,'audit',paste0(stem,'.collision-audit.json')))
    al<-jsonlite::fromJSON(file.path(p,'audit',paste0(stem,'.alignment.json')))
    stopifnot(a$verdict=='PASS',b$verdict=='PASS',b$summary$fail==0,b$summary$warn==0)
    stopifnot(al$verdict%in%c('PASS','NOT APPLICABLE'),al$summary$fail==0,al$summary$warn==0)
    qa[[paste(k,stem)]]<-data.frame(stage=k,figure=stem,formats='PDF;PNG;SVG;TIFF',minimum_pt=a$minimum_pt,
      collision_fail=b$summary$fail,collision_warn=b$summary$warn,alignment_verdict=al$verdict,raw_hashes=paste(unname(tools::md5sum(fs)),collapse=';'))
  }
  static<-readLines(file.path(p,'audit/figure_static_preflight.txt'),warn=FALSE)
  stopifnot(!any(grepl('[FAIL]',static,fixed=TRUE)))
  runs[[k]]<-p;res[[k]]<-r
}
log<-file.path(root,'outputs/family_four_stage_20261104');dir.create(log,recursive=TRUE,showWarnings=FALSE)
statuses<-lapply(c('discovery','strain','disease','functionality'),function(k){
  if(k%in%names(runs))return(data.frame(stage=k,status='COMPLETED',path=runs[[k]],error=''))
  mod<-file.path(root,'outputs',unname(cfg$output_modules[if(k=='strain')'family_strain'else'family_strain_disease']))
  paths<-list.dirs(file.path(mod,'_running'),recursive=FALSE,full.names=TRUE)
  paths<-paths[file.exists(file.path(paths,'run_status.txt'))];stopifnot(length(paths)>0)
  p<-paths[which.max(file.info(file.path(paths,'run_status.txt'))$mtime)]
  status<-readLines(file.path(p,'run_status.txt'),1);stopifnot(status=='BLOCKED_INPUT')
  data.frame(stage=k,status=status,path=p,error='Sequence-derived strain input is absent')
})
s<-do.call(rbind,statuses)
for(f in c('pipeline_status.csv','completion.json')) {
  old<-file.path(log,f);if(file.exists(old)&&!file.exists(file.path(log,paste0('initial_',f))))
    stopifnot(file.copy(old,file.path(log,paste0('initial_',f))))
}
write.csv(s,file.path(log,'pipeline_status.csv'),row.names=FALSE)
write.csv(do.call(rbind,qa),file.path(log,'final_figure_verification.csv'),row.names=FALSE)
record<-list(status='COMPLETED_AVAILABLE_INPUTS',all_four_stages_completed=FALSE,
  validation_passed=TRUE,figure_qa_passed=TRUE,shutdown_requested_by_user=TRUE,
  shutdown_condition='Full four-stage completion; clarification requested because two stages are blocked',
  stages=s,completed_runs=unname(unlist(runs)),checked_at=format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'))
jsonlite::write_json(record,file.path(log,'completion.json'),auto_unbox=TRUE,pretty=TRUE)
ss<-res$discovery$summary;ff<-res$functionality$summary
ix<-match(ss$model,ff$model);stopifnot(!anyNA(ix),all(ss$N==ff$N[ix]),all(ss$cases==ff$cases[ix]))
overview<-data.frame(Model=ss$model_display,Participants=ss$N,Steatosis_cases=ss$cases,Families=ss$families,
  Tested_SGBs=ss$tested,BH_supported_SGBs=ss$BH_lt_005,Tested_pathways=ff$tested[ix],BH_supported_pathways=ff$BH_lt_005[ix],
  SGB_global_omnibus_q=ss$global_P_BH_across_three_models,Pathway_global_omnibus_q=ff$global_P_BH_across_three_models[ix])
overview<-overview[match(c('All-age Model 3','Adult Model 3','Adult Model 4'),overview$Model),]
write.csv(overview,file.path(log,'Comparison_of_prespecified_models_EN.csv'),row.names=FALSE)
fl_write_csv(overview,file.path(runs$discovery,'tables/Comparison_of_prespecified_models_EN.csv'))
lines<-c('# Four-stage analysis: saved results and unresolved inputs',
  '',paste('Protocol:',f4_protocol),'',
  '1. De novo within-family SGB discovery: COMPLETED.',
  '2. Sequence-based strain typing and cross-generation sharing: BLOCKED_INPUT.',
  '3. True strain sharing versus both-member steatosis status: BLOCKED_INPUT.',
  '4. Within-family HUMAnN functional-potential associations: COMPLETED. This is not strain-specific functional attribution.',
  '', 'Adult Model 4 is primary. Both adult models have identical samples. All-age Model 3 is a separate sensitivity population.',
  '', '| Model | Participants | Families | SGBs tested / BH supported | Pathways tested / BH supported |',
  '|---|---:|---:|---:|---:|',
  vapply(seq_len(nrow(overview)),function(i){z<-overview[i,];sprintf('| %s | %d | %d | %d / %d | %d / %d |',
    z$Model,z$Participants,z$Families,z$Tested_SGBs,z$BH_supported_SGBs,z$Tested_pathways,z$BH_supported_pathways)},character(1)),
  '', 'Analysis populations overlap and counts are not additive. Adult adjustment comparisons use exactly the same people. Attenuation cannot be attributed to BMI alone because lifestyle variables are added simultaneously.',
  '', 'Primary SGB discoveries:',capture.output(print(res$discovery$effects[res$discovery$effects$model=='B2'&res$discovery$effects$BH_lt_005,
    c('feature_id','feature_label','beta_percentage_points','P_raw','P_BH')],row.names=FALSE)),
  '', 'Complete official-LDM statistics, conventional BH q, official LDM q, sample flow, covariate design and source annotations are retained. Raw inputs and previous analyses are unchanged.',
  'All current figures passed rendered text/collision QA and exist as PDF, SVG, PNG and TIFF. Four-format hashes are saved in final_figure_verification.csv.',
  'No discovery is labelled inherited, genetically determined or transmitted. Primary SGBs without exact HUMAnN species mapping receive no forced functional attribution.',
  'Shutdown was requested after completion. Because two strain-dependent stages cannot be completed, a clarification was requested; no affirmative answer is assumed from elapsed time.')
writeLines(lines,file.path(log,'Results_and_stage_status_EN.md'),useBytes=TRUE)
for(p in runs['discovery']) {
  file.copy(file.path(log,'Results_and_stage_status_EN.md'),file.path(p,'reports/Pipeline_stage_status_EN.md'),overwrite=TRUE)
  file.copy(file.path(root,'analysis_plan/strain_input_inventory_20261104.md'),file.path(p,'reports/Strain_input_inventory_EN.md'),overwrite=TRUE)
}
source(file.path(root,'R/result_gallery.R'),local=TRUE,encoding='UTF-8');fl_refresh_result_gallery(root)
print(s);cat('VERIFIED AVAILABLE WORK; strain-dependent stages remain blocked.\n')
