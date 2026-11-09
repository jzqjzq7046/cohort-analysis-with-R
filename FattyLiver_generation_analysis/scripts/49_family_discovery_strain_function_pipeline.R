.f4_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(.f4_root,'R/family_four_stage.R'),encoding='UTF-8')
f4_pipeline <- function(root=.f4_root,args=character()) {
  stopifnot(all(args%in%c('--run','--resume')))
  if(getRversion()<'4.5') {
    s<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',shQuote(file.path(root,'scripts/49_family_discovery_strain_function_pipeline.R')),vapply(args,shQuote,character(1))))
    if(s!=0)stop('Project-R pipeline failed');return(invisible(NULL))
  }
  log<-file.path(root,'outputs/family_four_stage_20261104');dir.create(log,recursive=TRUE,showWarnings=FALSE)
  result<-list();errors<-character()
  for(k in c('discovery','strain','disease','functionality')) {
    stage_args<-if(k%in%c('discovery','functionality')&&'--resume'%in%args)'--resume'else character()
    result[[k]]<-tryCatch(f4_run(k,args=stage_args,root=root),error=function(e){errors[k]<<-conditionMessage(e);NA_character_})
    saveRDS(list(runs=result,errors=errors),file.path(log,'pipeline_checkpoint.rds'))
  }
  status<-do.call(rbind,lapply(names(result),function(k){p<-result[[k]];data.frame(stage=k,
    status=if(is.na(p))'FAILED'else readLines(file.path(p,'run_status.txt'),1),path=p,
    error=if(k%in%names(errors))errors[[k]]else'')}))
  write.csv(status,file.path(log,'pipeline_status.csv'),row.names=FALSE)
  jsonlite::write_json(list(status=if(length(errors))'FAILED'else'COMPLETED_AVAILABLE_INPUTS',
    all_four_stages_completed=all(status$status=='COMPLETED'),shutdown_requested_by_user=TRUE,
    shutdown_requires_completion_or_user_scope_clarification=TRUE,
    validation_passed=!length(errors),figure_qa_passed=!length(errors),stages=status,
    completed_runs=unname(status$path[status$status=='COMPLETED'])),file.path(log,'completion.json'),auto_unbox=TRUE,pretty=TRUE)
  print(status);if(length(errors))stop('Some stages failed: inspect pipeline_status.csv')
  invisible(status)
}
if(sys.nframe()==0L)f4_pipeline(args=commandArgs(TRUE))else if(isTRUE(getOption('fatty_liver.autorun',TRUE)))f4_pipeline(args=if(isTRUE(getOption('fatty_liver.resume',FALSE)))'--resume'else character())
