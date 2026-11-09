# Complete pathway association -> function -> observed sources -> generations.
# Definitions are inert. Model fits use the official, audited project helpers.
pfs_protocol <- 'functional_profiles_sources_and_generation_20261102_v1'
pfs_setup <- function(root) {
  sys.source(file.path(root,'R/generation_ADE_extensions.R'),envir=.GlobalEnv)
  cfg <- ade_setup(root)
  sys.source(file.path(root,'R/generation_ADE_figures.R'),envir=.GlobalEnv)
  sys.source(file.path(root,'R/pathway_adult81_heatmaps.R'),envir=.GlobalEnv)
  cfg$output_module <- unname(cfg$output_modules['pathway_redesign'])
  cfg$run_id <- unname(cfg$run_ids['pathway_redesign']);cfg$run_id_explicit <- FALSE
  cfg$v3_permutations <- 9999L;cfg$v3_seed <- 20261102L
  cfg
}
pfs_save <- function(r,out)saveRDS(r,file.path(out,'private/analysis_result.rds'),compress='gzip')
pfs_recovery <- function(module,out) {
  dirs <- unlist(lapply(c('_running','history','current'),function(k)
    list.dirs(file.path(module,k),recursive=FALSE,full.names=TRUE)),use.names=FALSE)
  paths <- file.path(setdiff(dirs,out),'private/analysis_result.rds');paths<-paths[file.exists(paths)]
  if(!length(paths))stop('No completed module-44 statistical checkpoint to resume')
  path <- paths[which.max(file.info(paths)$mtime)];r<-readRDS(path)
  stopifnot(identical(r$display_protocol,pfs_protocol))
  previous<-dirname(dirname(path))
  for(folder in c('tables','reports','audit','figures')){
    fs<-list.files(file.path(previous,folder),full.names=TRUE)
    if(length(fs)&&!all(file.copy(fs,file.path(out,folder),recursive=TRUE,overwrite=TRUE)))stop('Cannot recover ',folder)
  }
  for(name in c('functional_input.rds','run_config.rds','source_standalone.rds')){
    p<-file.path(previous,'private',name)
    if(file.exists(p)&&!file.copy(p,file.path(out,'private',name),overwrite=TRUE))stop('Cannot recover ',name)
  }
  checkpoints<-list.files(file.path(previous,'private'),pattern='^composition_.*\\.rds$',full.names=TRUE)
  if(length(checkpoints)&&!all(file.copy(checkpoints,file.path(out,'private'),overwrite=TRUE)))
    stop('Cannot recover completed community checkpoints')
  r$recovered_from <- previous;r
}
pfs_main <- function(args=character(),root=Sys.getenv('FATTYLIVER_ROOT',
  unset='E:/YFS/workspace/FattyLiver_generation_analysis')) {
  if(getRversion()<'4.5'){
    s<-system2(file.path(root,'tools/R-4.5.2/bin/x64/Rscript.exe'),c('--vanilla',
      shQuote(file.path(root,'scripts/44_pathway_function_sources.R')),vapply(args,shQuote,character(1))))
    if(s!=0L)stop('Project-R pathway redesign failed');return(invisible(NULL))
  }
  stopifnot(length(args)<=1L,all(args%in%c('--run','--statistics-only','--resume','--figures-only','--help')))
  cfg<-pfs_setup(root)
  if('--help'%in%args){cat('44_pathway_function_sources.R [--run | --statistics-only | --resume | --figures-only]\n');return(invisible(NULL))}
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(k in c('private','audit','tables','figures','reports'))dir.create(file.path(out,k),showWarnings=FALSE)
  writeLines(c('RUNNING',pfs_protocol),file.path(out,'run_status.txt'))
  logold<-getOption('fatty_liver.ade_log');options(fatty_liver.ade_log=file.path(out,'audit/execution.log'))
  on.exit(options(fatty_liver.ade_log=logold),add=TRUE)
  tryCatch({
    ade_log('START pathway redesign ',cfg$run_id)
    src<-tf_source(root)
    if(any(c('--resume','--figures-only')%in%args)){
      r<-pfs_recovery(ctx$module,out);ade_protect(r$protected)
      ade_log('Recovered complete statistical estimates without fitting: ',r$recovered_from)
    }else{
      x<-tf_load(src);protected<-ade_fingerprint(src$files)
      imp<-ade_import_paths(cfg,x,out);x<-imp$x
      protected<-ov_bind(list(protected,imp$protected))
      cache<-file.path(ctx$module,'_checkpoints',pfs_protocol,'private');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
      # New module namespace: the first run cannot reuse any module-42 fit.
      r<-ade_analyze('D',x,cfg,out,cache,FALSE);r$protected<-protected
      r$display_protocol<-pfs_protocol;r$input_availability<-imp$coverage
      r$source_module24_run<-src$run_id;r$completed_stages<-'official_pathway_models'
      pfs_save(r,out);rm(x,imp);gc(FALSE)
    }
    saveRDS(cfg,file.path(out,'private/run_config.rds'));pfs_save(r,out)
    if('--statistics-only'%in%args){
      writeLines(c('STATISTICS_READY',pfs_protocol,'Staged; not adopted or published'),file.path(out,'run_status.txt'))
      ade_log('STATISTICS_READY: ',out);return(invisible(out))
    }
    # Extension functions are loaded after completed fits; recovery never refits.
    sys.source(file.path(root,'R/pathway_function_sources_details.R'),envir=.GlobalEnv)
    r<-pfs_details(r,cfg,out,'--figures-only'%in%args);pfs_save(r,out)
    sys.source(file.path(root,'R/pathway_function_sources_figures.R'),envir=.GlobalEnv)
    r<-pfs_prepare_adult81(r,cfg,out);pfs_save(r,out)
    pfs_reports(r,cfg,out);pfs_figures(r,cfg,out);pfs_validate(r,cfg,out)
    ade_protect(r$protected)
    writeLines(capture.output(sessionInfo()),file.path(out,'audit/sessionInfo.txt'))
    writeLines(c('COMPLETED',pfs_protocol),file.path(out,'run_status.txt'))
    target<-fl_output_complete_run(ctx)
    writeLines(basename(target),file.path(root,'outputs',cfg$output_module,'LATEST_RUN.txt'))
    options(fatty_liver.ade_log=file.path(target,'audit/execution.log'))
    ade_log('COMPLETED: ',target);invisible(target)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(conditionMessage(e),call.=FALSE)})
}
