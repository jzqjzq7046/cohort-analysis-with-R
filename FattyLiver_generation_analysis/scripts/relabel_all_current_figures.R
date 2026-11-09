invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
source(file.path(root,'R/literature_ABC_20261028.R'),encoding='UTF-8');cfg<-lr_setup(root)
source(file.path(root,'R/external_health_ranking.R'),encoding='UTF-8')
source(file.path(root,'R/generation_reference_phylogeny.R'),encoding='UTF-8')
audit<-file.path(root,'inventory/model_labels_20261030')
targets<-jsonlite::read_json(file.path(audit,'targets.json'),simplifyVector=TRUE)
# This script never calls any statistical fitting entry. Functions only draw
# saved data; tables/reports they generate stay in a private staging directory.
renderer_signature<-unname(tools::md5sum(file.path(root,'R/model_display_labels.R')))
allowed<-character();.nf_export<-nf_export;.tf_export<-tf_export;.ov_export<-ov_export
cached<-function(stem,out){p<-file.path(out,'audit/model_display_signature.txt');file.exists(p)&&identical(readLines(p,1),renderer_signature)&&all(file.exists(file.path(out,'figures',paste0(stem,c('.pdf','.png','.tiff','.svg')))))}
empty_export<-function()data.frame()
nf_export<-function(spec,out,cfg){if(!spec$name%in%allowed)return(empty_export());if(cached(spec$name,out))return(empty_export());.nf_export(spec,out,cfg)}
tf_export<-function(spec,stem,out,cfg){if(!stem%in%allowed)return(empty_export());if(cached(stem,out))return(empty_export());.tf_export(spec,stem,out,cfg)}
ov_export<-function(plot,stem,out,cfg,width=183,height=235){if(!stem%in%allowed)return(empty_export());if(cached(stem,out))return(empty_export());.ov_export(plot,stem,out,cfg,width,height)}
for(module in unique(targets$module)){
 if('--module24-only'%in%commandArgs(TRUE)&&module!='24_Two_generation_reanalysis')next
 from<-lr_adopted(root,module);out<-file.path(audit,'staging',module);allowed<-targets$stem[targets$module==module]
 for(folder in c('figures','tables','reports','audit','private'))dir.create(file.path(out,folder),recursive=TRUE,showWarnings=FALSE)
 tf_log('RELABEL SAVED FIGURES: ',module)
 if(module=='24_Two_generation_reanalysis'){
  # Native grid cohort diagram is handled separately after all ggplot exports.
  allowed<-setdiff(allowed,'Study_population_and_analysis_sets')
  for(f in list.files(file.path(from,'private'),pattern='^clinical_pairs_.*[.]rds$',full.names=TRUE))file.copy(f,file.path(out,'private'),overwrite=TRUE)
  load_stage<-function(n)readRDS(file.path(from,'private',paste0('stage_',n,'.rds')))
  v3_visualize(load_stage('prepare'),load_stage('clinical'),load_stage('alpha'),load_stage('similarity'),load_stage('composition'),load_stage('taxa'),cfg,out,only=allowed)
 }else{
  file<-if(module=='35_Generation_reference_phylogeny')'figure_data.rds'else'analysis_result.rds'
  r<-readRDS(file.path(from,'private',file))
  if(module=='25_Taxa_robustness')tf_figures25(r,cfg,out)
  if(module=='29_Within_family_steatosis')abc_figures(r,'B',cfg,out)
  if(module=='30_External_health_sets')ehr_figures(r,cfg,out)
  if(module=='32_Within_family_literature')lr_A_figures(r,cfg,out)
  if(module=='33_Steatosis_severity'){
    stopifnot(file.copy(file.path(from,'private/B_Grade_distribution_data.rds'),file.path(out,'private/B_Grade_distribution_data.rds'),overwrite=TRUE));lr_B_figures(r,cfg,out)
  }
  if(module=='35_Generation_reference_phylogeny'){
    lr_export(list(gp_detail(r,cfg)),'Supplement35_Thirty_most_abundant_SGBs',out,cfg,width=183,height=210)
  }
 }
 writeLines(renderer_signature,file.path(out,'audit/model_display_signature.txt'))
 tf_log('FINISHED STAGED: ',module)
}
# Native grid drawing reads only source counts; preserve its existing legend.
cohort_source<-file.path(root,'scripts/redraw_study_population.R')
cohort_code<-parse(cohort_source,encoding='UTF-8')
# Do not evaluate the concluding caption-writing expression in this task.
cohort_code<-Filter(function(e)!grepl('writeLines.*Study population',paste(deparse(e),collapse=' ')),as.list(cohort_code))
# Source script already uses established counts/layout. Its final writeLines is
# identified structurally rather than rewriting the caption text.
cohort_code<-Filter(function(e)!(is.call(e)&&identical(e[[1]],as.name('writeLines'))),cohort_code)
for(e in cohort_code)eval(e,envir=.GlobalEnv)
cfg<-fl_config(root)
cohort_out<-file.path(lr_adopted(root,'24_Two_generation_reanalysis'),'figures/Cohort_overview_two_generation_20261024')
cohort_audit<-file.path(lr_adopted(root,'24_Two_generation_reanalysis'),'audit/Cohort_overview_two_generation_20261024')
status<-system2(cfg$python,c(shQuote(file.path(root,'scripts/audit_nature_cairo_pdf.py')),shQuote(file.path(cohort_out,'Study_population_and_analysis_sets.pdf')),
 '--skill',shQuote(nf_skill_root()),'--prefix',shQuote(file.path(cohort_audit,'Study_population_and_analysis_sets')),'--max-body-pt','9.01'))
if(status!=0L)stop('Cohort diagram rendered audit failed')
writeLines('STAGED_RENDERING_COMPLETED',file.path(root,'inventory/model_labels_20261030/STATUS'))
tf_log('All requested saved-result redraws completed; awaiting verification/publication.')
