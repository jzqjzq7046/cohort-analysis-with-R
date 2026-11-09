# Only redraw the five user-specified saved-result figures; no model is run.
root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(root,'R/literature_ABC_20261028.R'),encoding='UTF-8')
cfg<-lr_setup(root)
source(file.path(root,'R/genus_main_text_figures.R'),encoding='UTF-8')
audit<-file.path(root,'inventory/no_panel_titles_20261104')
targets<-list('24_Two_generation_reanalysis'=c('Figure2_Clinical_cross_generation',
 'Figure3_Family_microbiome_similarity','Figure4_Shannon_diversity'),
 '25_Taxa_robustness'='Figure25_Overall_abundance_ratios',
 '39_Genus_main_text_figures'='Figure39_Main_All_genus_volcano')
parents<-setNames(lapply(names(targets),function(m)lr_adopted(root,m)),names(targets))
protected<-unlist(lapply(names(targets),function(m){
 p<-parents[[m]]
 c(list.files(file.path(p,'private'),pattern='(stage_.*|analysis_result|display_result)[.]rds$',full.names=TRUE),
   list.files(file.path(p,'tables'),pattern='[.]csv$',full.names=TRUE))
}))
hashes<-unname(tools::md5sum(protected));stopifnot(!anyNA(hashes))
actual_tf_export<-tf_export;allowed<-character()
tf_export<-function(spec,stem,out,cfg){if(!stem%in%allowed)return(data.frame());actual_tf_export(spec,stem,out,cfg)}
for(module in names(targets)){
 from<-parents[[module]];allowed<-targets[[module]]
 out<-file.path(audit,'staging',module)
 for(d in c('figures','tables','reports','private','audit'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
 backup<-file.path(audit,'previous_figures',module);dir.create(backup,recursive=TRUE,showWarnings=FALSE)
 for(stem in allowed)for(ext in c('pdf','png','tiff','svg')){
   f<-file.path(from,'figures',paste0(stem,'.',ext))
   if(!file.exists(file.path(backup,basename(f))))stopifnot(file.copy(f,backup))
 }
 if(module=='24_Two_generation_reanalysis'){
   for(f in list.files(file.path(from,'private'),pattern='^clinical_pairs_.*[.]rds$',full.names=TRUE))
     stopifnot(file.copy(f,file.path(out,'private'),overwrite=TRUE))
   load_stage<-function(n)readRDS(file.path(from,'private',paste0('stage_',n,'.rds')))
   v3_visualize(load_stage('prepare'),load_stage('clinical'),load_stage('alpha'),load_stage('similarity'),NULL,NULL,cfg,out,only=allowed)
 }else if(module=='25_Taxa_robustness'){
   tf_figures25(readRDS(file.path(from,'private/analysis_result.rds')),cfg,out)
 }else{
   a<-readRDS(file.path(from,'private/display_result.rds'))
   lr_export(list(sgm_volcano(a,cfg)),allowed,out,cfg,width=183,height=112)
 }
}
stopifnot(identical(hashes,unname(tools::md5sum(protected))))
fl_write_csv(data.frame(file=protected,unchanged=TRUE,md5=hashes),file.path(audit,'unchanged_statistics.csv'))
jsonlite::write_json(lapply(names(targets),function(m)list(module=m,parent=parents[[m]],stems=targets[[m]])),
 file.path(audit,'targets.json'),auto_unbox=TRUE,pretty=TRUE)
cat('FIVE_FIGURES_RENDERED_FOR_REVIEW:',audit,'\n')
