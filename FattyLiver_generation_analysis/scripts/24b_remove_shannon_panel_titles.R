# Display-only revision: remove both panel titles and retain all saved inference.
root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
source(file.path(root,'R/literature_ABC_20261028.R'),encoding='UTF-8')
cfg<-lr_setup(root)
module<-file.path(root,'outputs/24_Two_generation_reanalysis')
id<-readLines(file.path(module,'LATEST_RUN.txt'),1,warn=FALSE)
from<-file.path(module,'current',id)
stopifnot(readLines(file.path(from,'run_status.txt'),1)=='COMPLETED')
stem<-'Figure4_Shannon_diversity'
audit<-file.path(root,'inventory/figure4_no_panel_titles_20261104')
out<-file.path(audit,'staging')
for(d in c('figures','tables','reports','private','audit'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
backup<-file.path(audit,'previous_figures');dir.create(backup,recursive=TRUE,showWarnings=FALSE)
files<-file.path(from,'figures',paste0(stem,c('.pdf','.png','.tiff','.svg')))
for(f in files)if(!file.exists(file.path(backup,basename(f))))stopifnot(file.copy(f,backup))
protected<-list.files(file.path(from,'private'),pattern='^stage_.*[.]rds$',full.names=TRUE)
hashes<-unname(tools::md5sum(protected))
for(f in list.files(file.path(from,'private'),pattern='^clinical_pairs_.*[.]rds$',full.names=TRUE))
  stopifnot(file.copy(f,file.path(out,'private'),overwrite=TRUE))
load_stage<-function(n)readRDS(file.path(from,'private',paste0('stage_',n,'.rds')))
v3_visualize(load_stage('prepare'),load_stage('clinical'),load_stage('alpha'),
  load_stage('similarity'),NULL,NULL,cfg,out,only=stem)
stopifnot(identical(hashes,unname(tools::md5sum(protected))))
fl_write_csv(data.frame(file=basename(protected),unchanged=TRUE,md5=hashes),
  file.path(out,'audit/unchanged_saved_analysis.csv'))
writeLines(c('Display-only revision: remove a/b panel titles. Panel tags, data, inference, palette and axes retained.',
  'No model refit or P/q recalculation.'),file.path(out,'audit/display_revision.txt'))
cat('RENDERED_FOR_REVIEW:',out,'\n')
