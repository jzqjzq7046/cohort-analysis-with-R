root <- 'E:/YFS/workspace/FattyLiver_generation_analysis'
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),.libPaths()))
for(f in c('config.R','R/common.R','R/nature_export_approved.R','R/reanalysis_numerical_helpers.R','R/reanalysis_20261024.R','R/reanalysis_clinical_mi.R'))source(file.path(root,f),encoding='UTF-8')
cfg<-fl_config(root);Sys.setenv(PYTHONPATH=file.path(root,'.figure_qa_python'))
run<-readLines(file.path(root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),n=1)
out<-file.path(root,'outputs/24_Two_generation_reanalysis/_running',run)
clinical<-readRDS(file.path(out,'private/stage_clinical.rds'));mi<-readRDS(file.path(out,'private/stage_clinical_mi.rds'))
print(v3_mi_visualize(clinical,mi,cfg,out))
