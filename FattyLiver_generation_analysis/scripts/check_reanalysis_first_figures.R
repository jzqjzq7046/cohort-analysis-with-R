# Render and inspect only already completed clinical outputs during long model runs.
root<-'E:/YFS/workspace/FattyLiver_generation_analysis'
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),.libPaths()))
for(f in c('config.R','R/common.R','R/entry_helpers.R','R/nature_export_approved.R','R/reanalysis_numerical_helpers.R','R/reanalysis_20261024.R'))source(file.path(root,f),encoding='UTF-8')
cfg<-fl_config(root);Sys.setenv(PYTHONPATH=file.path(root,'.figure_qa_python'))
run<-readLines(file.path(root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),n=1)
out<-file.path(root,'outputs/24_Two_generation_reanalysis/_running',run)
x<-readRDS(file.path(out,'private/stage_prepare.rds'));cl<-readRDS(file.path(out,'private/stage_clinical.rds'));alpha<-readRDS(file.path(out,'private/stage_alpha.rds'))
similarity_file<-file.path(out,'private/stage_similarity.rds')
similarity<-if(file.exists(similarity_file))readRDS(similarity_file)else NULL
taxa_file<-file.path(out,'private/stage_taxa.rds');taxa<-if(file.exists(taxa_file))readRDS(taxa_file)else NULL
args<-commandArgs(TRUE)
only<-if('--taxa-only'%in%args)c('Figure6_B1_SGB_volcano','Figure6_B2_SGB_volcano','Figure7_B1_Representative_abundance','Figure7_B2_Representative_abundance')else NULL
v3_visualize(x,cl,alpha,similarity,NULL,taxa,cfg,out,only=only)
