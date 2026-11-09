root <- 'E:/YFS/workspace/FattyLiver_generation_analysis'
.libPaths(c(file.path(root,'.r-library-sgb'),.libPaths()))
source(file.path(root,'R/confounders.R'))
for (v in c('processed_microbiome_v1','processed_microbiome_three_generation_v1_20261019_111036')) {
 d <- readRDS(file.path(root,'data/processed',v,'microbiome_processed.rds'))
 cat('\nVERSION',v,'N',nrow(d$metadata),'FIELDS',paste(names(d),collapse=','),'\n')
 cat('metadata:',paste(names(d$metadata),collapse=','),'\n')
 cat('counts',dim(d$full$counts),'filtered',dim(d$matrices$counts),'\n')
 cat('features:',paste(names(d$features),collapse=','),'\n')
 print(head(d$features,2)); print(d$source_paths)
 rm(d);gc(FALSE)
}
cd <- fl_cv_coldata('E:/YFS/documents-export-2026-08-31/tse_metaphlan_2025.rds')
cat('raw metadata',nrow(cd),'rownames',head(rownames(cd)),'\n')
for(v in c('sample_id','tutknom18','tutkno80perhe','sp','smoke18','ANTBV18','PAI18','RASVAMAKSA18','sukupolvi','prep_batch','seq_pool')) if(v %in% names(cd)){cat(v,':');print(head(cd[[v]]));print(table(cd[[v]],useNA='ifany'))}
cat('packages\n');print(sapply(c('Maaslin2','lmerTest','geepack','vegan','detectseparation','openxlsx','ggplot2','patchwork','svglite','digest'),function(p) if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else 'MISSING'))
