# Read-only source inspection for the user-authorized A/D/E extension.
root <- 'E:/YFS/workspace/FattyLiver_generation_analysis'
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),file.path(R.home(),'library')))
pkgs <- c('SummarizedExperiment','TreeSummarizedExperiment','HDF5Array','rhdf5','DelayedArray','lme4','lmerTest','Maaslin2','jsonlite')
print(setNames(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE),pkgs))
x <- readRDS(file.path(root,'outputs/24_Two_generation_reanalysis/current/two_generation_20261024/private/stage_prepare.rds'))
cat('Prepared fields:',paste(names(x),collapse=', '),'\n')
cat('Metadata fields:',paste(names(x$metadata),collapse=', '),'\n')
cat('Relative matrix:',paste(dim(x$relative),collapse=' x '),'\n')
source(file.path(root,'R/reanalysis_numerical_helpers.R'))
source(file.path(root,'R/reanalysis_20261024.R'))
sc <- v3_micro_samples(x)
for(nm in names(sc)) {
 d<-v3_design(x,sc[[nm]]$mask,sc[[nm]]$fields)
 cat(nm,': N=',nrow(d$metadata),' eligible=',sum(d$filter$retained),' nuisance=',length(d$nuisance),'\n',sep='')
 print(table(d$metadata$generation,d$metadata$fatty_liver))
 cat('Bases: ',paste(d$bases,collapse=' '),'; nuisance: ',paste(d$nuisance,collapse=' '),'\n',sep='')
}
{
 for(kind in c('pathabundance','pathcoverage')) {
  p<-paste0('E:/YFS/documents-export-2026-08-31/tse_humann_',kind,'_h5se')
  z<-if(requireNamespace('HDF5Array',quietly=TRUE))HDF5Array::loadHDF5SummarizedExperiment(p)else readRDS(file.path(p,'se.rds'))
  cat('HUMAnN ',kind,': class=',class(z)[1],'; dimensions=',paste(dim(z),collapse=' x '),'\n',sep='')
  cat('Assays:',paste(SummarizedExperiment::assayNames(z),collapse=', '),'\n')
  cat('Row names:',paste(head(rownames(z),12),collapse='; '),'\n')
  cat('Column metadata:',paste(names(SummarizedExperiment::colData(z)),collapse='; '),'\n')
  cat('Row metadata:',paste(names(SummarizedExperiment::rowData(z)),collapse='; '),'\n')
  cat('Direct sample-name matches:',sum(x$metadata$sample_id%in%colnames(z)),'\n')
  cat('Unstratified rows:',sum(!grepl('|',rownames(z),fixed=TRUE)),'\n')
  cat('Metadata names:',paste(names(S4Vectors::metadata(z)),collapse=', '),'\n')
  cat('Row annotation preview:\n');print(head(as.data.frame(SummarizedExperiment::rowData(z)),3))
  if(requireNamespace('HDF5Array',quietly=TRUE))print(SummarizedExperiment::assay(z)[1:4,1:3])
 }
}
