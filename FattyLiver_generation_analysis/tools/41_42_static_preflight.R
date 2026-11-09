# Validate the actual R plotting source together with its delegated exporter.
root <- Sys.getenv('FATTYLIVER_ROOT','E:/YFS/workspace/FattyLiver_generation_analysis')
source(file.path(root,'R/generation_ADE_extensions.R'),encoding='UTF-8')
cfg<-ade_setup(root)
target<-file.path(root,'inventory/ADE_20261101/plotting_source_with_exporter.R')
f<-function(name)paste0(name,' <- ',paste(deparse(get(name)),collapse='\n'))
writeLines(c(readLines(file.path(root,'R/generation_ADE_figures.R'),warn=FALSE),
  readLines(file.path(root,'R/generation_A_top10_figure.R'),warn=FALSE),f('tf_export'),
  f('ov_export'),f('tf_plot_theme'),f('nf_theme')),target,useBytes=TRUE)
script<-file.path(nf_skill_root(),'scripts/validate_figure.py')
res<-system2(cfg$python,c(shQuote(script),shQuote(target),'--backend','r','--json'),stdout=TRUE)
writeLines(res,file.path(root,'inventory/ADE_20261101/static_figure_preflight.json'))
z<-jsonlite::fromJSON(paste(res,collapse='\n'))
cat(paste(res,collapse='\n'),'\n')
stopifnot(z$summary$counts$FAIL==0)
