# Audit the actual plotting/export helpers, avoiding statistics-only source text.
root<-'E:/YFS/workspace/FattyLiver_generation_analysis'
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),.libPaths()))
for(f in c('config.R','R/common.R','R/nature_export_approved.R','R/reanalysis_numerical_helpers.R','R/reanalysis_20261024.R','R/reanalysis_clinical_mi.R'))source(file.path(root,f),encoding='UTF-8')
run<-readLines(file.path(root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),n=1)
audit<-file.path(root,'outputs/24_Two_generation_reanalysis/_running',run,'audit')
names<-c('v3_style','v3_visualize','v3_mi_visualize','nf_theme','nf_skill_root','nf_alignment','nf_export','nf_qa_command','ov_export')
lines<-unlist(lapply(names,function(n)c(paste0(n,' <- '),deparse(get(n),width.cutoff=120L),'')))
path<-file.path(audit,'plot_and_export_source.R');writeLines(lines,path)
invisible(parse(path));cfg<-fl_config(root)
nf_qa_command(cfg,'validate_figure.py',c(shQuote(path),'--json'),file.path(audit,'plot_source_preflight.json'))
