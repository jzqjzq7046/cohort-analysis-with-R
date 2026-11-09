# Display-only repair requested 2026-11-05. No model fitting; never reads/writes PPTX.
root <- Sys.getenv('FATTYLIVER_ROOT', unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', '.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),file.path(R.home(),'library')))
Sys.setenv(FATTYLIVER_ROOT=root,PYTHONPATH=file.path(root,'.figure_qa_python'))
for (f in c('config.R','R/common.R','R/entry_helpers.R','R/nature_export_approved.R',
            'R/reanalysis_numerical_helpers.R','R/reanalysis_20261024.R','R/result_gallery.R'))
  source(file.path(root,f),encoding='UTF-8')
cfg <- fl_config(root)
run <- readLines(file.path(root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),n=1,warn=FALSE)
current <- file.path(root,'outputs/24_Two_generation_reanalysis/current',run)
stopifnot(readLines(file.path(current,'run_status.txt'),n=1,warn=FALSE)=='COMPLETED')
stem <- 'Figure4_Shannon_diversity'
revision <- file.path(root,'outputs/24_Two_generation_reanalysis/_figure_revisions',
                      paste0(format(Sys.time(),'%Y%m%d_%H%M%S'),'_Shannon_legend'))
stage <- file.path(revision,'staging')
for (d in c('figures','audit','previous/figures','previous/audit'))
  dir.create(file.path(stage,d),recursive=TRUE,showWarnings=FALSE)
protected <- c(list.files(file.path(current,'private'),full.names=TRUE,recursive=TRUE),
               list.files(file.path(current,'tables'),full.names=TRUE,recursive=TRUE),
               list.files(file.path(current,'reports'),full.names=TRUE,recursive=TRUE),
               list.files(file.path(current,'figures'),full.names=TRUE))
protected <- protected[!grepl(paste0('[/\\\\]',stem,'[.]'),protected)]
fingerprints <- data.frame(path=protected,md5=unname(tools::md5sum(protected)))
check_protected <- function() stopifnot(identical(unname(tools::md5sum(fingerprints$path)),fingerprints$md5))
write.csv(fingerprints,file.path(stage,'audit/protected_fingerprints.csv'),row.names=FALSE)

# Execute just the existing five-line Shannon plotting block, avoiding
# every analysis stage and every unrelated plotting/table-writing block.
x <- readRDS(file.path(current,'private/stage_prepare.rds'))
alpha <- readRDS(file.path(current,'private/stage_alpha.rds'))
binary <- setNames(unname(cfg$binary_colors),c('No steatosis','Steatosis'))
patch <- function(plots,ncol=2) patchwork::wrap_plots(lapply(seq_along(plots),function(i)
  plots[[i]]+v3_style()+ggplot2::labs(tag=letters[i])),ncol=ncol)
captured <- NULL
export <- function(p,name,width=183,height=115,multi=FALSE,panels=1L) {
  stopifnot(identical(name,stem),multi,panels==2L)
  captured <<- list(name=name,plot=p,width=width,height=height,panels=panels,export_tiff=TRUE)
}
lines <- readLines(file.path(root,'R/reanalysis_20261024.R'),encoding='UTF-8',warn=FALSE)
finish <- grep("export(patch(list(p1,p2)),'Figure4_Shannon_diversity'",lines,fixed=TRUE)
stopifnot(length(finish)==1L)
block <- lines[(finish-4L):finish]
stopifnot(grepl('m<-x$metadata[v3_micro_samples(x)$B1$mask,]',block[1],fixed=TRUE),
          length(parse(text=block))==7L)
eval(parse(text=block),envir=.GlobalEnv)
stopifnot(!is.null(captured),nrow(m)==3193L,nrow(a)==3L)
built_panel_b <- ggplot2::ggplot_build(p2)
stopifnot(identical(as.character(built_panel_b$plot$scales$get_scales('colour')$get_breaks()),c('G0','G1','G2')))
saved <- alpha$effects[alpha$effects$model=='B1' & alpha$effects$generation!='Overall',]
stopifnot(identical(a$beta,saved$beta),identical(a$CI_low,saved$CI_low),identical(a$CI_high,saved$CI_high))
plot_source <- file.path(stage,'audit/Shannon_legend_plot_source.R')
helpers <- c('v3_style','nf_theme','patch','nf_alignment','nf_export','nf_qa_command','nf_skill_root')
writeLines(c('# Exact source used for the requested legend-only revision.',
             unlist(lapply(helpers,function(n)c(paste0(n,' <- '),deparse(get(n)),''))),block),plot_source)
invisible(parse(plot_source))
nf_qa_command(cfg,'validate_figure.py',c(shQuote(plot_source),'--json'),
              file.path(stage,'audit/Shannon_legend_source_preflight.json'))
nf_export(captured,stage,cfg)
check_protected()

# Archive the four old formats and matching QA before updating the current image.
for (ext in c('pdf','png','svg','tiff')) {
  old <- file.path(current,'figures',paste0(stem,'.',ext))
  backup <- file.path(stage,'previous/figures',basename(old))
  stopifnot(file.exists(old),file.copy(old,backup),
            identical(unname(tools::md5sum(old)),unname(tools::md5sum(backup))))
}
qa <- list.files(file.path(stage,'audit'),full.names=TRUE)
for (f in qa) {
  old <- file.path(current,'audit',basename(f))
  if (file.exists(old)) stopifnot(file.copy(old,file.path(stage,'previous/audit',basename(f))))
}
for (f in list.files(file.path(stage,'figures'),full.names=TRUE))
  stopifnot(file.copy(f,file.path(current,'figures',basename(f)),overwrite=TRUE))
for (f in qa) stopifnot(file.copy(f,file.path(current,'audit',basename(f)),overwrite=TRUE))
check_protected()
fl_refresh_result_gallery(root)
check_protected()
writeLines(c('PASS: generation legend restored in panel b; order G0, G1, G2.',
             'PASS: original source sample, estimates, confidence intervals and all protected files unchanged.',
             'PASS: unchanged strict rendered alignment, effective-font and collision gates.',
             'No statistical models refitted. No PPTX file accessed or modified.',
             paste0('Previous images and QA: ',file.path(stage,'previous'))),
           file.path(current,'audit/Shannon_legend_revision_20261105.txt'))
cat('COMPLETED:',file.path(current,'figures',paste0(stem,'.png')),'\n')
cat('REVISION:',revision,'\n')
