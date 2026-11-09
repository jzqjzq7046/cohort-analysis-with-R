# Visualization-only revision; never calls a statistical analysis stage.
.revision_root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(.revision_root,'.r-library-sgb'),file.path(.revision_root,'.r-library'),.libPaths()))
for(f in c('config.R','R/common.R','R/entry_helpers.R','R/nature_export_approved.R','R/reanalysis_numerical_helpers.R','R/reanalysis_20261024.R','R/reanalysis_clinical_mi.R','R/result_gallery.R'))source(file.path(.revision_root,f),encoding='UTF-8')
revise_publish_gallery<-function(current,cfg,backup){
 root<-cfg$project_root;central<-file.path(root,'Results_表格与图片');manifest<-file.path(central,'file_manifest.csv')
 if(file.exists(manifest)){
  old<-read.csv(manifest,fileEncoding='UTF-8-BOM',check.names=FALSE)
  retire<-old$module=='24_Two_generation_reanalysis'&grepl('^Figure6_(B1|B2)_SGB_volcano[.](pdf|png|tiff|svg)$',old$file)
  if(any(retire)){
   dir.create(file.path(backup,'central_gallery'),recursive=TRUE,showWarnings=FALSE)
   stopifnot(file.copy(manifest,file.path(backup,'central_gallery/file_manifest.csv'),overwrite=TRUE))
   base<-paste0(tolower(normalizePath(central,winslash='/',mustWork=TRUE)),'/')
   for(p in old$viewing_path[retire]){
    path<-file.path(central,p)
    if(file.exists(path)){
     stopifnot(startsWith(tolower(normalizePath(path,winslash='/',mustWork=TRUE)),base))
     dest<-file.path(backup,'central_gallery',p);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
     stopifnot(file.copy(path,dest,overwrite=TRUE),identical(unname(tools::md5sum(path)),unname(tools::md5sum(dest))))
     unlink(path)
     if(file.exists(path))message('Archived previous viewing copy remains open; retain the file but exclude it from the current index: ',basename(path))
    }
   }
   fl_write_csv(old[!retire,,drop=FALSE],manifest)
  }
 }
 tryCatch(fl_refresh_result_gallery(root),error=function(e){
  saved<-file.path(backup,'central_gallery/file_manifest.csv')
  if(file.exists(saved)){
   original<-read.csv(saved,fileEncoding='UTF-8-BOM',check.names=FALSE)
   fixed<-original[original$module=='24_Two_generation_reanalysis'&grepl('^Figure6_(B1|B2)_SGB_volcano[.](pdf|png|tiff|svg)$',original$file),,drop=FALSE]
   base<-paste0(tolower(normalizePath(central,winslash='/',mustWork=TRUE)),'/')
   for(p in fixed$viewing_path){
    src<-file.path(backup,'central_gallery',p);dest<-file.path(central,p)
    stopifnot(startsWith(paste0(tolower(normalizePath(dirname(dest),winslash='/',mustWork=TRUE)),'/'),base))
    if(file.exists(src)&&(!file.exists(dest)||!identical(unname(tools::md5sum(src)),unname(tools::md5sum(dest)))))stopifnot(file.copy(src,dest,overwrite=TRUE))
   }
   stopifnot(file.copy(saved,manifest,overwrite=TRUE))
  }
  message('Central viewer PDFs are open; previous central viewing copies are retained. The dedicated revised-result viewer is published independently.')
 })
 v3_publish_view(current,cfg)
}
revise_two_generation_figures<-function(root=.revision_root){
 cfg<-fl_config(root);Sys.setenv(PYTHONPATH=file.path(root,'.figure_qa_python'))
 run<-readLines(file.path(root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),n=1)
 current<-file.path(root,'outputs/24_Two_generation_reanalysis/current',run)
 stopifnot(readLines(file.path(current,'run_status.txt'),n=1)=='COMPLETED')
 revision<-file.path(root,'outputs/24_Two_generation_reanalysis/_figure_revisions',format(Sys.time(),'%Y%m%d_%H%M%S'))
 stage<-file.path(revision,'staging');backup<-file.path(revision,'previous')
 for(p in c('figures','tables','reports','audit','private'))dir.create(file.path(stage,p),recursive=TRUE,showWarnings=FALSE)
 dir.create(backup,recursive=TRUE,showWarnings=FALSE)
 protected<-c(list.files(file.path(current,'private'),pattern='^stage_.*[.]rds$',full.names=TRUE),setdiff(list.files(file.path(current,'tables'),pattern='[.]csv$',full.names=TRUE),file.path(current,'tables/Figure_export_manifest_EN.csv')))
 before<-unname(tools::md5sum(protected))
 for(folder in c('tables','reports'))file.copy(list.files(file.path(current,folder),full.names=TRUE),file.path(stage,folder),overwrite=TRUE)
 file.copy(list.files(file.path(current,'private'),pattern='^clinical_pairs_.*[.]rds$',full.names=TRUE),file.path(stage,'private'))
 loadstage<-function(n)readRDS(file.path(current,'private',paste0('stage_',n,'.rds')))
 selected<-c('Figure2_Clinical_cross_generation','Figure3_Family_microbiome_similarity','Figure6_Overall_and_adult_SGB_volcano')
 v3_visualize(loadstage('prepare'),loadstage('clinical'),loadstage('alpha'),loadstage('similarity'),loadstage('composition'),loadstage('taxa'),cfg,stage,only=selected)
 # Keep the existing workbook; update only the figure-manifest worksheet.
 workbook<-file.path(stage,'tables/Hepatic_steatosis_two_generation_tables_EN.xlsx')
 wb<-openxlsx::loadWorkbook(workbook);contents<-openxlsx::read.xlsx(wb,'Contents')
 sheet<-contents$sheet[contents$source=='Figure_export_manifest_EN.csv'];stopifnot(length(sheet)==1L)
 manifest<-read.csv(file.path(stage,'tables/Figure_export_manifest_EN.csv'))
 openxlsx::deleteData(wb,sheet,cols=1:6,rows=1:40,gridExpand=TRUE)
 openxlsx::writeData(wb,sheet,manifest,headerStyle=openxlsx::createStyle(fontName='Times New Roman',fontSize=11,textDecoration='bold',border=c('top','bottom'),wrapText=TRUE))
 openxlsx::saveWorkbook(wb,workbook,overwrite=TRUE)
 for(f in list.files(file.path(stage,'reports'),full.names=TRUE)){
  lines<-readLines(f,warn=FALSE,encoding='UTF-8')
  lines<-gsub('13张图','12张图',lines,fixed=TRUE)
  lines<-gsub('Panel a shows mean Bray-Curtis distances','Panel a connects the IF and BF group-mean Bray-Curtis distances',lines,fixed=TRUE)
  lines<-gsub('Figure 6 (B1 and B2).','Figure 6. Panel a: all-age B1; panel b: adult B2. Panels use separate axis ranges.',lines,fixed=TRUE)
  lines<-gsub('Similarity bars are means or adjusted coefficients','Similarity points are group means and bars are adjusted coefficients; connecting lines join group means, not individual observations',lines,fixed=TRUE)
  writeLines(enc2utf8(lines),f,useBytes=TRUE)
 }
 writeLines(c('# Figure revision','', 'Figure 2: filled circle/diamond symbols have black outlines.','Figure 3: connected IF/BF group means with numeric labels, and labelled horizontal adjusted-difference bars.','Figure 6: all-age B1 (left) and adult B2 (right) are assembled into one two-panel figure.','Statistical inputs, estimates, confidence intervals and P values are unchanged. Previous figure versions are archived outside the current viewing folder.'),file.path(stage,'reports/Figure_revision_EN.md'))
 helpers<-c('v3_style','v3_visualize','v3_mi_visualize','nf_theme','nf_skill_root','nf_alignment','nf_export','nf_qa_command','ov_export')
 plot_source<-file.path(stage,'audit/plot_and_export_source.R');writeLines(unlist(lapply(helpers,function(n)c(paste0(n,' <- '),deparse(get(n),width.cutoff=120L),''))),plot_source)
 invisible(parse(plot_source));nf_qa_command(cfg,'validate_figure.py',c(shQuote(plot_source),'--json'),file.path(stage,'audit/plot_source_preflight.json'))
 stopifnot(identical(before,unname(tools::md5sum(protected))))
 # All requested exports passed rendered QA before publication starts.
 archive_file<-function(path,relative){if(file.exists(path)){dest<-file.path(backup,relative);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(path,dest,overwrite=TRUE),identical(unname(tools::md5sum(path)),unname(tools::md5sum(dest))))}}
 obsolete<-c('Figure6_B1_SGB_volcano','Figure6_B2_SGB_volcano')
 for(stem in c(selected,obsolete))for(ext in c('pdf','png','tiff','svg'))archive_file(file.path(current,'figures',paste0(stem,'.',ext)),file.path('figures',paste0(stem,'.',ext)))
 for(folder in c('tables','reports'))for(f in list.files(file.path(stage,folder),full.names=TRUE)){
  relative<-file.path(folder,basename(f));dest<-file.path(current,relative);archive_file(dest,relative)
  stopifnot(file.copy(f,dest,overwrite=TRUE))
 }
 for(f in list.files(file.path(stage,'figures'),full.names=TRUE))stopifnot(file.copy(f,file.path(current,'figures',basename(f)),overwrite=TRUE))
 for(stem in obsolete)for(ext in c('pdf','png','tiff','svg')){
  unlink(file.path(current,'figures',paste0(stem,'.',ext)))
  viewed<-file.path(root,'Results_表格与图片/Reanalysis_two_generation_20261024',run,'figures',paste0(stem,'.',ext))
  archive_file(viewed,file.path('previous_viewer',paste0(stem,'.',ext)));unlink(viewed)
 }
 for(f in list.files(file.path(stage,'audit'),full.names=TRUE)){dest<-file.path(current,'audit',basename(f));archive_file(dest,file.path('audit',basename(f)));file.copy(f,dest,overwrite=TRUE)}
 archive_file(file.path(current,'audit/reanalysis_20261024.R'),'audit/reanalysis_20261024.R')
 file.copy(file.path(root,'R/reanalysis_20261024.R'),file.path(current,'audit/reanalysis_20261024.R'),overwrite=TRUE)
 stopifnot(identical(before,unname(tools::md5sum(protected))))
 writeLines(c('PASS: visualization-only revision; protected stage RDS and statistical CSV MD5 unchanged.',paste0('Previous versions: ',revision)),file.path(current,'audit/figure_revision_validation.txt'))
 view<-revise_publish_gallery(current,cfg,backup)
 cat('COMPLETED VISUALIZATION REVISION:',view,'\n')
 invisible(view)
}
if(sys.nframe()==0L)revise_two_generation_figures()else if(isTRUE(getOption('fatty_liver.autorun',TRUE)))revise_two_generation_figures()
