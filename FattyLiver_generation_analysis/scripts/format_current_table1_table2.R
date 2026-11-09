# Descriptive summaries and presentation only; existing regression fits are read.
root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),.libPaths()))
for(f in c('config.R','R/common.R','R/nature_export_approved.R','R/reanalysis_numerical_helpers.R','R/reanalysis_20261024.R'))source(file.path(root,f),encoding='UTF-8')
source(file.path(root,'R/table2_presentation.R'),encoding='UTF-8')
.table2_only<-'--table2-only'%in%commandArgs(TRUE)
cfg<-fl_config(root);Sys.setenv(PYTHONPATH=file.path(root,'.figure_qa_python'))
run<-readLines(file.path(root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),n=1)
current<-file.path(root,'outputs/24_Two_generation_reanalysis/current',run)
out<-file.path(current,'tables/Manuscript_Table1_Table2_20261024');dir.create(out,recursive=TRUE,showWarnings=FALSE)
audit<-file.path(current,'audit','Table1_Table2_20261024');dir.create(audit,recursive=TRUE,showWarnings=FALSE)
models_path<-file.path(current,'private/stage_clinical.rds');protected<-c(file.path(current,'private/stage_prepare.rds'),models_path);before<-unname(tools::md5sum(protected))
x<-readRDS(protected[1]);m<-x$clinical
stopifnot(nrow(m)==5660L,!anyDuplicated(m$participant_id))
r<-read.csv(file.path(current,'private/clinical_selected.csv'),colClasses='character',fileEncoding='UTF-8-BOM')
j<-match(m$participant_id,r$tutknom18);stopifnot(!anyNA(j));m$smoking<-suppressWarnings(as.numeric(r$smoke18[j]));stopifnot(all(is.na(m$smoking)|m$smoking%in%0:1))
smoke_check<-match(x$metadata$participant_id,m$participant_id);stopifnot(!anyNA(smoke_check),all(m$smoking[smoke_check]==x$metadata$smoking,na.rm=TRUE))
m$abx_recent<-x$metadata$abx_recent[match(m$participant_id,x$metadata$participant_id)]
labels<-read.csv(file.path(out,'Table1_Cohort.csv'),check.names=FALSE)[[1]]
continuous<-c('age','BMI18','WAIST18','syst18','DKVAI18','GLUK18','HbA1c18','CRP18','TOTKOL18','TRIGLY18','HDLKOL18','LDLKOL18','ALAT18','ASAT18','alkodpd18','dietscore_g18')
row_for<-c(2,5,6,10:22)
stopifnot(length(row_for)==length(continuous),length(labels)==24L)
t1<-data.frame(Characteristic=labels,check.names=FALSE);detail<-list()
for(g in c('Overall','G0','G1','G2')){
 z<-if(g=='Overall')m else m[m$generation==g,,drop=FALSE];v<-vector('list',24);v[[1]]<-nrow(z)
 for(i in seq_along(continuous)){
  field<-continuous[i];a<-z[[field]];a<-a[is.finite(a)];digits<-if(field%in%c('CRP18','alkodpd18'))2L else 1L;q<-quantile(a,c(.25,.5,.75),names=FALSE)
  v[[row_for[i]]]<-sprintf(paste0('%.',digits,'f [%.',digits,'f, %.',digits,'f]; n=%d'),q[2],q[1],q[3],length(a))
  detail[[paste(g,field)]]<-data.frame(group=g,variable=field,N=length(a),median=q[2],Q1=q[1],Q3=q[3])
 }
 binary<-list(`3`=z$sex=='Female',`4`=z$age<18,`7`=z$fatty_liver==1,`8`=z$fat_code==2,`9`=z$fat_code==3,`23`=z$smoking==1,`24`=z$abx_recent==1)
 for(i in names(binary)){a<-binary[[i]];den<-sum(!is.na(a));v[[as.integer(i)]]<-sprintf('%d/%d (%.1f%%)',sum(a,na.rm=TRUE),den,100*mean(a,na.rm=TRUE))}
 t1[[paste0(g,' (N = ',format(nrow(z),big.mark=','),')')]]<-unlist(v)
}
cl<-readRDS(models_path)$models
prev<-read.csv(file.path(current,'tables/Figure2_prevalence_source_EN.csv'),check.names=FALSE)
raw<-cl[cl$model%in%c('Unadjusted','C1'),];stopifnot(nrow(raw)==6L,all(grepl('^ESTIMATED',raw$status)))
cells<-function(comparison,exposure){a<-prev[prev$comparison==comparison&prev$exposure==exposure,];stopifnot(nrow(a)==1L);sprintf('%d/%d (%.1f%%)',a$cases,a$N,a$percent)}
t2<-fl_table2_data(raw,prev)
stopifnot(identical(t1$Characteristic,labels),sum(m$fatty_liver)==1054L,all(raw$n%in%c(1097,2262,1329)))
notes1<-c('Values are median [25th, 75th percentiles] with available n, or n/N (%). N is the number with non-missing data for each characteristic. No regression adjustment is applied.',
 'The clinical cohort includes 5,660 participants in 1,693 primary families with at least two observed eligible generations. All-age analyses retain minors. A stool sample is not required for clinical inclusion.',
 'Hepatic steatosis is defined by ultrasound (mild or clear fat accumulation); MASLD is not classified. G0, G1 and G2 denote the recorded generations.',
 'Recent antibiotic use is available only in the nested microbiome cohort (3,217 participants). Its denominator refers to available stool-questionnaire data in that subset, not all 5,660 clinical participants.',
 'BMI, body mass index; CRP, C-reactive protein; HDL, high-density lipoprotein; LDL, low-density lipoprotein; ALT, alanine aminotransferase; AST, aspartate aminotransferase.')
notes2<-fl_table2_notes()
if(!.table2_only){fl_write_csv(t1,file.path(out,'Table1_Cohort.csv'));fl_write_csv(ov_bind(detail),file.path(audit,'Table1_continuous_summary.csv'))}
fl_write_csv(t2,file.path(out,'Table2_Clinical_association.csv'))
spec<-list(Table1=list(title='Table 1. Clinical characteristics and hepatic steatosis by generation',columns=names(t1),rows=unname(lapply(seq_len(nrow(t1)),function(i){z<-as.list(t1[i,]);if(i==1)z[-1]<-lapply(z[-1],as.numeric);unname(z)})),notes=notes1),Table2=list(title='Table 2. Cross-generation association of hepatic steatosis',columns=names(t2),rows=unname(lapply(seq_len(nrow(t2)),function(i)unname(as.list(t2[i,])))),notes=notes2))
spec$Table2$group_header<-fl_table2_group;spec$Table2$group_columns<-c(4L,5L)
jsonlite::write_json(spec,file.path(audit,'table_spec.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
wrap_mm<-function(text,width,font){
 if(!nzchar(text))return('');lines<-character();for(part in strsplit(text,'\n',fixed=TRUE)[[1]]){words<-strsplit(part,' ',fixed=TRUE)[[1]];line<-'';for(word in words){test<-trimws(paste(line,word));if(nzchar(line)&&grid::convertWidth(grid::stringWidth(test),'mm',valueOnly=TRUE)>width){lines<-c(lines,line);line<-word}else line<-test};lines<-c(lines,line)};paste(lines,collapse='\n')
}
draw_table<-function(d,title,notes,width,height,font,cols,headerheight){
 if(ncol(d)==8L)return(fl_draw_table2(d,title,notes,width,height,font,cols,headerheight))
 grid::grid.newpage();grid::pushViewport(grid::viewport(gp=grid::gpar(fontfamily='Times New Roman',fontsize=font)))
 margin<-8;available<-width-2*margin;cw<-available*cols;edges<-c(margin,margin+cumsum(cw));top<-18
 tx<-function(text,x,y,just='centre',size=font,face='plain')grid::grid.text(text,x=grid::unit(x,'mm'),y=grid::unit(height-y,'mm'),just=just,gp=grid::gpar(fontsize=size,fontface=face),default.units='mm')
 rule<-function(y,lwd=0.7)grid::grid.lines(x=grid::unit(c(margin,width-margin),'mm'),y=grid::unit(rep(height-y,2),'mm'),gp=grid::gpar(lwd=lwd))
 tx(title,width/2,8,size=font+2,face='bold');rule(top,.9)
 for(j in seq_len(ncol(d))){text<-wrap_mm(names(d)[j],cw[j]-2.2,font);tx(text,(edges[j]+edges[j+1])/2,top+headerheight/2,size=font,face='bold')}
 y<-top+headerheight;rule(y,.6)
 for(i in seq_len(nrow(d))){
  values<-vapply(seq_len(ncol(d)),function(j){v<-d[[j]][i];if(is.numeric(v)){if(j==3&&ncol(d)==9)format(v,big.mark=',',scientific=FALSE)else if(j==9&&ncol(d)==9){if(v<.0001)format(v,scientific=TRUE,digits=3)else sprintf('%.4f',v)}else sprintf('%.2f',v)}else as.character(v)},character(1))
  if(ncol(d)==5)values[-1]<-gsub('; n=','\nn=',values[-1],fixed=TRUE)
  texts<-vapply(seq_along(values),function(j)wrap_mm(values[j],cw[j]-2.4,font),character(1));lines<-max(lengths(strsplit(texts,'\n',fixed=TRUE)));rh<-max(5.3,lines*font*25.4/72*1.1+1.7)
  for(j in seq_along(values))tx(texts[j],if(j==1)edges[j]+1 else (edges[j]+edges[j+1])/2,y+rh/2,just=if(j==1)'left'else'centre')
  y<-y+rh
 };rule(y,.9);y<-y+5
 grid::pushViewport(grid::viewport(gp=grid::gpar(fontsize=font-.5)))
 for(note in notes){text<-wrap_mm(note,available,font-.5);n<-length(strsplit(text,'\n',fixed=TRUE)[[1]]);tx(text,margin,y,just=c('left','top'),size=font-.5);y<-y+n*(font-.5)*25.4/72*1.15+2}
 grid::popViewport(2);stopifnot(y<height-3)
}
export_table<-function(d,spec,stem,width,height,font,cols,headerheight){
 draw<-function()draw_table(d,spec$title,spec$notes,width,height,font,cols,headerheight)
 for(ext in c('pdf','png','svg','tiff')){
  path<-file.path(out,paste0(stem,'.',ext))
  switch(ext,pdf=grDevices::cairo_pdf(path,width=width/25.4,height=height/25.4,family='Times New Roman'),png=grDevices::png(path,width=width,height=height,units='mm',res=600,type='cairo',bg='white'),svg=svglite::svglite(path,width=width/25.4,height=height/25.4,bg='white'),tiff=grDevices::tiff(path,width=width,height=height,units='mm',res=600,type='cairo',compression='lzw',bg='white'))
  tryCatch(draw(),finally=grDevices::dev.off())
 }
 raw_text<-suppressWarnings(system2(cfg$python,c(shQuote(file.path(nf_skill_root(),'scripts/audit_pdf_text.py')),shQuote(file.path(out,paste0(stem,'.pdf'))),'--min-pt','5','--json'),stdout=TRUE,stderr=TRUE));writeLines(raw_text,file.path(audit,paste0(stem,'_raw_Tf.json')))
 status<-system2(cfg$python,c(shQuote(file.path(root,'scripts/audit_nature_cairo_pdf.py')),shQuote(file.path(out,paste0(stem,'.pdf'))),'--skill',shQuote(nf_skill_root()),'--prefix',shQuote(file.path(audit,stem)),'--max-body-pt','11.01'))
 if(status!=0)stop('Table rendered QA failed: ',stem)
}
if(!.table2_only)export_table(t1,spec$Table1,'Table1_Cohort',210,270,8,c(.34,rep(.165,4)),12)
export_table(t2,spec$Table2,'Table2_Clinical_association',267,140,9,c(.10,.07,.12,.20,.20,.06,.15,.10),25)
stopifnot(identical(before,unname(tools::md5sum(protected))))
writeLines('PASS: original cohort RDS and clinical model RDS unchanged; clinical N=5660; clinical Table2 values match completed C1 models.',file.path(audit,'validation.txt'))
writeLines(out,file.path(root,'inventory/current_table1_table2_output.txt'))
cat('TABLES PREPARED:',out,'\n')
