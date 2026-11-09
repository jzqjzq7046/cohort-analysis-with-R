# Source-backed two-branch cohort overview; drawing/export uses R only.
# No statistical models are fitted and no original analysis files are modified.
root <- Sys.getenv('FATTYLIVER_ROOT', unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if (.Platform$OS.type=='windows') invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),.libPaths()))
source(file.path(root,'R/model_display_labels.R'),encoding='UTF-8')
library(grid)
run <- readLines(file.path(root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),n=1)
current <- file.path(root,'outputs/24_Two_generation_reanalysis/current',run)
folder <- 'Cohort_overview_two_generation_20261024'
out <- file.path(current,'figures',folder)
audit <- file.path(current,'audit',folder)
# Preserve the adopted diagram and its audit before this appearance-only revision.
if (file.exists(file.path(out,'Study_population_and_analysis_sets.pdf'))) {
  backup <- file.path(root,'outputs/24_Two_generation_reanalysis/history/figure_revisions',
    paste0('cohort_overview_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid()))
  dir.create(backup,recursive=TRUE,showWarnings=FALSE)
  for (p in c(out,audit)) if (dir.exists(p)) {
    target <- file.path(backup,if(identical(p,out))'figures' else 'audit')
    dir.create(target,recursive=TRUE,showWarnings=FALSE)
    fs <- list.files(p,full.names=TRUE,all.files=TRUE,no..=TRUE)
    if(length(fs) && !all(file.copy(fs,target,recursive=TRUE,copy.date=TRUE)))
      stop('Cannot preserve the preceding cohort diagram')
  }
}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(audit,recursive=TRUE,showWarnings=FALSE)
read_tab <- function(name) read.csv(file.path(current,'tables',paste0(name,'.csv')),check.names=FALSE)
co <- read_tab('00_Cohort_counts_EN'); flow <- read_tab('00_Sample_flow_EN')
cl <- read_tab('Clinical_pair_sample_flow_EN'); mi <- read_tab('Adult_clinical_MI_sample_flow_EN')
pa <- read_tab('Pair_sample_flow_EN'); sg <- read_tab('SGB_model_summary_EN')
protected <- list.files(file.path(current,'tables'),pattern='\\.csv$',full.names=TRUE)
before <- tools::md5sum(protected)
get_n <- function(stage) { z<-flow$N[flow$stage==stage];stopifnot(length(z)==1L);z }
get_co <- function(kind,group='Overall',field='N') {
  z<-co[co$cohort==kind & co$group==group,field];stopifnot(length(z)==1L);z
}
clinical <- 'Clinical: >=2 observed generations'; microbial <- 'Microbiome: >=2 observed generations'
comparisons <- c('G0-G1','G1-G2','G0-G2')
clinical_N <- vapply(comparisons,function(k){z<-cl$n[cl$comparison==k & cl$stage=='Eligible younger individuals'];stopifnot(length(z)==1L);z},numeric(1))
mi <- mi[match(comparisons,mi$comparison),]
pN <- function(model) vapply(comparisons,function(k){z<-pa$complete_pair_eligible_N[pa$comparison==k & pa$model==model];stopifnot(length(z)==1L);z},numeric(1))
cn <- get_co(clinical); mn <- get_co(microbial)
cgen <- vapply(c('G0','G1','G2'),function(g)get_co(clinical,g),numeric(1))
mgen <- vapply(c('G0','G1','G2'),function(g)get_co(microbial,g),numeric(1))
stopifnot(cn==5660,mn==3217,sum(cgen)==cn,sum(mgen)==mn,
          identical(unname(clinical_N),c(1097L,2262L,1329L)) || all(clinical_N==c(1097,2262,1329)))
# Independently checked against the source 2025 microbiome colData (4,639 rows).
source_microbiome_N <- 4639L
payload <- list(clinical_source=get_n('Clinical source records'),clinical_ultrasound=get_n('Clinical valid ultrasound/generation'),
  clinical_cohort=cn,clinical_families=get_co(clinical,field='families'),clinical_generations=cgen,
  microbiome_matched=source_microbiome_N,microbiome_QC=get_n('Microbiome valid ultrasound/QC'),
  microbiome_cohort=mn,microbiome_families=get_co(microbial,field='families'),microbiome_generations=mgen,
  clinical_C1=clinical_N,adult_CC=mi$complete_case_N,adult_MI=mi$adult_age_sex_complete_N,
  pair_P1=pN('P1'),pair_P2=pN('P2'),B1=sg$N[sg$model=='B1'],adult_B1_B2=sg$N[sg$model=='B2'])
jsonlite::write_json(payload,file.path(audit,'verified_counts.json'),auto_unbox=TRUE,pretty=TRUE)
source_counts <- do.call(rbind,lapply(names(payload),function(k)data.frame(item=k,group=if(length(payload[[k]])==3)c('G0/G0-G1','G1/G1-G2','G2/G0-G2') else 'Overall',N=unname(payload[[k]]))))
write.csv(source_counts,file.path(out,'Source_counts_EN.csv'),row.names=FALSE,fileEncoding='UTF-8')
contract <- list(claim='Two source-data paths lead to distinct clinical and microbiome cohorts and three established analysis directions.',
  archetype='schematic-led composite',backend='R native grid',width_mm=183,height_mm=128,dpi=600,
  font='Arial',body_pt=6.5,minimum_pt=5.8,title_pt=9,maximum_contract_pt=9.01,
  alignment='Single unified flow diagram, no quantitative plot panels. Column headers share a baseline; cohort boxes share width/height and generation anchors.',
  changes='2026-11-01 marked revision: remove both cohort footnotes; replace the dense right column with three concise analysis blocks. No models, sample-count pairs, extra concepts or statistical reanalysis.',
  right_column='Established analysis names only. The three clinical generation comparisons are retained; detailed model-specific sample sizes remain in Source_counts_EN.csv.',
  count_units='Clinical comparisons: younger participants. Similarity: unique participants, not dyads. Analysis sets overlap.')
jsonlite::write_json(contract,file.path(audit,'figure_contract.json'),auto_unbox=TRUE,pretty=TRUE)
fmt <- function(x) format(x,big.mark=',',scientific=FALSE,trim=TRUE)
W<-183; H<-128
ink <- '#171717'; grey <- '#989CA0'; pale <- '#D4D7DA'
gen_cols <- c('#F1B656','#397FC7','#040676')
text <- function(label,x,y,size=6.5,bold=FALSE,just='left',colour=ink,width=NULL) {
  label<-fl_model_labels(label)
  gp<-gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=colour)
  if(!is.null(width)) {
    pushViewport(viewport(gp=gp));actual<-convertWidth(stringWidth(label),'mm',valueOnly=TRUE);popViewport()
    if(actual>width)stop('Text exceeds column: ',label,' (',actual,' > ',width,')')
  }
  grid.text(label,x=unit(x,'mm'),y=unit(H-y,'mm'),just=c(just,'centre'),gp=gp)
}
seg <- function(x,y,x2,y2,col=grey,lwd=.6,arrow=FALSE) grid.lines(unit(c(x,x2),'mm'),unit(H-c(y,y2),'mm'),
  gp=gpar(col=col,lwd=lwd,lineend='butt'),arrow=if(arrow)grid::arrow(length=unit(1.1,'mm'),angle=25,type='closed') else NULL)
cohort_box <- function(top,title,N,families,gens) {
  grid.rect(x=unit(86,'mm'),y=unit(H-(top+21.5),'mm'),width=unit(52,'mm'),height=unit(43,'mm'),gp=gpar(fill='white',col=grey,lwd=.7))
  text(title,86,top+6,7,TRUE,'centre')
  text(paste(fmt(N),'participants'),86,top+14,7,TRUE,'centre')
  text(paste(fmt(families),'primary families'),86,top+19,6.5,FALSE,'centre')
  for(j in 1:3) {
    x<-68+(j-1)*18
    seg(x-3,top+27,x+3,top+27,col=gen_cols[j],lwd=2.6)
    text(paste0('G',j-1),x,top+32,6.5,TRUE,'centre')
    text(paste0('n = ',fmt(gens[j])),x,top+37,6.2,FALSE,'centre')
  }
}
analysis_box <- function(top,height,lines,comparisons_line=NULL) {
  grid.rect(x=unit(151.5,'mm'),y=unit(H-(top+height/2),'mm'),
    width=unit(53,'mm'),height=unit(height,'mm'),
    gp=gpar(fill='white',col=pale,lwd=.7))
  ys <- if(is.null(comparisons_line)) top+height/2+c(-2.5,2.5) else top+c(7,12)
  for(j in seq_along(lines)) text(lines[j],151.5,ys[j],6.8,TRUE,'centre',width=49)
  if(!is.null(comparisons_line)) text(comparisons_line,151.5,top+21,6.4,FALSE,'centre',width=49)
}
draw <- function() {
  grid.newpage()
  text('Study populations and analyses',5,7,9,TRUE)
  text('Source data',5,20,7,TRUE)
  text('Study populations',60,20,7,TRUE)
  text('Analyses',125,20,7,TRUE)
  seg(54.5,25,54.5,123,col=pale,lwd=.6)
  seg(119,25,119,123,col=pale,lwd=.6)
  # Clinical branch: stool availability is not a selection requirement.
  text('Clinical records',26,34,6.8,TRUE,'centre')
  text(paste0('n = ',fmt(payload$clinical_source)),26,40,6.5,FALSE,'centre')
  seg(26,44,26,50,arrow=TRUE)
  text('Valid ultrasound + generation',26,56,6.1,FALSE,'centre')
  text(paste0('n = ',fmt(payload$clinical_ultrasound)),26,62,6.5,FALSE,'centre')
  seg(48,59,57.5,59,arrow=TRUE)
  cohort_box(27,'Clinical cohort',cn,payload$clinical_families,cgen)
  seg(112.5,48.5,123,48.5,arrow=TRUE)
  analysis_box(34.5,28,c('Intergenerational','steatosis associations'),
    'G0–G1   ·   G1–G2   ·   G0–G2')
  # Microbiome branch has its own ultrasound/QC and family eligibility filters.
  text('Clinical–microbiome matched',26,86,6.2,TRUE,'centre')
  text(paste0('n = ',fmt(source_microbiome_N)),26,92,6.5,FALSE,'centre')
  seg(26,96,26,102,arrow=TRUE)
  text('Valid ultrasound + microbiome QC',26,108,5.8,FALSE,'centre')
  text(paste0('n = ',fmt(payload$microbiome_QC)),26,114,6.5,FALSE,'centre')
  seg(48,111,57.5,111,arrow=TRUE)
  cohort_box(79,'Microbiome cohort',mn,payload$microbiome_families,mgen)
  seg(112.5,100.5,116,100.5);seg(116,90.5,116,114.5)
  seg(116,90.5,123,90.5,arrow=TRUE);seg(116,114.5,123,114.5,arrow=TRUE)
  analysis_box(81.5,18,c('Family microbiome','similarity'))
  analysis_box(105.5,18,c('Steatosis–microbiome','associations'))
  # Definitions are retained in Figure_legend_EN.txt, outside the figure.
}
stem <- file.path(out,'Study_population_and_analysis_sets')
grDevices::cairo_pdf(paste0(stem,'.pdf'),width=W/25.4,height=H/25.4,family='Arial');draw();dev.off()
svglite::svglite(paste0(stem,'.svg'),width=W/25.4,height=H/25.4);draw();dev.off()
grDevices::png(paste0(stem,'.png'),width=W,height=H,units='mm',res=600,type='cairo',bg='white');draw();dev.off()
grDevices::tiff(paste0(stem,'.tiff'),width=W,height=H,units='mm',res=600,type='cairo',compression='lzw',bg='white');draw();dev.off()
stopifnot(identical(before,tools::md5sum(protected)))
writeLines(c('Study population and analysis overview.',
 'Clinical and microbiome cohorts require at least two eligible observed generations within a primary-family label. The clinical cohort does not require a stool sample. Minors are retained in all-age analyses.',
 'The right-hand column identifies the established analysis directions. Detailed model-specific sample sizes are provided separately in Source_counts_EN.csv; clinical comparison counts denote younger participants and family-similarity counts denote unique eligible participants, not pairs.',
 'Hepatic steatosis is defined by ultrasound. Analysis samples overlap and cannot be added.'),file.path(out,'Figure_legend_EN.txt'))
message('DRAWING COMPLETE: ',out)
