# Presentation-only definitions for the approved Table 2 header, 2026-10-30.
fl_table2_group <- 'Hepatic steatosis in the younger generation, n/N (%)'
fl_table2_data <- function(raw,prev){
 cells<-function(comparison,exposure){a<-prev[prev$comparison==comparison&prev$exposure==exposure,];stopifnot(nrow(a)==1L);sprintf('%d/%d (%.1f%%)',a$cases,a$N,a$percent)}
 data.frame(`Generational comparison`=raw$comparison,Model=ifelse(raw$model=='C1','Model 1','Unadjusted'),
  `Participants of younger generation, N`=raw$n,
  `No affected older-generation member observed`=vapply(raw$comparison,cells,character(1),exposure='No older case'),
  `≥1 affected older-generation member observed`=vapply(raw$comparison,cells,character(1),exposure='Older case'),
  OR=raw$OR,`95% CI`=sprintf('%.2f–%.2f',raw$CI_low,raw$CI_high),`P value`=raw$p_value,check.names=FALSE)
}
fl_table2_notes <- function()c(
 'Participants refer to G1 in the G0–G1 comparison and G2 in the G1–G2 and G0–G2 comparisons. Values are cases/participants assessed in the younger generation, stratified by whether at least one observed older-generation member in the same primary family had hepatic steatosis.',
 'Unadjusted: No covariate adjustment. Model 1: Adjusted for age and sex of the younger participant and the number of observed older-generation members when this number varied; the latter was constant and omitted for G1–G2.',
 'Both models used logistic generalized estimating equations with primary-family-clustered robust standard errors and the same participants within each comparison. ORs compare participants with versus without an observed affected older-generation member. P values are two-sided and unadjusted for multiple testing.',
 'Hepatic steatosis was defined by ultrasound. Only the two relevant generations were required; minors were retained. Primary family labels do not establish direct parent–child or grandparent–grandchild relationships. OR, odds ratio; CI, confidence interval.')
fl_draw_table2 <- function(d,title,notes,width,height,font,cols,headerheight){
 stopifnot(ncol(d)==8L,abs(sum(cols)-1)<1e-10)
 grid::grid.newpage();grid::pushViewport(grid::viewport(gp=grid::gpar(fontfamily='Times New Roman',fontsize=font)))
 margin<-8;available<-width-2*margin;cw<-available*cols;edges<-c(margin,margin+cumsum(cw));top<-17
 tx<-function(text,x,y,just='centre',size=font,face='plain')grid::grid.text(text,x=grid::unit(x,'mm'),y=grid::unit(height-y,'mm'),just=just,gp=grid::gpar(fontsize=size,fontface=face),default.units='mm')
 line<-function(y,x1=margin,x2=width-margin,lwd=.7)grid::grid.lines(x=grid::unit(c(x1,x2),'mm'),y=grid::unit(rep(height-y,2),'mm'),gp=grid::gpar(lwd=lwd))
 tx(title,width/2,7,size=font+2,face='bold');line(top,lwd=.9)
 split<-top+10
 # The shared heading spans the two disease-status columns only.
 tx(wrap_mm(fl_table2_group,sum(cw[4:5])-3,font),(edges[4]+edges[6])/2,top+5,face='bold')
 line(split,edges[4]+1,edges[6]-1,lwd=.55)
 for(j in seq_len(8)){
  y<-if(j%in%4:5)(split+top+headerheight)/2 else top+headerheight/2
  tx(wrap_mm(names(d)[j],cw[j]-3,font),(edges[j]+edges[j+1])/2,y,face='bold')
 }
 y<-top+headerheight;line(y,lwd=.6)
 for(i in seq_len(nrow(d))){
  values<-vapply(seq_len(8),function(j){v<-d[[j]][i];if(j==3)format(v,big.mark=',',scientific=FALSE)else if(j==6)sprintf('%.2f',v)else if(j==8){if(v<.0001)format(v,scientific=TRUE,digits=3)else sprintf('%.4f',v)}else as.character(v)},character(1))
  rh<-6.1
  for(j in seq_len(8))tx(values[j],(edges[j]+edges[j+1])/2,y+rh/2)
  y<-y+rh
 };line(y,lwd=.9);y<-y+4
 grid::pushViewport(grid::viewport(gp=grid::gpar(fontsize=font-.5)))
 for(note in notes){s<-wrap_mm(note,available,font-.5);n<-length(strsplit(s,'\n',fixed=TRUE)[[1]]);tx(s,margin,y,just=c('left','top'),size=font-.5);y<-y+n*(font-.5)*25.4/72*1.18+2}
 grid::popViewport(2);stopifnot(y<height-3)
}
