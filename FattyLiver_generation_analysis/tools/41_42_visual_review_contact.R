# R-only private contact sheets for final human layout review; no statistics.
root <- Sys.getenv('FATTYLIVER_ROOT','E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),file.path(R.home(),'library')))
stopifnot(requireNamespace('png',quietly=TRUE))
args <- commandArgs(TRUE); kinds <- c(A='41_Generation_heterogeneity',D='42_Functional_pathways')
if(length(args))kinds<-kinds[names(kinds)%in%args]
for(k in names(kinds)){
  p<-list.files(file.path(root,'outputs',kinds[[k]],'current'),'^Figure_export_manifest_EN.csv$',full.names=TRUE,recursive=TRUE)
  stopifnot(length(p)==1L);out<-dirname(dirname(p));m<-read.csv(p,fileEncoding='UTF-8-BOM')
  pages<-split(m$figure,ceiling(seq_len(nrow(m))/6L))
  for(i in seq_along(pages)){
    target<-file.path(out,'audit',sprintf('Visual_review_contact_%02d.png',i))
    grDevices::png(target,width=2000,height=1600,res=160,type='cairo',bg='white')
    grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,3)))
    for(j in seq_along(pages[[i]])){
      stem<-pages[[i]][j];a<-png::readPNG(file.path(out,'figures',paste0(stem,'.png')))
      grid::pushViewport(grid::viewport(layout.pos.row=ceiling(j/3),layout.pos.col=(j-1)%%3+1))
      ratio<-dim(a)[2]/dim(a)[1]
      grid::grid.raster(a,x=.5,y=.47,width=if(ratio>1.2).96 else .75*ratio,
        height=if(ratio>1.2).96/ratio else .75,interpolate=TRUE)
      grid::grid.text(sub('^Figure','',stem),y=.98,gp=grid::gpar(fontfamily='Arial',fontsize=8))
      grid::popViewport()
    }
    grid::popViewport();grDevices::dev.off()
    cat(normalizePath(target,winslash='/'),'\n')
  }
}
