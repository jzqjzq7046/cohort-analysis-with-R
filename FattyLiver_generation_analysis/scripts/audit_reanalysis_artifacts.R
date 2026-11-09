# Final-size sources remain untouched; these sheets are audit-only previews.
root <- 'E:/YFS/workspace/FattyLiver_generation_analysis'
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),.libPaths()))
run <- readLines(file.path(root,'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt'),n=1)
dirs <- file.path(root,'outputs/24_Two_generation_reanalysis',c('_running','current'),run)
out <- dirs[dir.exists(dirs)];stopifnot(length(out)==1L)
files <- list.files(file.path(out,'figures'),pattern='[.]png$',full.names=TRUE)
pages <- split(files,ceiling(seq_along(files)/4))
for(i in seq_along(pages)){
 grDevices::png(file.path(out,'audit',sprintf('contact_sheet_%02d.png',i)),width=1800,height=1400,type='cairo',bg='white')
 grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,2)))
 for(j in seq_along(pages[[i]])){
  f<-pages[[i]][j];im<-png::readPNG(f,native=TRUE);grid::pushViewport(grid::viewport(layout.pos.row=ceiling(j/2),layout.pos.col=1+(j-1)%%2))
  grid::grid.text(basename(f),y=.98,gp=grid::gpar(fontsize=11))
  aspect<-nrow(im)/ncol(im);w<-.96;h<-.90
  if(aspect<.90/.96*1400/1800)h<-w*aspect*1800/1400 else w<-h/aspect*1400/1800
  grid::grid.raster(im,x=.5,y=.47,width=grid::unit(w,'npc'),height=grid::unit(h,'npc'),interpolate=TRUE)
  grid::popViewport()
 };grid::popViewport();grDevices::dev.off()
}
cat('Audit contact sheets:',length(pages),'for',length(files),'figures\n')
