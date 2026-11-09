# Read-only ABI/package inspection; no data fitting or dependency installation.
root<-Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis')
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),file.path(root,'.r-library'),file.path(R.home(),'library'),.libPaths()))
cat(R.version.string,'\n');print(.libPaths())
for(p in c('Matrix','lme4','Maaslin2')){
  cat(p,as.character(packageVersion(p)),find.package(p),'\n')
  print(packageDescription(p)[['Built']])
}
warnings_seen<-character()
withCallingHandlers(requireNamespace('lme4',quietly=TRUE),warning=function(w){
  warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')})
print(Matrix::Matrix.Version())
cat('lme4 compiled dependency ABI:\n');print(lme4:::.Matrix.abi.build.version)
cat('Namespace-load warnings:\n');print(warnings_seen)
stopifnot(identical(as.character(Matrix::Matrix.Version()$abi),as.character(lme4:::.Matrix.abi.build.version)),
  !any(grepl('ABI version mismatch',warnings_seen)))
cat('PASS: project Matrix runtime and lme4 compiled ABI match.\n')
