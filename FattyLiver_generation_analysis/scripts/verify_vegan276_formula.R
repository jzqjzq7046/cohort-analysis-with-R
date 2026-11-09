# Diagnose the protocol difference using the saved full-data reduced residuals.
root<-'E:/YFS/workspace/FattyLiver_generation_analysis'
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
.libPaths(c(file.path(root,'.r-library-sgb'),.libPaths()))
out<-file.path(root,'outputs/24_Two_generation_reanalysis/_running/two_generation_20261024');prefix<-file.path(out,'private/work_kernel_B1_family');N<-3193L;B<-9999L;k<-as.integer(file.info(paste0(prefix,'_Q.bin'))$size/(8*N))
rb<-function(s,w,n,size){con<-file(paste0(prefix,s),'rb');on.exit(close(con));readBin(con,w,n=n,size=size,endian='little')}
E<-matrix(rb('_E.bin','double',N*N,8),N,N);Q<-matrix(rb('_Q.bin','double',N*k,8),N,k);P<-matrix(rb('_P.bin','integer',N*B,4),B,N)
m<-readRDS(file.path(out,'private/stage_prepare.rds'))$metadata;m<-m[complete.cases(m[,c('generation','age','sex','prep_batch','seq_pool','season_sin','season_cos','abx_recent')]),];m<-m[order(m$family_id,m$sample_id),];r<-m$fatty_liver-ave(m$fatty_liver,m$family_id)
Qc<-Q[,-1,drop=FALSE];qx<-Q[,1];num<-as.numeric(crossprod(qx,E%*%qx));den<-sum(diag(E))-num
fit<-list(Ybar=E,CCA=list(rank=1L,qrank=1L,tot.chi=num,QR=qr(cbind(Qc,r)),u=matrix(0,N,1)),CA=list(tot.chi=den,rank=N-k-1,u=matrix(0,N,1)),pCCA=list(rank=ncol(Qc),QR=qr(Qc)),method='dbrda',call=quote(NULL));class(fit)<-c('dbrda','rda','cca')
test<-vegan::permutest(fit,permutations=P[1:10,,drop=FALSE],model='reduced',parallel=1L);df<-N-length(unique(m$family_id))-k
formula<-vapply(1:10,function(i){xi<-r[order(P[i,])];norm<-sum(r^2)-sum(crossprod(Qc,xi)^2);nu<-as.numeric(crossprod(xi,E%*%xi))/norm;nu/((sum(diag(E))-nu)/df)},numeric(1))
z<-data.frame(draw=1:10,current_vegan=as.numeric(test$F.perm)*df/test$df[2],new_formula=formula);z$absolute_error<-abs(z$current_vegan-z$new_formula);print(z,digits=16);stopifnot(max(z$absolute_error)<1e-7)
write.csv(z,file.path(out,'audit/vegan276_formula_verification.csv'),row.names=FALSE);cat('PASS: current vegan re-residualizes exposure against permuted conditions; corrected formula verified.\n')
