# Independent GLS and published-output verification; does not refit discovery.
invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
root<-'E:/YFS/workspace/FattyLiver_generation_analysis'
options(fatty_liver.autorun=FALSE)
source(file.path(root,'R/family_candidate_followup.R'),encoding='UTF-8')
source(file.path(root,'R/family_four_stage.R'),encoding='UTF-8');cfg<-f4_setup(root)
base<-file.path(root,'outputs/51_Fig6_family_mixed_models')
if(file.exists(file.path(base,'LATEST_RUN.txt')))out<-fc_adopted(root,'51_Fig6_family_mixed_models')else {
  out<-list.dirs(file.path(base,'current'),recursive=FALSE,full.names=TRUE)
  stopifnot(length(out)==1L,readLines(file.path(out,'run_status.txt'),1)=='COMPLETED')
}
inp<-readRDS(file.path(out,'private/input_population.rds'))
r<-readRDS(file.path(out,'private/analysis_result.rds'))
ade_protect(inp$protected)
tests<-read.csv(file.path(out,'tables/Fig6_tests_EN.csv'),fileEncoding='UTF-8-BOM',check.names=FALSE)
stopifnot(nrow(tests)==2*372*9,nrow(inp$metadata)==2792,nlevels(inp$metadata$family_id)==1104)
checks<-list();j<-0L
for(mode in c('Primary','Steatosis-conditional'))for(field in c('BMI18','ALAT18')) {
  id<-inp$candidates$feature_id[1]
  z<-tests[tests$mode==mode&tests$phenotype==field&tests$feature_id==id,];stopifnot(nrow(z)==1)
  mm<-inp$metadata;y0<-mm[[field]]
  ok<-is.finite(y0)&if(mode=='Steatosis-conditional')!is.na(mm$fatty_liver)else TRUE
  v<-y0[ok];v<-qnorm((rank(v,ties.method='average')-.5)/length(v));yy<-(v-mean(v))/sd(v)
  mm<-droplevels(mm[ok,]);mm$microbe<-inp$abundance_transformed[id,mm$sample_id]
  form<-reformulate(c('microbe','age_z','sex',if(field!='BMI18')'BMI_z',if(mode=='Steatosis-conditional')'fatty_liver'))
  W<-model.matrix(form,mm)
  # Sherman-Morrison inverse for random-intercept covariance in every family.
  inv<-function(B) {
    if(is.null(dim(B)))B<-matrix(B,ncol=1)
    ans<-B/z$residual_variance
    for(ii in split(seq_len(nrow(mm)),mm$family_id)) {
      adjustment<-z$family_variance/(z$residual_variance*(z$residual_variance+length(ii)*z$family_variance))
      ans[ii,]<-ans[ii,,drop=FALSE]-matrix(rep(adjustment*colSums(B[ii,,drop=FALSE]),each=length(ii)),nrow=length(ii))
    }
    ans
  }
  covb<-solve(crossprod(W,inv(W)));bb<-covb%*%crossprod(W,inv(yy))
  est<-bb['microbe',1];se<-sqrt(covb['microbe','microbe'])
  p<-2*pt(-abs(est/se),df=z$df);ci<-est+c(-1,1)*qt(.975,z$df)*se
  delta<-max(abs(c(est-z$beta,se-z$SE,p-z$P_raw,ci[1]-z$CI_low,ci[2]-z$CI_high)))
  stopifnot(delta<1e-8,nrow(mm)==z$N,nlevels(mm$family_id)==z$families)
  j<-j+1L;checks[[j]]<-data.frame(mode=mode,phenotype=field,feature_id=id,
    check='Independent random-intercept GLS beta/SE, t CI and P',max_abs_difference=delta,passed=TRUE)
}
for(mode in unique(tests$mode)) {
  z<-tests[tests$mode==mode,];good<-z$status=='ESTIMATED'&is.finite(z$P_raw)
  pp<-rep(1,nrow(z));pp[good]<-z$P_raw[good]
  stopifnot(max(abs(p.adjust(pp,'BH')[good]-z$P_BH[good]))<1e-12)
}
v3_table(ov_bind(checks),'Fig6_independent_numerical_validation_EN',out)
ade_protect(inp$protected)
cat('NUMERICAL_VALIDATION_PASSED\n');print(ov_bind(checks))
print(read.csv(file.path(out,'tables/Fig6_candidate_association_counts_EN.csv'),fileEncoding='UTF-8-BOM'))
print(read.csv(file.path(out,'tables/Fig6_flow_EN.csv'),fileEncoding='UTF-8-BOM'))
cat('DIAGNOSTICS\n')
dg<-read.csv(file.path(out,'tables/Fig6_candidate_model_diagnostics_EN.csv'),fileEncoding='UTF-8-BOM')
print(summary(dg[,c('residual_skewness','residual_excess_kurtosis','residual_normal_QQ_correlation','abs_residual_fitted_spearman','residual_age_squared_correlation','residual_BMI_squared_correlation')]))
cat('SINGULAR/CONVERGENCE\n');print(with(tests,table(mode,singular,converged)))
cat('FAILED\n');print(tests[tests$status!='ESTIMATED',c('feature_id','phenotype','mode','warning')])
