# Diagnostic-triggered, bounded effect-stability check; no new discovery claims.
invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
root<-'E:/YFS/workspace/FattyLiver_generation_analysis'
options(fatty_liver.autorun=FALSE)
source(file.path(root,'R/family_candidate_followup.R'),encoding='UTF-8')
source(file.path(root,'R/family_four_stage.R'),encoding='UTF-8');cfg<-f4_setup(root)
source(file.path(root,'R/fig6_family_LMM.R'),encoding='UTF-8')
out<-fc_adopted(root,'51_Fig6_family_mixed_models')
inp<-readRDS(file.path(out,'private/input_population.rds'));r<-readRDS(file.path(out,'private/analysis_result.rds'))
ade_protect(inp$protected)
# Correct descriptive diagnostic labels only; saved fitted estimates are reused.
f6_tables_reports(r,inp,out)
rows<-list();k<-0L
for(mode in c('Primary','Steatosis-conditional'))for(id in inp$candidates$feature_id) {
  original<-r$tests[r$tests$mode==mode&r$tests$phenotype=='WAIST18'&r$tests$feature_id==id,]
  if(!nrow(original))next
  mm<-inp$metadata;y<-f6_transform(mm$WAIST18)
  ok<-is.finite(y)&if(mode=='Steatosis-conditional')!is.na(mm$fatty_liver)else TRUE
  d<-droplevels(mm[ok,]);d$y<-y[ok];d$microbe<-inp$abundance_transformed[id,d$sample_id]
  form<-as.formula(paste('y ~ microbe + age_z + sex + BMI_z + I(BMI_z^2)',
    if(mode=='Steatosis-conditional')'+ fatty_liver'else'', '+ (1 | family_id)'))
  fit<-lmerTest::lmer(form,d,REML=TRUE,control=lme4::lmerControl(check.conv.singular='ignore'))
  stopifnot(!length(fit@optinfo$conv$lme4$messages))
  b<-coef(summary(fit,ddf='Satterthwaite'))['microbe',]
  k<-k+1L;rows[[k]]<-data.frame(mode=mode,feature_id=id,N=nrow(d),
    original_beta=original$beta,quadratic_BMI_beta=b['Estimate'],
    original_CI_low=original$CI_low,original_CI_high=original$CI_high,
    quadratic_CI_low=b['Estimate']-qt(.975,b['df'])*b['Std. Error'],
    quadratic_CI_high=b['Estimate']+qt(.975,b['df'])*b['Std. Error'],
    beta_difference=b['Estimate']-original$beta,same_direction=sign(b['Estimate'])==sign(original$beta),
    original_residual_BMI2_correlation=cor(r$diagnostics$residual_z[r$diagnostics$mode==mode&r$diagnostics$trait_id=='WAIST18'&r$diagnostics$feature_id==id],d$BMI_z^2),
    quadratic_residual_BMI2_correlation=cor(residuals(fit),d$BMI_z^2),
    interpretation='Post hoc diagnostic sensitivity for effect stability only; no new discovery or significance count')
}
tab<-ov_bind(rows);v3_table(tab,'Fig6_waist_BMI_curvature_sensitivity_EN',out)
delta<-max(abs(tab$beta_difference));same<-sum(tab$same_direction)
note<-c('Diagnostic-triggered waist sensitivity: the primary linear-BMI model showed residual association with squared BMI (approximately -0.33). A quadratic BMI term was added in the same participants for the eight eligible displayed SGBs in each panel (16 fits). This check was motivated by model diagnostics after the primary analysis, not prespecified. It assesses coefficient stability only; no additional FDR/significance claims are made.',
  sprintf('Same direction: %s/%s. Maximum absolute beta change: %.4f standardized units. Main Fig.6 and its prespecified complete-background BH values are retained.',same,nrow(tab),delta),
  'Primary model assumptions therefore warrant caution, particularly for waist circumference. Transforming outcomes does not guarantee linear BMI adjustment; the sensitivity is reported alongside the agreed primary model, without silently replacing it.')
writeLines(note,file.path(out,'reports/Diagnostic_sensitivity_EN.txt'),useBytes=TRUE)
zh<-sprintf('诊断补充：腰围模型存在BMI非线性残差；对8个合格候选菌×2个面板加入BMI平方项进行事后系数稳定性检查，%s/%s方向一致，β最大绝对变化%.4f。主图及原全背景BH结果保留；腰围结果需结合该限制解释。此检查不新增显著性结论。',same,nrow(tab),delta)
cat(zh,'\n',file=file.path(out,'reports/Fig6_legend_and_results_ZH.txt'),append=TRUE)
cat(paste(note,collapse='\n'),'\n',file=file.path(out,'reports/Methods_results_limits_EN.md'),append=TRUE)
ade_protect(inp$protected);print(tab);cat(zh,'\n')
