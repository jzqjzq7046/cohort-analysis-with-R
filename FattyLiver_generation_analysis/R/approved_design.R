# User-approved revision; numerical kernels retain their independent checks.
fl_ap_design <- "approved_design_v1_20261020"
fl_ap_labels <- c(core_base="Technical reference",core_adjusted="Main adjusted",adult_core="Adult reference",
  adult_lifestyle="Adult lifestyle",adult_BMI="Adult + BMI",bowel_base="Bowel sample reference",
  bowel_adjusted="Bowel + preparation",unperturbed_core="No recent perturbation",DNA_base="DNA sample reference",core_DNA="DNA concentration")
fl_cv_scenarios <- function(m,module="diversity"){
  stopifnot(module %in% c("diversity","similarity"))
  core<-c("seq_pool","survey_center","abx_recent","season_sin","season_cos","freeze_delay_log",if(module=="diversity")"stool_gap_log")
  life<-c("smoking","alcohol","diet_score","physical_activity","energy_log");bowel<-c("constipation","loose_stools","probiotic")
  baseline<-c("age","sex","generation","prep_batch","log_total")
  cc<-function(v)complete.cases(m[,unique(c(baseline,v)),drop=FALSE])
  cm<-cc(core);am<-!is.na(m$age)&m$age>=18&cc(c(core,life,"BMI18"));bm<-cc(c(core,bowel));dm<-cc(c(core,"log_DNA"))
  um<-cm&!is.na(m$gastro_recent)&m$gastro_recent==0&m$abx_recent==0
  make<-function(mask,cov,name)list(mask=mask,covariates=cov,label=unname(fl_ap_labels[name]))
  list(core_base=make(cm,"seq_pool","core_base"),core_adjusted=make(cm,core,"core_adjusted"),
    adult_core=make(am,core,"adult_core"),adult_lifestyle=make(am,c(core,life),"adult_lifestyle"),adult_BMI=make(am,c(core,life,"BMI18"),"adult_BMI"),
    bowel_base=make(bm,core,"bowel_base"),bowel_adjusted=make(bm,c(core,bowel),"bowel_adjusted"),
    unperturbed_core=make(um,core,"unperturbed_core"),DNA_base=make(dm,core,"DNA_base"),core_DNA=make(dm,c(core,"log_DNA"),"core_DNA"))
}
fl_cg_review_specs <- function(nutrition=FALSE){
  core<-c("age","sex","n_older_observed","older_mean_age","older_female_fraction","survey_center","urbanicity")
  life<-c("alcohol","smoking","diet_score","physical_activity","energy_log");full<-c(core,life,"BMI18")
  models<-list(list("R0_crude_same_sample","review_core",character(),FALSE,"all_age","crude_same_sample",core),
    list("R1_review_primary","review_core",core,TRUE,"all_age","all_age_clinical_core",core),
    list("R2_adult_core_same_sample","adult_review",core,FALSE,"adult_nested","adult_same_sample_reference",full),
    list("R3_adult_lifestyle","adult_review",c(core,life),FALSE,"adult_nested","adult_lifestyle_extension",full),
    list("R4_adult_lifestyle_BMI","adult_review",full,FALSE,"adult_nested","adiposity_conditioned_sensitivity",full))
  do.call(rbind,lapply(models,function(x)data.frame(model=x[[1]],sample_set=x[[2]],covariates=paste(x[[3]],collapse="+"),
    primary=x[[4]],comparison_pair=x[[5]],analysis_role=x[[6]],complete_case_fields=paste(x[[7]],collapse="+"),
    references="Long2019;VanderWeele2019;LiangZeger1986",study_adaptation="Approved revision; fixed three-generation primary families; any observed older case; no verified pedigree; family GEE; common samples; no default medication layer")))
}
# Descriptive-only computation: superseded 18-model analysis is not rerun.
fl_cg_compute<-function(m,cfg){
  scope<-fl_cg_prepare_metadata(m);m<-scope$cohort;flows<-desc<-coverage<-list()
  for(g in cfg$cg_pairs){key<-paste(g,collapse="-");b<-fl_cg_make_comparison(m,g);d<-b$data
    flows[[key]]<-b$flow;desc[[key]]<-fl_cg_describe(d,key)
    coverage[[key]]<-do.call(rbind,lapply(sort(unique(d$n_older_observed)),function(k){z<-d[d$n_older_observed==k,]
      data.frame(comparison=key,n_older_observed=k,n=nrow(z),families=length(unique(z$family_id)),older_any_case_n=sum(z$older_any_case),younger_cases=sum(z$fatty_liver))}))}
  list(family_eligibility=scope$summary,cohort_summary=data.frame(n=nrow(m),families=scope$contract$families,cases=sum(m$fatty_liver),rule=scope$contract$rule),
    flow=do.call(rbind,flows),descriptive=do.call(rbind,desc),older_coverage=do.call(rbind,coverage),models=data.frame(status=character()),private=list())
}
fl_ap_table<-function(d){
  esc<-function(x){x<-gsub("&","&amp;",as.character(x),fixed=TRUE);gsub("<","&lt;",x,fixed=TRUE)}
  if(!nrow(d))return("<p>No rows.</p>")
  cols<-lapply(d,function(v)if(is.numeric(v))ifelse(is.na(v),"",format(signif(v,5),trim=TRUE))else ifelse(is.na(v),"",as.character(v)))
  rows<-vapply(seq_len(nrow(d)),function(i)paste0("<tr>",paste0("<td>",esc(vapply(cols,`[`,character(1),i)),"</td>",collapse=""),"</tr>"),character(1))
  paste0("<div class='scroll'><table><tr>",paste0("<th>",esc(names(d)),"</th>",collapse=""),"</tr>",paste(rows,collapse=""),"</table></div>")
}
fl_ap_page<-function(path,title,notes,sections){
  lines<-c("<!doctype html><html lang='zh-CN'><meta charset='utf-8'>",
    "<style>body{font:15px/1.7 Arial,'Microsoft YaHei';max-width:1450px;margin:28px auto;padding:12px;color:#243746}.scroll{overflow:auto}table{border-collapse:collapse;font-size:12px}td,th{border:1px solid #ddd;padding:6px;white-space:nowrap}th{background:#eef2f7}</style>",
    paste0("<title>",title,"</title><h1>",title,"</h1>"),paste0("<p>",notes,"</p>"))
  for(n in names(sections))lines<-c(lines,paste0("<h2>",n,"</h2>"),fl_ap_table(sections[[n]]))
  writeLines(enc2utf8(c(lines,"</html>")),path,useBytes=TRUE)
}
fl_ap_clinical_wide<-function(r){
  do.call(rbind,lapply(c("G0-G1","G1-G2","G0-G2"),function(cmp){
    z<-r[r$comparison==cmp,];a<-z[z$model=="R0_crude_same_sample",];b<-z[z$model=="R1_review_primary",]
    stopifnot(nrow(a)==1L,nrow(b)==1L,a$n==b$n,a$cases==b$cases,a$families==b$families)
    data.frame(comparison=cmp,n=a$n,cases=a$cases,families=a$families,minors=a$minors_n,
      crude_OR=a$OR,crude_CI_low=a$CI_low,crude_CI_high=a$CI_high,crude_P_raw=a$p_value,
      adjusted_OR=b$OR,adjusted_CI_low=b$CI_low,adjusted_CI_high=b$CI_high,adjusted_P_raw=b$p_value,
      crude_status=a$status,adjusted_status=b$status,adjusted_covariates=b$conditioned_fields)
  }))
}
fl_ap_clinical_display<-function(r){
  z<-fl_ap_clinical_wide(r)
  interval<-function(e,l,h)sprintf("%.2f (%.2f–%.2f)",e,l,h)
  p<-function(v)ifelse(v<.0001,formatC(v,format="e",digits=2),formatC(v,format="f",digits=4))
  data.frame(Comparison=z$comparison,`N / cases / families`=paste(z$n,z$cases,z$families,sep=" / "),
    `Crude OR (95%CI)`=interval(z$crude_OR,z$crude_CI_low,z$crude_CI_high),`Crude raw P`=p(z$crude_P_raw),
    `Adjusted OR (95%CI)`=interval(z$adjusted_OR,z$adjusted_CI_low,z$adjusted_CI_high),`Adjusted raw P`=p(z$adjusted_P_raw),check.names=FALSE)
}
fl_cg_html<-function(results,path){
  r<-results$review_models;main<-r[r$sample_set=="review_core",];cols<-c("comparison","model","n","cases","families","minors_n","OR","CI_low","CI_high","p_value","status")
  fl_ap_page(path,"06 脂肪肝跨代关联：确认设计",c("原始P；所有模型家系稳健GEE。粗与主调整严格同样本；成人三层另在共同样本内比较。",
    "主家系标签不证明直接亲子。分离或不可估如实保留，不能视为零关联。"),list("主要表：粗与调整"=fl_ap_clinical_display(r),"原始病例比例"=results$descriptive,
    "补充成人共同样本"=r[r$sample_set=="adult_review",cols],"同样本核对"=results$review_same_sample_audit,"完整模型诊断"=r,"基线"=results$cohort_summary))
  fl_write_csv(main,file.path(dirname(path),"Main_clinical_crude_adjusted.csv"));fl_write_csv(r[r$sample_set=="adult_review",],file.path(dirname(path),"Supplementary_clinical_adult.csv"))
  fl_write_csv(fl_ap_clinical_wide(r),file.path(dirname(path),"Main_clinical_comparison_wide.csv"))
}
fl_cg_review_note<-function(results,path){writeLines(enc2utf8(c("# 06 确认设计结果","",
  "粗模型仅暴露；主要调整年轻年龄/性别、年长人数/平均年龄/女性比例、调查地区和城乡。所有模型家系稳健GEE。",
  "粗/主要模型严格同样本；成人核心、生活方式、BMI严格同样本。不再由药物缺失决定成人样本。药物不默认叠加。",
  "原始双侧P；点95%CI。少事件、分离、不稳定诊断均保留，不作因果/中介解释。",
  "主表见Main_clinical_crude_adjusted.csv；成人见Supplementary_clinical_adult.csv。",
  "完整确认设计见config/approved_design_20261020.md。")),path,useBytes=TRUE)}
fl_ms_report<-function(r,out_dir,version){
  ar<-r$adjusted_results;main<-ar[ar$scenario=="core_adjusted",];show<-c("scenario","comparison","n_people","families","within_pairs","observed_within_mean","observed_between_mean","adjusted_distance_difference","p_raw_two_sided","status")
  notes<-c("主要调整双方年龄/性别、调查地区、建库批次、测序池、输入矩阵总量、季节、记录冷冻间隔、近期抗生素及同地区/同批次/同测序池。",
    "9999次Freedman–Lane学生化t残差MRQAP；角色内性别×5岁块节点置换；原始双侧P。配对数非独立样本量，不提供独立配对CI。",
    "负系数代表家系内距离较小，不推断菌株传播。成人、排便、DNA均有同样本参照。",
    paste0("输入：",version,"；",fl_ap_design))
  fl_ap_page(file.path(out_dir,"microbiome_family_similarity.html"),"07 家系菌群相似性：确认设计",notes,
    list("主要调整结果"=main[,show],"主要同样本IF/BF原始描述"=r$adjusted_pair_summary[r$adjusted_pair_summary$scenario=="core_adjusted",],
      "补充调整"=ar[ar$scenario!="core_adjusted",show],"完整数值与诊断"=ar,"可交换范围"=r$adjusted_exchangeability))
  fl_write_csv(main,file.path(out_dir,"Main_similarity_adjusted.csv"));fl_write_csv(ar[ar$scenario!="core_adjusted",],file.path(out_dir,"Supplementary_similarity.csv"))
  writeLines(enc2utf8(c("# 07 确认设计结果","",notes)),file.path(out_dir,"results_zh.md"),useBytes=TRUE)
}
fl_dv_review_report<-function(review,out,covariate_version){
  a<-review$Shannon_models;b<-review$composition_association;main<-a[a$scenario=="core_adjusted",];bm<-b[b$scenario=="core_adjusted",]
  notes<-c("主要调整各代年龄样条3df、性别、代际、调查地区、建库批次、测序池、log10输入矩阵总量、近期抗生素、季节周期、记录冷冻间隔和粪便/超声时间差。",
    "Shannon家系稳健GEE，分代差值及95%CI，原始双侧P。组成：Lingoes Bray–Curtis partial dbRDA，Condition家系和协变量，家系内残差置换，原始F上尾P。",
    "Shannon边际均值差与组成家庭内条件关联目标不同。PCoA为原始距离坐标，未声称去除协变量。",
    "全龄保留未成年人；成人生活方式/BMI严格同样本；排便及DNA另有同样本基线；排除干扰为限制样本检验。",
    "排便及乳酸菌制剂移出主模型。不默认叠加糖尿病/药物/替代营养模型。矩阵总量不称已核实reads数。",
    paste0("协变量版本：",covariate_version,"；",fl_ap_design))
  fl_ap_page(file.path(out,"microbiome_steatosis_diversity.html"),"08 脂肪肝与菌群多样性：确认设计",notes,
    list("主要Shannon"=main,"主要条件组成"=bm,"分代样本"=review$sample_counts,"补充Shannon"=a[a$scenario!="core_adjusted",],
      "补充组成"=b[b$scenario!="core_adjusted",],"异质性"=review$generation_heterogeneity,"置换信息"=review$exchangeability))
  fl_write_csv(main,file.path(out,"Main_Shannon_adjusted.csv"));fl_write_csv(bm,file.path(out,"Main_composition_adjusted.csv"))
  writeLines(enc2utf8(c("# 08 确认设计结果","",notes)),file.path(out,"results_zh.md"),useBytes=TRUE)
}
fl_dv_run<-function(cfg){
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-file.path(ctx$run_dir,"diversity");audit<-file.path(ctx$run_dir,"audit");private<-file.path(ctx$run_dir,"private")
  for(p in c(out,audit,private))dir.create(p)
  writeLines("RUNNING",file.path(ctx$run_dir,"run_status.txt"))
  tryCatch({
    d<-fl_load_processed_microbiome(cfg$project_root,version=cfg$dv_data_version);reviewed<-fl_cv_load(cfg,d$metadata)
    m<-fl_dv_prepare(reviewed$metadata,d$full$relative,d$full$counts);D<-vegan::vegdist(t(d$full$relative),method="bray")
    scenarios<-fl_cv_scenarios(m);mask<-scenarios$core_adjusted$mask
    fl_write_csv(fl_dv_descriptive(m[mask,]),file.path(out,"Table2_diversity_descriptive.csv"))
    coord<-vegan::wcmdscale(as.dist(as.matrix(D)[mask,mask]),k=2,add="lingoes",eig=TRUE)
    saveRDS(list(metadata=m[mask,],coordinates=coord$points,eigenvalues=coord$eig),file.path(private,"pcoa_display.rds"))
    review<-fl_dv_review_compute(m,D,cfg,private,out,audit);fl_dv_review_report(review,out,reviewed$version)
    fl_write_csv(data.frame(covariate_version=reviewed$version),file.path(audit,"covariate_version.csv"))
    paths<-unique(c(file.path(cfg$processed_microbiome_root,d$version,"microbiome_processed.rds"),reviewed$source_files,fl_dv_code_paths(cfg),
      file.path(cfg$project_root,c("R/approved_design.R","config/approved_design_20261020.md"))))
    fl_write_csv(data.frame(path=paths,md5=unname(tools::md5sum(paths))),file.path(audit,"input_code_fingerprints.csv"))
    fl_write_csv(data.frame(scenario=names(scenarios),label=vapply(scenarios,`[[`,character(1),"label"),
      covariates=vapply(scenarios,function(z)paste(z$covariates,collapse=";"),character(1))),file.path(out,"model_definitions.csv"))
    capture.output(sessionInfo(),file=file.path(audit,"sessionInfo.txt"));writeLines("COMPLETED",file.path(ctx$run_dir,"run_status.txt"))
    final<-fl_output_complete_run(ctx);cat("Completed approved 08: ",final,"\n",sep="");invisible(review)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(e)})
}

# Manuscript drawing is centralized in entry 11; do not emit duplicate legacy plots.
fl_dv_figures<-function(m,primary,main_fit,out,audit){
  e<-data.frame(family=m$family_id,generation=m$generation,residual=m$shannon-as.numeric(fitted(main_fit)))
  z<-aggregate(e$residual,list(generation=e$generation),function(v)c(n=length(v),mean=mean(v),SD=sd(v)))
  fl_write_csv(z,file.path(audit,"primary_residual_summary.csv"))
  invisible(NULL)
}
