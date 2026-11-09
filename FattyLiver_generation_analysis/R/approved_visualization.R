# R-only quantitative grids. Real adopted results; no model refits or resampling.
# Figure contract and variable mapping: config/approved_design_20261020.md.
fl_viz_palette<-c("#CC247C","#E95351","#F7A24F","#FBEB66","#4EA660","#79CAFB","#5292F7","#AA77E9")
fl_viz_generation<-c(G0="#CC247C",G1="#5292F7",G2="#4EA660")
fl_viz_comparison<-setNames(unname(fl_viz_generation),c("G0-G1","G1-G2","G0-G2"))
fl_viz_comparisons<-names(fl_viz_comparison)
nf_status_colors<-c("No steatosis"="#AA77E9","Steatosis"="#E95351")
nf_shapes<-c(G0=21,G1=22,G2=24)
nf_palette_scale<-function()ggplot2::scale_fill_manual(values=fl_viz_generation)
fl_av_read<-function(p)read.csv(p,fileEncoding="UTF-8-BOM",check.names=FALSE,stringsAsFactors=FALSE,na.strings=c("NA",""))
fl_av_csv<-function(d,p)write.csv(d,p,row.names=FALSE,na="",fileEncoding="UTF-8")
fl_av_p<-function(p)ifelse(is.na(p),"NE",ifelse(p<.0001,formatC(p,format="e",digits=2),formatC(p,format="f",digits=4)))
fl_av_landing<-function(root,module,id,target="index.html"){
  url<-paste0("current/",id,"/",target)
  writeLines(paste0("<!doctype html><meta charset='utf-8'><meta http-equiv='refresh' content='0;url=",url,"'><a href='",url,"'>打开最新结果</a>"),file.path(root,"outputs",module,"index.html"),useBytes=TRUE)
}
fl_av_current<-function(root,module,approved=FALSE){
  runs<-list.dirs(file.path(root,"outputs",module,"current"),recursive=FALSE,full.names=TRUE)
  if(length(runs)!=1L||!file.exists(file.path(runs,"run_status.txt"))||readLines(file.path(runs,"run_status.txt"),n=1)!="COMPLETED")stop("No adopted completed run: ",module)
  if(approved&&!grepl("approved_design_v1",basename(runs),fixed=TRUE))stop("Superseded design in ",module)
  runs
}
fl_av_inputs<-function(root){
  paths<-list();add<-function(key,module,file,approved=FALSE){p<-file.path(fl_av_current(root,module,approved),file)
    if(!file.exists(p))stop("Required result missing: ",p);paths[[key]]<<-p}
  add("table1","01_Table1","tables/Table1_by_generation.csv")
  add("proportions","02_Figure1","figures/Figure1_proportion_source_data.csv")
  add("flow","02_Figure1","figures/sample_flow.csv");add("age","02_Figure1","figures/Figure1_age_source_data.csv")
  add("clinical","06_Cross_generation_clinical","clinical/review_models.csv",TRUE)
  add("cells","06_Cross_generation_clinical","clinical/review_model_cells.csv",TRUE)
  add("clinical_samples","06_Cross_generation_clinical","clinical/review_same_sample_audit.csv",TRUE)
  add("distance","07_Cross_generation_microbiome_similarity","similarity/adjusted_results.csv",TRUE)
  add("pair","07_Cross_generation_microbiome_similarity","similarity/adjusted_pair_summary.csv",TRUE)
  add("shannon","08_Microbiome_steatosis_diversity","diversity/review_Shannon_models.csv",TRUE)
  add("composition","08_Microbiome_steatosis_diversity","diversity/review_composition_association.csv",TRUE)
  add("samples","08_Microbiome_steatosis_diversity","diversity/review_sample_counts.csv",TRUE)
  add("heterogeneity","08_Microbiome_steatosis_diversity","diversity/review_generation_heterogeneity.csv",TRUE)
  add("distribution","08_Microbiome_steatosis_diversity","private/review_core_adjusted/alpha_model.rds",TRUE)
  add("pcoa","08_Microbiome_steatosis_diversity","private/pcoa_display.rds",TRUE)
  d<-lapply(paths[!names(paths)%in%c("distribution","pcoa")],fl_av_read)
  alpha<-readRDS(paths$distribution);m<-alpha$metadata
  d$raw<-m[,c("generation","fatty_liver","shannon")];d$raw$status<-factor(ifelse(m$fatty_liver==1,"Steatosis","No steatosis"),levels=names(nf_status_colors))
  co<-readRDS(paths$pcoa);d$pcoa<-data.frame(co$coordinates,generation=co$metadata$generation,fatty_liver=co$metadata$fatty_liver)
  names(d$pcoa)[1:2]<-c("PC1","PC2");d$pcoa$status<-factor(ifelse(d$pcoa$fatty_liver==1,"Steatosis","No steatosis"),levels=names(nf_status_colors))
  d$pcoa_percent<-100*head(co$eigenvalues,2)/sum(co$eigenvalues[co$eigenvalues>0])
  stopifnot(nrow(d$raw)==sum(d$samples$n[d$samples$scenario=="core_adjusted"]),
    sum(d$raw$fatty_liver)==sum(d$samples$cases[d$samples$scenario=="core_adjusted"]),all(d$clinical_samples$same_participants))
  for(key in c("shannon","composition","distance")){
    z<-d[[key]];if(!setequal(unique(z$scenario),names(fl_ap_labels)))stop("Unexpected approved scenarios: ",key)
    p<-z[[if(key=="distance")"p_raw_two_sided" else "p_raw"]]
    if(any(!is.finite(p)|p<=0|p>1))stop("Invalid permutation/model P: ",key)
  }
  for(v in c("shannon","clinical")){
    z<-d[[v]];i<-grepl("^ESTIMATED",z$status)
    est<-if(v=="clinical")z$OR else z$estimate
    if(any(!is.finite(est[i])|!is.finite(z$CI_low[i])|!is.finite(z$CI_high[i])|z$CI_low[i]>est[i]|z$CI_high[i]<est[i]))stop("Invalid estimated CI: ",v)
    if(v=="clinical"&&any(est[i]<=0|z$CI_low[i]<=0|z$CI_high[i]<=0))stop("Log OR axes require strictly positive estimates and intervals")
  }
  list(data=d,manifest=data.frame(key=names(paths),path=unlist(paths),md5=unname(tools::md5sum(unlist(paths)))))
}

fl_av_main_clinical<-function(d){
  c<-d$cells[d$cells$model=="R0_crude_same_sample",];raw<-do.call(rbind,lapply(split(c,interaction(c$comparison,c$older_any_case,drop=TRUE)),function(z)
    data.frame(comparison=z$comparison[1],older=if(z$older_any_case[1]==1)"Older case" else "No older case",n=sum(z$n),cases=sum(z$n[z$younger_case==1]),pct=100*sum(z$n[z$younger_case==1])/sum(z$n))))
  raw$comparison<-factor(raw$comparison,levels=fl_viz_comparisons);raw$older<-factor(raw$older,levels=c("No older case","Older case"))
  a<-ggplot(raw,aes(comparison,pct,fill=older,shape=older))+geom_point(position=position_dodge(width=.35),size=2,stroke=.3)+
    scale_fill_manual(values=setNames(c("#AA77E9","#E95351"),levels(raw$older)))+scale_shape_manual(values=c(21,24))+
    scale_y_continuous(limits=c(0,40),breaks=c(0,10,20,30,40))+labs(x=NULL,y="Younger steatosis (%)",title="Observed prevalence")
  z<-d$clinical[d$clinical$review_primary,];z$comparison<-factor(z$comparison,levels=rev(fl_viz_comparisons))
  b<-ggplot(z,aes(OR,comparison,colour=comparison,fill=comparison))+geom_vline(xintercept=1,linetype=2,linewidth=.25,colour="#777777")+
    geom_segment(aes(x=CI_low,xend=CI_high,yend=comparison),linewidth=.4)+geom_point(shape=23,size=2,stroke=.3)+
    scale_colour_manual(values=fl_viz_comparison,guide="none")+scale_fill_manual(values=fl_viz_comparison,guide="none")+
    scale_x_log10(limits=c(.2,16),breaks=c(.2,1,5,10))+labs(x="Adjusted OR (log scale)",y=NULL,title="Main association")
  list(plot=nf_patch(list(a,b),"AB",collect_guides=TRUE),source=raw)
}
fl_av_family<-function(d){
  z<-d$pair[d$pair$scenario=="core_adjusted",];z$comparison<-factor(z$comparison,levels=rev(fl_viz_comparisons));z$pair_type<-factor(z$pair_type,levels=c("BF","IF"))
  a<-ggplot(z,aes(mean_distance,comparison,fill=pair_type,shape=pair_type))+geom_point(position=position_dodge(width=.38),size=2,stroke=.3)+
    scale_fill_manual(values=c(BF="#AA77E9",IF="#F7A24F"))+scale_shape_manual(values=c(BF=21,IF=24))+
    labs(x="Mean Bray-Curtis distance",y=NULL,title="Observed similarity")
  q<-d$distance[d$distance$scenario=="core_adjusted",];q$comparison<-factor(q$comparison,levels=rev(fl_viz_comparisons))
  b<-ggplot(q,aes(adjusted_distance_difference,comparison,fill=comparison))+geom_vline(xintercept=0,linetype=2,linewidth=.25,colour="#777777")+
    geom_point(shape=23,size=2,stroke=.3)+scale_fill_manual(values=fl_viz_comparison,guide="none")+
    scale_x_continuous(expand=expansion(mult=c(.15,.15)))+labs(x="Adjusted IF - BF distance",y=NULL,title="Conditional similarity")
  nf_patch(list(a,b),"AB",collect_guides=TRUE)
}
fl_av_shannon<-function(d){
  a<-ggplot(d$raw,aes(generation,shannon,fill=status))+geom_boxplot(width=.58,outlier.shape=NA,linewidth=.25,position=position_dodge(width=.7))+
    geom_point(aes(shape=status),position=position_jitterdodge(jitter.width=.12,dodge.width=.7,seed=82),size=.35,stroke=.1,alpha=.22)+
    scale_fill_manual(values=nf_status_colors)+scale_shape_manual(values=c(21,24))+labs(x=NULL,y="Shannon",title="Observed diversity")
  z<-d$shannon[d$shannon$scenario=="core_adjusted",];z$generation<-factor(z$generation,levels=c("G2","G1","G0"))
  b<-ggplot(z,aes(estimate,generation,colour=generation,fill=generation))+geom_vline(xintercept=0,linetype=2,linewidth=.25,colour="#777777")+
    geom_segment(aes(x=CI_low,xend=CI_high,yend=generation),linewidth=.4)+geom_point(shape=23,size=2,stroke=.3)+
    scale_colour_manual(values=fl_viz_generation,guide="none")+scale_fill_manual(values=fl_viz_generation,guide="none")+
    labs(x="Adjusted Shannon difference",y=NULL,title="Main association")
  nf_patch(list(a,b),"AB",collect_guides=TRUE)
}
fl_av_composition<-function(d){
  a<-ggplot(d$pcoa,aes(PC1,PC2,fill=status,shape=generation))+geom_point(size=.8,stroke=.14,alpha=.45)+
    scale_fill_manual(values=nf_status_colors)+scale_shape_manual(values=nf_shapes)+
    labs(x=sprintf("PCoA 1 (%.1f%%)",d$pcoa_percent[1]),y=sprintf("PCoA 2 (%.1f%%)",d$pcoa_percent[2]),title="Observed composition")+
    guides(fill=guide_legend(override.aes=list(shape=21,size=1.5,alpha=1)),shape=guide_legend(override.aes=list(fill="white",size=1.5,alpha=1)))
  z<-d$composition[d$composition$scenario=="core_adjusted",];z$scheme<-factor(z$scheme,levels=c("family_generation","family"),labels=c("Within family + generation","Within family"))
  b<-ggplot(z,aes(100*R2_total_corrected,scheme))+geom_col(fill="#F7A24F",width=.45,colour="#333333",linewidth=.25)+
    scale_x_continuous(expand=expansion(mult=c(0,.12)))+labs(x="Explained total inertia (%)",y=NULL,title="Conditional association")
  nf_patch(list(a,b),"AB",collect_guides=TRUE)
}
fl_av_supp_clinical<-function(d){
  models<-c("R2_adult_core_same_sample","R3_adult_lifestyle","R4_adult_lifestyle_BMI");labels<-c("Adult reference","Adult lifestyle","Adult + BMI")
  ps<-lapply(fl_viz_comparisons,function(cmp){z<-d$clinical[d$clinical$comparison==cmp&d$clinical$model%in%models,];z<-z[match(models,z$model),]
    z$label<-labels;z$label[!grepl("^ESTIMATED",z$status)]<-paste0(z$label[!grepl("^ESTIMATED",z$status)]," [NE]")
    z$label<-factor(z$label,levels=rev(z$label));v<-z[grepl("^ESTIMATED",z$status),]
    ggplot(z,aes(OR,label))+geom_vline(xintercept=1,linetype=2,colour="#777777",linewidth=.25)+
      geom_segment(data=v,aes(x=CI_low,xend=CI_high,yend=label),linewidth=.35,colour=fl_viz_comparison[cmp])+
      geom_point(data=v,shape=21,size=1.8,fill=fl_viz_comparison[cmp],stroke=.3)+
      scale_y_discrete(limits=rev(as.character(z$label)),drop=FALSE)+
      scale_x_log10(limits=c(.1,30),breaks=c(.1,1,10))+labs(x="OR (log scale)",y=NULL,title=paste0(cmp,"; n=",z$n[1],", cases=",z$cases[1]))})
  nf_patch(ps,"ABC")
}
fl_av_supp_family<-function(d){
  lim<-range(c(0,d$distance$adjusted_distance_difference));lim<-lim+c(-1,1)*diff(lim)*.04
  ps<-lapply(fl_viz_comparisons,function(cmp){z<-d$distance[d$distance$comparison==cmp,];z$label<-factor(fl_ap_labels[z$scenario],levels=rev(unname(fl_ap_labels)))
    ggplot(z,aes(adjusted_distance_difference,label))+geom_vline(xintercept=0,linetype=2,linewidth=.25,colour="#777777")+
      geom_point(shape=21,fill=fl_viz_comparison[cmp],size=1.6,stroke=.3)+scale_x_continuous(limits=lim,breaks=pretty(lim,n=3))+
      labs(x="Adjusted IF - BF",y=NULL,title=cmp)})
  nf_patch(ps,"ABC")
}
fl_av_supp_shannon<-function(d){
  lim<-range(c(0,d$shannon$CI_low,d$shannon$CI_high));lim<-lim+c(-1,1)*diff(lim)*.04
  ps<-lapply(c("G0","G1","G2"),function(g){z<-d$shannon[d$shannon$generation==g,];z$label<-factor(fl_ap_labels[z$scenario],levels=rev(unname(fl_ap_labels)))
    ggplot(z,aes(estimate,label))+geom_vline(xintercept=0,linetype=2,linewidth=.25,colour="#777777")+
      geom_segment(aes(x=CI_low,xend=CI_high,yend=label),linewidth=.3,colour=fl_viz_generation[g])+
      geom_point(shape=nf_shapes[g],fill=fl_viz_generation[g],size=1.6,stroke=.3)+scale_x_continuous(limits=lim,breaks=pretty(lim,n=3))+
      labs(x="Shannon difference",y=NULL,title=g)})
  nf_patch(ps,"ABC")
}
fl_av_supp_composition<-function(d){
  ps<-lapply(c("family","family_generation"),function(s){z<-d$composition[d$composition$scheme==s,];z$label<-factor(fl_ap_labels[z$scenario],levels=rev(unname(fl_ap_labels)))
    ggplot(z,aes(100*R2_total_corrected,label))+geom_col(width=.5,fill=if(s=="family")"#F7A24F" else "#79CAFB",colour="#333333",linewidth=.2)+
      scale_x_continuous(expand=expansion(mult=c(0,.1)))+labs(x="Explained total inertia (%)",y=NULL,title=if(s=="family")"Within family" else "Within family and generation")})
  nf_patch(ps,"AB")
}

fl_av_summary<-function(d){
  c<-d$clinical[d$clinical$review_primary,];a<-d$shannon[d$shannon$scenario=="core_adjusted",];b<-d$composition[d$composition$scenario=="core_adjusted",];s<-d$distance[d$distance$scenario=="core_adjusted",]
  line<-c("# 确认设计重分析结果 / Approved-design results","",paste0("版本：",fl_ap_design),
    "基线319个主家系、1275人。主家系三代资格在子集分析前固定，保留未成年人。所有P为原始P；横断面关联不作因果、传播或中介结论。","",
    "## 脂肪肝跨代关联（主要GEE）","","| 比较 | 人数/病例 | OR（95%CI） | 原始P |","|---|---:|---|---:|")
  for(i in seq_len(nrow(c)))line<-c(line,sprintf("| %s | %d/%d | %.3f（%.3f–%.3f） | %s |",c$comparison[i],c$n[i],c$cases[i],c$OR[i],c$CI_low[i],c$CI_high[i],fl_av_p(c$p_value[i])))
  line<-c(line,"","粗与主调整同样本；成人核心、生活方式、BMI另用共同样本。不可估状态完整保留，不将NE解释为无关联。",
    sprintf("成人扩展共%d项，其中%d项因分离等原因不可估；不能据此确认充分生活方式/BMI调整后的跨代关联。",sum(d$clinical$sample_set=="adult_review"),sum(d$clinical$sample_set=="adult_review"&!grepl("^ESTIMATED",d$clinical$status))),"",
    "## 家系内外菌群相似性（主要MRQAP）","","| 比较 | 人数 | IF均值 | BF均值 | 调整IF−BF | 原始P |","|---|---:|---:|---:|---:|---:|")
  for(i in seq_len(nrow(s)))line<-c(line,sprintf("| %s | %d | %.4f | %.4f | %.4f | %s |",s$comparison[i],s$n_people[i],s$observed_within_mean[i],s$observed_between_mean[i],s$adjusted_distance_difference[i],fl_av_p(s$p_raw_two_sided[i])))
  line<-c(line,"","配对数不是独立样本数；无独立配对OLS置信区间。9999次置换的add-one P最小为0.0001；若多项等于此值，表示均达到本次置换分辨率，并非真实P完全相同。","","## 脂肪肝与Shannon","","| 代际 | 人数/病例 | 差值（95%CI） | 原始P |","|---|---:|---|---:|")
  for(i in seq_len(nrow(a)))line<-c(line,sprintf("| %s | %d/%d | %.3f（%.3f–%.3f） | %s |",a$generation[i],a$n[i],a$cases[i],a$estimate[i],a$CI_low[i],a$CI_high[i],fl_av_p(a$p_raw[i])))
  line<-c(line,"","## 整体菌群组成（家庭内条件关联）","")
  for(i in seq_len(nrow(b)))line<-c(line,sprintf("- %s：R²_total=%.5f%%，原始P=%s。",b$scheme[i],100*b$R2_total_corrected[i],fl_av_p(b$p_raw[i])))
  c(line,"","Shannon与组成估计目标不同。PCoA展示原始距离，不代表坐标已去除协变量。模型、样本与敏感性完整表见同目录。",
    "正文图仅放原始情况和主要调整；成人、排便、DNA和近期扰动检验分到Supplementary。此设计在旧结果已知后修订，未按P挑选模型。")
}

fl_ap_summary_main<-function(args=character(),root){
  cfg<-fl_config(root);cfg$run_ids['results_visualization']<-"approved_results_summary_v1";opt<-fl_entry_options(args,cfg,"results_visualization");cfg<-opt$cfg
  if(opt$mode=="--help"){cat("10: --run collects approved tables and numerical summaries; --check validates adopted inputs. No models are refitted.\n");return(invisible(NULL))}
  if(opt$mode=="--self-test")stop("Use tests/approved_design_contract.R and tests/validate_approved_results.R")
  input<-fl_av_inputs(root);d<-input$data
  if(opt$mode=="--check"){cat("10 approved input and schema checks passed\n");return(invisible(input$manifest))}
  ctx<-fl_output_start_run(cfg);out<-ctx$run_dir;for(p in c("tables","reports","audit"))dir.create(file.path(out,p))
  writeLines("RUNNING",file.path(out,"run_status.txt"))
  tryCatch({
    for(k in names(d)[vapply(d,is.data.frame,logical(1))&!names(d)%in%c("raw","pcoa")])fl_av_csv(d[[k]],file.path(out,"tables",paste0(k,".csv")))
    fl_av_csv(d$clinical[d$clinical$sample_set=="review_core",],file.path(out,"tables/Main_clinical_crude_adjusted.csv"))
    fl_av_csv(fl_ap_clinical_wide(d$clinical),file.path(out,"tables/Main_clinical_comparison_wide.csv"))
    for(k in c("distance","shannon","composition"))fl_av_csv(d[[k]][d[[k]]$scenario=="core_adjusted",],file.path(out,"tables",paste0("Main_",k,".csv")))
    writeLines(enc2utf8(fl_av_summary(d)),file.path(out,"reports/results_summary_zh.md"),useBytes=TRUE)
    fl_ap_page(file.path(out,"reports/results_summary.html"),"脂肪肝确认设计：重点结果",c("主要与补充结果已分开；全部原始P。Nature合图见11模块。"),
      list("临床：粗与主要调整"=fl_ap_clinical_display(d$clinical),
        "主要相似性"=d$distance[d$distance$scenario=="core_adjusted",c("comparison","n_people","families","observed_within_mean","observed_between_mean","adjusted_distance_difference","p_raw_two_sided")],
        "主要Shannon"=d$shannon[d$shannon$scenario=="core_adjusted",c("generation","n","cases","families","estimate","CI_low","CI_high","p_raw")],
        "主要条件组成"=d$composition[d$composition$scenario=="core_adjusted",c("scheme","n","cases","families","R2_total_corrected","R2_partial","p_raw")]))
    fl_av_csv(input$manifest,file.path(out,"audit/input_manifest.csv"));file.copy(file.path(root,"config/approved_design_20261020.md"),file.path(out,"reports/approved_design_zh.md"))
    writeLines("COMPLETED",file.path(out,"run_status.txt"));final<-fl_output_complete_run(ctx)
    fl_av_landing(root,"10_Results_visualization",ctx$cfg$run_id,"reports/results_summary.html")
    cat("Completed result summary: ",final,"\n",sep="");invisible(final)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(e)})
}

fl_ap_figure_main<-function(args=character(),root){
  library(ggplot2);library(patchwork)
  cfg<-fl_config(root);cfg$run_ids['nature_composites']<-"nature_approved_design_v1";opt<-fl_entry_options(args,cfg,"nature_composites");cfg<-opt$cfg
  if(opt$mode=="--help"){cat("11: --run renders five main and four Supplementary composites from adopted approved results; --check checks inputs. R-only, no model refits.\n");return(invisible(NULL))}
  if(opt$mode=="--self-test")stop("Use the real-data layout/export QA; no simulated figures are generated")
  input<-fl_av_inputs(root);d<-input$data
  if(opt$mode=="--check"){cat("11 approved inputs checked; R backend\n");return(invisible(input$manifest))}
  oldpath<-Sys.getenv("PYTHONPATH",unset="");on.exit(Sys.setenv(PYTHONPATH=oldpath),add=TRUE)
  Sys.setenv(PYTHONPATH=paste(c(file.path(root,".figure_qa_python"),if(nzchar(oldpath))oldpath),collapse=.Platform$path.sep))
  oldtheme<-theme_set(nf_theme());on.exit(theme_set(oldtheme),add=TRUE)
  ctx<-fl_output_start_run(cfg);out<-ctx$run_dir;for(p in c("figures","source_data","private","audit","reports"))dir.create(file.path(out,p))
  writeLines("RUNNING",file.path(out,"run_status.txt"))
  tryCatch({
    bundle<-file.path(out,"audit/approved_figure_source_bundle.R")
    writeLines(c(readLines(file.path(root,"R/nature_export_approved.R"),encoding="UTF-8"),
      readLines(file.path(root,"R/approved_visualization.R"),encoding="UTF-8")),bundle,useBytes=TRUE)
    nf_qa_command(cfg,"validate_figure.py",c(shQuote(bundle),"--json"),file.path(out,"audit/source-preflight.json"))
    cl<-fl_av_main_clinical(d)
    specs<-list(list(name="Fig01_Cohort",plot=nf_baseline(d),width=183,height=150,panels=4),
      list(name="Fig02_Clinical",plot=cl$plot,width=183,height=82,panels=2),list(name="Fig03_Family_microbiome",plot=fl_av_family(d),width=183,height=82,panels=2),
      list(name="Fig04_Shannon",plot=fl_av_shannon(d),width=183,height=90,panels=2),list(name="Fig05_Composition",plot=fl_av_composition(d),width=183,height=95,panels=2),
      list(name="Supp01_Clinical_adult",plot=fl_av_supp_clinical(d),width=183,height=85,panels=3),
      list(name="Supp02_Similarity_sensitivity",plot=fl_av_supp_family(d),width=183,height=125,panels=3),
      list(name="Supp03_Shannon_sensitivity",plot=fl_av_supp_shannon(d),width=183,height=125,panels=3),
      list(name="Supp04_Composition_sensitivity",plot=fl_av_supp_composition(d),width=183,height=115,panels=2))
    fl_av_csv(input$manifest,file.path(out,"audit/input_manifest.csv"));fl_av_csv(cl$source,file.path(out,"source_data/clinical_prevalence_same_sample.csv"))
    for(i in seq_len(nrow(input$manifest))){r<-input$manifest[i,];if(!r$key%in%c("distribution","pcoa"))file.copy(r$path,file.path(out,"source_data",paste0(r$key,".csv")))}
    saveRDS(list(shannon=d$raw,pcoa=d$pcoa),file.path(out,"private/deidentified_plot_values.rds"))
    registry<-list();for(i in seq_along(specs)){cat("Approved Nature figure: ",specs[[i]]$name,"\n",sep="");flush.console();registry[[i]]<-nf_export(specs[[i]],out,cfg)}
    fl_av_csv(do.call(rbind,registry),file.path(out,"reports/figure_registry.csv"));fl_av_legends(d,out,specs)
    writeLines(enc2utf8(fl_av_summary(d)),file.path(out,"reports/results_summary_zh.md"),useBytes=TRUE)
    file.copy(file.path(root,"config/approved_design_20261020.md"),file.path(out,"reports/approved_design_zh.md"))
    writeLines(enc2utf8(c("# 图形QA记录","R-only；183mm；Arial；正文5–7pt，面板标签8pt；PDF/SVG/600dpi PNG。",
      "复用原项目图的风格及导出/对齐代码，依据新设计重建数据映射。无模拟、抽样、显著性筛选。",
      "每图最终对齐清单、原始及Cairo实际字号/碰撞审计见audit。所有WARN须结合成图审阅。",
      "Main仅呈现主要调整；补充模型明文命名；不可估不绘制假零效应。",
      "统计区间：临床与Shannon为家系sandwich逐项95%CI；相似性无独立配对CI；柱图为描述或R²无误差条。",
      "个体原始点与PCoA数值仅存private，公开Source Data为汇总统计，身份字段不进入图或公共表。")),file.path(out,"reports/QA_notes_zh.md"),useBytes=TRUE)
    if(any(unname(tools::md5sum(input$manifest$path))!=input$manifest$md5))stop("Source changed while rendering")
    saveRDS(lapply(specs,function(s)s$plot),file.path(out,"private/composite_plots.rds"))
    capture.output(sessionInfo(),file=file.path(out,"audit/sessionInfo.txt"));writeLines("COMPLETED",file.path(out,"run_status.txt"))
    final<-fl_output_complete_run(ctx);fl_av_landing(root,"11_Nature_composites",ctx$cfg$run_id)
    cat("Completed approved Nature figures: ",final,"\n",sep="");invisible(final)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(e)})
}
fl_av_legends<-function(d,out,specs){
  n<-sum(d$samples$n[d$samples$scenario=="core_adjusted"]);cases<-sum(d$samples$cases[d$samples$scenario=="core_adjusted"])
  legends<-c(
    Fig01_Cohort="Cohort and generation structure. Flow and analysis subsets; age median, IQR and range; observed steatosis proportions, including G2 minors and adults. Baseline: 1,275 people in 319 primary families. Age intervals are descriptive, not CIs. No inferential test. Family labels do not establish direct pedigrees or cohabitation.",
    Fig02_Clinical="Cross-generation clinical association. a, Observed prevalence on the same complete-case participants as the crude/main comparison. b, Main logistic GEE ORs with family-sandwich 95% CIs (G0-G1: 317 people/88 cases; G1-G2 and G0-G2: 484/17). Adjusted for younger age/sex, observed older count/mean age/female fraction, survey region and urbanicity; constants omitted. Crude ORs are in Source Data. All raw two-sided P; no multiplicity adjustment. Adult extensions are Supplementary.",
    Fig03_Family_microbiome="Familial microbiome similarity. a, Observed mean IF/BF Bray-Curtis distances on each main complete-case comparison. b, Adjusted IF minus BF distance; no independent-pair CI. Both-endpoint age ns3, sex, region, preparation batch, sequencing pool, input total, season, freezing interval and recent antibiotics; same-region/batch/pool indicators. Inference: 9,999 Freedman-Lane node residual permutations within generation roles and sex/5-year age blocks, absolute studentized t and add-one raw two-sided P. Pair counts are not independent sample sizes. Similarity does not establish strain transmission.",
    Fig04_Shannon=paste0("Steatosis and Shannon diversity. a, All individual observations in the main complete-case sample (n=",n,"; cases=",cases,"). Boxes: median/IQR, whiskers 1.5 IQR. b, Generation-specific steatosis minus no-steatosis means from Gaussian GEE; family-sandwich pointwise 95% CIs. Age ns3 by generation, sex, generation, region, preparation batch, sequencing pool, log10 matrix total, season cycle, freezing interval, stool-ultrasound interval and recent antibiotics. Raw two-sided P, no multiplicity correction. BMI/lifestyle are Supplementary conditional analyses."),
    Fig05_Composition=paste0("Community composition and conditional association. a, PCoA of Lingoes-corrected full-denominator Bray-Curtis distances for the main sample (n=",n,"). Coordinates are descriptive and not covariate-residualized. b, Steatosis explained total corrected inertia from partial dbRDA conditioning on family and the main adjustment set used for Shannon. Permutation schemes: within family, or within family and generation. Raw upper-tail F-test P from 9,999 reduced residual permutations; no multiplicity correction. Both schemes use the same fitted effect, so bar heights can coincide. Conditional within-family composition and marginal GEE Shannon are distinct estimands."),
    Supp01_Clinical_adult="Adult clinical extensions: core, lifestyle (alcohol, smoking, diet score, activity, energy) and additional BMI on a common complete-case sample within each comparison. Family GEE, pointwise robust 95% CIs, raw two-sided P. [NE] marks non-estimable logistic models; no fake zero/null points are drawn. Counts/diagnostic reasons and all missingness are in Source Data.",
    Supp02_Similarity_sensitivity="Familial similarity sensitivity analyses. Main, technical reference, adult nested lifestyle/BMI, bowel/preparation, DNA and no-recent-perturbation comparisons. Matched-sample baselines are identified in the labels and tables. Each dot is a distance coefficient without independent-pair CI; 9,999 constrained Freedman-Lane MRQAP permutations, raw two-sided P. Results from different samples cannot be interpreted solely as covariate attenuation.",
    Supp03_Shannon_sensitivity="Shannon sensitivity analyses. Generation-specific Gaussian GEE estimates and family-sandwich 95% CIs across the documented scenarios. Adult core/lifestyle/BMI share participants; bowel and DNA extensions have their own matched references. Exclusion of recent antibiotics/gastroenteritis changes the sample. Raw two-sided P in Source Data; no multiple-testing correction or mediation interpretation.",
    Supp04_Composition_sensitivity="Conditional composition sensitivity analyses. Steatosis explained total Lingoes-corrected inertia under the documented scenarios, with within-family and within-family-generation permutation inference. Identical bars across schemes reflect the same fitted model, not replicated samples. Source Data give raw upper-tail F-test P and effective exchangeability. No independent confidence interval is implied.")
  lines<-c("# Figure legends / 图注","","Main与Supplementary分开。完整协变量及中文解释见approved_design_zh.md和results_summary_zh.md。","")
  final_base<-sub("/_running/","/current/",gsub("\\","/",out,fixed=TRUE),fixed=TRUE)
  html<-c(paste0("<!doctype html><html lang='zh-CN'><meta charset='utf-8'><base href='file:///",final_base,"/'><title>脂肪肝确认设计合图</title><style>body{font:15px/1.6 Arial,'Microsoft YaHei';max-width:1200px;margin:25px auto}img{max-width:100%}p{max-width:1100px}</style>"),"<h1>脂肪肝：确认设计后的主要图与补充图</h1><p><a href='reports/results_summary_zh.md'>中文结果摘要</a> · <a href='reports/approved_design_zh.md'>完整设计</a> · <a href='source_data/clinical.csv'>全部临床模型</a></p>")
  for(s in specs){l<-legends[[s$name]];lines<-c(lines,paste0("## ",s$name),"",l,"");html<-c(html,paste0("<h2>",s$name,"</h2><a href='figures/",s$name,".pdf'>PDF</a> · <a href='figures/",s$name,".svg'>SVG</a><br><img src='figures/",s$name,".png'><p>",l,"</p>"))}
  writeLines(enc2utf8(lines),file.path(out,"reports/figure_legends.md"),useBytes=TRUE);writeLines(enc2utf8(c(html,"</html>")),file.path(out,"index.html"),useBytes=TRUE)
}
