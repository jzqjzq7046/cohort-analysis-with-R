# Joint technical control now belongs to the main models; no repeated fits.
fl_joint_technical_verification_main<-function(args=character(),root){
 cfg<-fl_config(root);cfg$run_ids['technical_covariates']<-"technical_approved_design_verification_v1"
 opt<-fl_entry_options(args,cfg,"technical_covariates");cfg<-opt$cfg
 if(opt$mode=="--help"){cat("09 verifies joint technical control in adopted 07/08; no model refits. Use --run or --check.\n");return(invisible(NULL))}
 if(opt$mode=="--self-test")stop("Use tests/validate_approved_results.R")
 input<-fl_av_inputs(root);d<-input$data
 checks<-data.frame(module=c("07 MRQAP","08 Shannon","08 partial dbRDA"),
   sequencing_pool_in_every_scenario=c(all(grepl("seq_pool",d$distance$covariates,fixed=TRUE)),
     all(grepl("seq_pool",d$shannon$covariates_added,fixed=TRUE)),all(grepl("seq_pool",d$composition$covariates_added,fixed=TRUE))))
 if(!all(checks$sequencing_pool_in_every_scenario))stop("Joint technical control missing")
 if(opt$mode=="--check"){print(checks);return(invisible(checks))}
 ctx<-fl_output_start_run(cfg);out<-ctx$run_dir;for(p in c("technical","audit"))dir.create(file.path(out,p))
 writeLines("RUNNING",file.path(out,"run_status.txt"))
 tryCatch({
  fl_av_csv(checks,file.path(out,"technical/joint_technical_verification.csv"));fl_av_csv(input$manifest,file.path(out,"audit/input_manifest.csv"))
  writeLines(enc2utf8(c("# 联合技术控制已纳入主模型","09为只读核查，不重复拟合旧药物/测序池扩展模型。",
    "所有07/08方案包含建库批次与测序池；共线项按QR处理，并记录别名。结果直接阅读07/08主要和补充表。")),file.path(out,"technical/README_zh.md"),useBytes=TRUE)
  writeLines("COMPLETED",file.path(out,"run_status.txt"));fl_output_complete_run(ctx)
 },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(e)})
}
