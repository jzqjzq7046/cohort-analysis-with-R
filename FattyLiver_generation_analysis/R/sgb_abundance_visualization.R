# R-only figures from adopted effect estimates; no models are refitted.
fl_ab_viz_forest <- function(z,selected,groups,colours,title,x_label) {
  z<-z[z$feature_id%in%selected,,drop=FALSE]
  z$Group<-factor(z$Group,levels=groups)
  offsets<-if(length(groups)==1L)0 else seq(.26,-.26,length.out=length(groups))
  names(offsets)<-groups
  order_y<-rev(selected);z$y<-match(z$feature_id,order_y)+offsets[as.character(z$Group)]
  valid<-startsWith(z$status,"ESTIMATED")&is.finite(z$estimate)&is.finite(z$CI_low)&is.finite(z$CI_high)
  labels<-sub("^.*?(SGB[0-9]+).*$","\\1",order_y,perl=TRUE)
  p<-ggplot2::ggplot(z,ggplot2::aes(estimate,y,colour=Group))+
    ggplot2::geom_vline(xintercept=0,linetype=2,linewidth=.25,colour="#888888")+
    ggplot2::geom_errorbar(data=z[valid,,drop=FALSE],ggplot2::aes(xmin=CI_low,xmax=CI_high),
      orientation="y",width=0,linewidth=.3)+
    ggplot2::geom_point(data=z[valid,,drop=FALSE],size=1.55,shape=18)+
    ggplot2::scale_colour_manual(values=colours,drop=FALSE)+
    ggplot2::scale_y_continuous(breaks=seq_along(order_y),labels=labels,
      limits=c(.4,length(order_y)+.6),expand=ggplot2::expansion(mult=0))+
    ggplot2::labs(title=title,x=x_label,y=NULL)
  if(any(!valid))p<-p+ggplot2::geom_text(data=z[!valid,,drop=FALSE],ggplot2::aes(x=0,label="NE"),
    size=1.95,family="Arial",show.legend=FALSE)
  if(length(groups)==1L)p<-p+ggplot2::guides(colour="none")
  p
}

fl_ab_viz_effects <- function(r,selected) {
  cols<-c(G0="#F1B656",G1="#397FC7",G2="#040676")
  a<-r$results$primary$effects;a$Group<-"Overall"
  b<-r$results$generation$effects;b$Group<-sub("^case_","",b$term)
  p1<-fl_ab_viz_forest(a,selected,"Overall",c(Overall="#39877F"),"Overall association","Adjusted difference (SD)")
  p2<-fl_ab_viz_forest(b,selected,names(cols),cols,"Generation-specific associations","Adjusted difference (SD)")
  nf_patch(list(p1,p2),"AB",collect_guides=TRUE)+patchwork::plot_annotation(theme=nf_theme())
}

fl_ab_viz_extensions <- function(r,selected) {
  a<-do.call(rbind,lapply(c("primary","no_antibiotics"),function(v){z<-r$results[[v]]$effects;z$Group<-v;z}))
  a$Group<-ifelse(a$Group=="primary","Primary","No recent antibiotics")
  b<-do.call(rbind,lapply(c("adult_core","adult_lifestyle","adult_BMI"),function(v){z<-r$results[[v]]$effects;z$Group<-v;z}))
  map<-c(adult_core="Core",adult_lifestyle="Lifestyle",adult_BMI="Lifestyle + BMI")
  b$Group<-unname(map[b$Group])
  p1<-fl_ab_viz_forest(a,selected,c("Primary","No recent antibiotics"),
    c(Primary="#39877F","No recent antibiotics"="#B2556B"),"Antibiotic sensitivity","Adjusted difference (SD)")
  p2<-fl_ab_viz_forest(b,selected,unname(map),
    c(Core="#39877F",Lifestyle="#B2556B","Lifestyle + BMI"="#8370A6"),
    "Adult complete-case extensions","Adjusted difference (SD)")
  nf_patch(list(p1,p2),"AB",collect_guides=FALSE)
}

fl_ab_viz_composition <- function(r,selected) {
  a<-r$results$primary$effects;a$Group<-"Primary"
  b<-r$linda;b$estimate<-b$estimate_log2;b$Group<-"LinDA"
  p1<-fl_ab_viz_forest(a,selected,"Primary",c(Primary="#39877F"),"Primary relative abundance","Adjusted difference (SD)")
  p2<-fl_ab_viz_forest(b,selected,"LinDA",c(LinDA="#B2556B"),"Composition-adjusted sensitivity","Bias-corrected effect (log2)")
  nf_patch(list(p1,p2),"AB")
}

fl_ab_viz_run <- function(cfg) {
  module<-file.path(cfg$project_root,"outputs/18_SGB_abundance")
  pointer<-readLines(file.path(module,"CURRENT.txt"),n=1)
  if(!grepl("^current/[A-Za-z0-9_.-]+$",pointer))stop("Unsafe analysis pointer")
  source_run<-file.path(module,pointer)
  if(readLines(file.path(source_run,"run_status.txt"),n=1)!="COMPLETED")stop("Analysis is not adopted/completed")
  input<-file.path(source_run,"private/derived_results.rds");r<-readRDS(input)
  selected<-r$complete$feature_id[r$complete$shown_in_figure]
  selected<-selected[order(r$complete$primary_p_raw[match(selected,r$complete$feature_id)],selected)]
  ctx<-fl_output_start_run(cfg);cfg<-ctx$cfg;out<-ctx$run_dir
  for(v in c("figures","tables","reports","audit","private"))dir.create(file.path(out,v),showWarnings=FALSE)
  writeLines("RUNNING",file.path(out,"run_status.txt"))
  oldpath<-Sys.getenv("PYTHONPATH",unset="");on.exit(Sys.setenv(PYTHONPATH=oldpath),add=TRUE)
  Sys.setenv(PYTHONPATH=paste(c(file.path(cfg$project_root,".figure_qa_python"),if(nzchar(oldpath))oldpath),collapse=.Platform$path.sep))
  tryCatch({
    specs<-list(
      list(name="Fig06_SGB_abundance",plot=fl_ab_viz_effects(r,selected),width=183,height=132,panels=2,export_tiff=TRUE),
      list(name="Supp05_SGB_abundance_extensions",plot=fl_ab_viz_extensions(r,selected),width=183,height=140,panels=2,export_tiff=TRUE),
      list(name="Supp06_SGB_composition_sensitivity",plot=fl_ab_viz_composition(r,selected),width=183,height=132,panels=2,export_tiff=TRUE))
    bundle<-file.path(out,"audit/figure_source_bundle.R")
    writeLines(c(readLines(file.path(cfg$project_root,"R/nature_export_approved.R")),
      readLines(file.path(cfg$project_root,"R/sgb_abundance_visualization.R"))),bundle)
    nf_qa_command(cfg,"validate_figure.py",c(shQuote(bundle),"--json"),file.path(out,"audit/source_preflight.json"))
    exports<-do.call(rbind,lapply(specs,function(z)nf_export(z,out,cfg)))
    fl_write_csv(exports,file.path(out,"audit/export_manifest.csv"))
    fl_write_csv(r$complete[r$complete$feature_id%in%selected,],file.path(out,"tables/SGB_figure_source_EN.csv"))
    display<-do.call(rbind,lapply(names(r$results),function(v){
      z<-r$results[[v]]$effects;z<-z[z$feature_id%in%selected,,drop=FALSE]
      data.frame(SGB_id=sub("^t__","",z$feature_id),analysis=v,
        group=ifelse(grepl("^case_",z$term),sub("^case_","",z$term),"Overall"),
        estimate=z$estimate,CI_low=z$CI_low,CI_high=z$CI_high,P_value=z$p_raw,
        status=z$status,unit="SD of core-reference log10 relative abundance")
    }))
    z<-r$linda;z<-z[z$feature_id%in%selected,,drop=FALSE]
    display<-rbind(display,data.frame(SGB_id=sub("^t__","",z$feature_id),analysis="LinDA",
      group="Overall",estimate=z$estimate_log2,CI_low=z$CI_low,CI_high=z$CI_high,
      P_value=z$p_raw,status=z$status,unit="LinDA bias-corrected log2 effect"))
    fl_write_csv(display,file.path(out,"tables/SGB_display_effects_EN.csv"))
    fl_write_csv(data.frame(path=input,md5=unname(tools::md5sum(input))),file.path(out,"audit/source_analysis_manifest.csv"))
    saveRDS(lapply(specs,function(z)z$plot),file.path(out,"private/figure_plots.rds"))
    legends<-c("# Figure legends", "",
      "Fig. 6 | Hepatic steatosis-associated SGB abundance across generations.",
      "a, Overall adjusted associations. b, Generation-specific adjusted associations. Points show steatosis minus no-steatosis differences; lines show 95% Satterthwaite confidence intervals. All models use family random intercepts and the prespecified core adjustment set. Outcomes are log10 relative abundance with the fixed pseudocount, standardized once in the 1267-person core sample; all panels/generations use the same feature-specific scale. G0/G1/G2 colours are gold/blue/navy. The 12 displayed SGBs are selected by the prespecified primary raw-P and detection-support rule; selection is for display and does not establish multiplicity-controlled discoveries. The complete table includes all 695 SGBs. G2 has 17 cases; NE denotes not estimable.","",
      "Supplementary Fig. 5 | Sensitivity to recent antibiotics and adult covariate adjustment.",
      "a, Primary sample (1267 people, 232 cases) and exclusion of recent antibiotics (1165 people, 210 cases). b, Core, lifestyle and lifestyle-plus-BMI models fitted to the identical 736 adult complete cases (158 cases). Points and lines indicate standardized differences and 95% confidence intervals. Adult results are complete-case sensitivities, not pooled clinical multiple-imputation results. The fixed full-core standardization is retained; population changes between a and b must not be attributed solely to adjustment.","",
      "Supplementary Fig. 6 | Sensitivity to compositional modelling.",
      "a, Primary standardized log10 relative-abundance differences. b, LinDA bias-corrected log2 effects and 95% intervals in the same core sample, with family random intercepts and the same nuisance column space. LinDA uses all 695 proportions jointly, its official zero handling and default 3% winsorization. Coefficients have different units and null interpretations; panels assess sensitivity rather than numerical equivalence or external validation. Raw P values are exploratory.")
    writeLines(legends,file.path(out,"reports/Figure_legends_EN.md"))
    writeLines(c("# Figure QA record","","Style-only inheritance from the adopted project figures: Arial, white background, 183-mm width, restrained palette, generation-specific gold/blue/navy and separate teal/berry comparison colours.",
      "Both panels in each figure display the identical prespecified SGB list. All finite estimable effects and their full 95% intervals are retained. Different x-axis units in the compositional sensitivity are explicitly labelled.",
      "Final rendered alignment, PDF text and collision reports are in audit/. Manual panel-by-panel inspection is recorded separately after rendering."),file.path(out,"audit/visual_QA_notes.md"))
    capture.output(sessionInfo(),file=file.path(out,"audit/sessionInfo.txt"))
    writeLines("COMPLETED",file.path(out,"run_status.txt"));final<-fl_output_complete_run(ctx)
    cat("COMPLETED FIGURES: ",final,"\n",sep="");invisible(final)
  },error=function(e){fl_output_fail_run(ctx,conditionMessage(e));stop(e)})
}
