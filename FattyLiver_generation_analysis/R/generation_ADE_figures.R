# R-only complete, paginated displays. No re-fitting or q-value recomputation.
ade_theme <- function() tf_plot_theme()+ggplot2::theme(
  axis.text=ggplot2::element_text(size=5.7,colour='black'),
  axis.text.y=ggplot2::element_text(size=5.4),
  strip.text=ggplot2::element_text(size=6),strip.background=ggplot2::element_blank(),
  panel.spacing.x=grid::unit(4,'mm'),legend.title=ggplot2::element_text(size=6.3),
  legend.text=ggplot2::element_text(size=5.7))
ade_model_label <- function(x) unname(c(B1='Model 3: all ages',Adult_B1_same_sample='Model 3: adults',B2='Model 4: adults')[x])
ade_model_factor <- function(x) factor(ade_model_label(x),levels=ade_model_label(c('B1','Adult_B1_same_sample','B2')))
ade_order <- function(t) t$feature_id[order(t$P_BH,t$P_raw,t$feature_id,na.last=TRUE)]
ade_ylabels <- function(anno,ids) {
  sapply(ids,function(id){
    j<-match(id,anno$feature_id)
    if('scientific_name'%in%names(anno)&&!is.na(j)&&!is.na(anno$scientific_name[j]))
      paste0('italic(',encodeString(anno$scientific_name[j],quote='"'),')~',encodeString(paste0('(',sub('^t__','',id),')'),quote='"'))else
        encodeString(sub('^t__','',id),quote='"')
  },USE.NAMES=FALSE) |> parse(text=_)
}
ade_ratio_breaks <- function(lo,hi) {
  if(log2(hi/lo)<2){
    b<-pretty(c(lo,hi),n=5)
    return(sort(unique(c(1,b[b>0&b>=lo&b<=hi]))))
  }
  powers<-seq(floor(log2(lo)),ceiling(log2(hi)))
  if(length(powers)>7)powers<-powers[unique(round(seq(1,length(powers),length.out=7)))]
  sort(unique(c(1,2^powers)))
}
ade_ratio_text <- function(v) vapply(v,function(z){
  if(!is.finite(z))''else if(z>=.01&&z<1000)format(signif(z,3),scientific=FALSE,trim=TRUE)else format(z,digits=2,scientific=TRUE,trim=TRUE)
},character(1))
ade_export <- function(p,stem,out,cfg,width=183,height=150) {
  ov_export(p,stem,out,cfg,width=width,height=height)
}
ade_overview <- function(r,cfg,out) {
  if(r$kind=='A'){
    z<-r$heterogeneity[r$heterogeneity$model=='B1'&tf_bool(r$heterogeneity$analysis_valid),,drop=FALSE]
    a<-r$effects[r$effects$model=='B1',,drop=FALSE]
    z$spread<-vapply(z$feature_id,function(id)diff(range(a$beta_log2[a$feature_id==id],finite=TRUE)),numeric(1))
    z$x<-2^z$spread;z$y<--log10(pmax(z$P_BH,.Machine$double.xmin))
    z$state<-factor(ifelse(z$FDR_lt_005,'BH q < 0.05','Other tested SGBs'),levels=c('Other tested SGBs','BH q < 0.05'))
    p<-ggplot2::ggplot(z,ggplot2::aes(x,y,colour=state))+
      ggplot2::geom_hline(yintercept=-log10(.05),linetype=2,colour='#888888',linewidth=.25)+
      ggplot2::geom_point(size=1.35,alpha=.85,shape=16)+
      ggplot2::scale_colour_manual(values=setNames(unname(cfg$binary_colors),levels(z$state)))+
      ggplot2::scale_x_log10(breaks=ade_ratio_breaks(min(z$x),max(z$x)),labels=ade_ratio_text)+
      ggplot2::labs(title=sprintf('Generation heterogeneity: %d of %d valid SGBs at BH q < 0.05',sum(z$FDR_lt_005),nrow(z)),
        x='Largest / smallest generation abundance ratio',y='−log10(BH q)')+ade_theme()
    stem<-'Figure41_All_SGB_generation_heterogeneity'
  }else{
    z<-r$effects[r$effects$model=='B1'&r$effects$generation=='Overall'&tf_bool(r$effects$analysis_valid),,drop=FALSE]
    z$x<-z$abundance_ratio;z$y<--log10(pmax(z$P_BH,.Machine$double.xmin))
    z$state<-factor(ifelse(!z$FDR_lt_005,'Other tested pathways',ifelse(z$beta_log2<0,'Lower with steatosis','Higher with steatosis')),
      levels=c('Other tested pathways','Lower with steatosis','Higher with steatosis'))
    p<-ggplot2::ggplot(z,ggplot2::aes(x,y,colour=state))+
      ggplot2::geom_vline(xintercept=1,linetype=2,colour='#888888',linewidth=.25)+
      ggplot2::geom_hline(yintercept=-log10(.05),linetype=2,colour='#888888',linewidth=.25)+
      ggplot2::geom_point(size=1.35,alpha=.85,shape=16)+
      ggplot2::scale_colour_manual(values=c('Other tested pathways'='#B8B8B8','Lower with steatosis'=unname(cfg$binary_colors[1]),'Higher with steatosis'=unname(cfg$binary_colors[2])))+
      ggplot2::scale_x_log10(breaks=ade_ratio_breaks(min(z$x),max(z$x)),labels=ade_ratio_text)+
      ggplot2::labs(title=sprintf('Steatosis-associated pathways: %d of %d valid tests at BH q < 0.05',sum(z$FDR_lt_005),nrow(z)),
        x='Adjusted abundance ratio',y='−log10(BH q)')+ade_theme()
    stem<-'Figure42_All_pathway_associations'
  }
  ade_table(z,paste0(stem,'_source_EN'),out)
  ade_export(p,stem,out,cfg,width=160,height=115)
}

ade_heatmaps <- function(r,ids,cfg,out,prefix,models=c('B1','Adult_B1_same_sample','B2'),generation='all') {
  if(!length(ids))return(data.frame())
  e<-r$effects[r$effects$feature_id%in%ids&r$effects$model%in%models,,drop=FALSE]
  e<-e[if(generation=='Overall')e$generation=='Overall'else e$generation!='Overall',,drop=FALSE]
  valid<-tf_bool(e$analysis_valid)&is.finite(e$beta_log2)
  # One common linear log2 colour range for the complete display family.
  # Round outward, retaining small pathway differences without saturation.
  lim<-max(.25,ceiling(max(abs(e$beta_log2[valid]),na.rm=TRUE)*4)/4)
  pages<-split(ids,ceiling(seq_along(ids)/32L));exports<-rows<-list()
  for(i in seq_along(pages)){
    keep<-pages[[i]]
    grid<-expand.grid(feature_id=keep,generation=if(generation=='Overall')'Overall'else c('G0','G1','G2'),model=models,stringsAsFactors=FALSE)
    z<-merge(grid,e,by=c('feature_id','generation','model'),all.x=TRUE,sort=FALSE)
    z$feature_id<-factor(z$feature_id,levels=rev(keep));z$model_label<-factor(ade_model_label(z$model),levels=ade_model_label(models))
    z$generation<-factor(z$generation,levels=if(generation=='Overall')'Overall'else c('G0','G1','G2'))
    z$fill<-ifelse(tf_bool(z$analysis_valid),z$beta_log2,NA_real_)
    z$star<-ifelse(tf_bool(z$FDR_lt_005),'*','')
    p<-ggplot2::ggplot(z,ggplot2::aes(generation,feature_id,fill=fill))+
      ggplot2::geom_tile(width=.96,height=.94,colour='white',linewidth=.12)+
      ggplot2::geom_text(ggplot2::aes(label=star),size=5.7/(72.27/25.4),colour='black',na.rm=TRUE)+
      ggplot2::facet_grid(.~model_label,drop=FALSE)+
      ggplot2::scale_fill_gradient2(low=cfg$binary_colors[1],mid='white',high=cfg$binary_colors[2],midpoint=0,
        limits=c(-lim,lim),na.value='#E4E4E4',breaks=c(-lim,0,lim),labels=ade_ratio_text(2^c(-lim,0,lim)),
        name='Adjusted abundance ratio',guide=ggplot2::guide_colourbar(direction='horizontal',title.position='top',
          display='rectangles',nbin=64,ticks=FALSE,barwidth=grid::unit(52,'mm'),barheight=grid::unit(2.5,'mm')))+
      ggplot2::scale_y_discrete(labels=function(v)ade_ylabels(r$annotation,v),drop=FALSE,
        expand=if(generation=='Overall')ggplot2::expansion(add=.6)else ggplot2::expansion(add=c(.6,1.3)))+
      ggplot2::labs(title=sprintf('%s: page %d/%d',if(r$kind=='A')'All heterogeneity-supported SGBs'else if(generation=='Overall')'All steatosis-associated pathways'else'All generation-heterogeneous pathways',i,length(pages)),x=NULL,y=NULL)+
      ade_theme()+ggplot2::theme(axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),legend.position='bottom')
    if(generation!='Overall'){
      band<-expand.grid(generation=c('G0','G1','G2'),model_label=levels(z$model_label),stringsAsFactors=FALSE)
      band$generation<-factor(band$generation,levels=c('G0','G1','G2'))
      band$model_label<-factor(band$model_label,levels=levels(z$model_label))
      p<-p+ggplot2::geom_segment(data=band,ggplot2::aes(x=as.numeric(generation)-.32,xend=as.numeric(generation)+.32,
        y=length(keep)+.8,yend=length(keep)+.8,colour=generation),inherit.aes=FALSE,linewidth=1.1,show.legend=FALSE)+
        ggplot2::scale_colour_manual(values=cfg$generation_colors)
    }
    stem<-sprintf('%s_%02d',prefix,i)
    exports[[i]]<-ade_export(p,stem,out,cfg,width=183,height=max(55,40+length(keep)*3.6))
    rows[[i]]<-data.frame(figure=stem,page=i,display_order=seq_along(keep)+(i-1)*32,feature_id=keep,
      selection='All primary BH-supported features; no count cap',common_log2_colour_limit=lim)
  }
  ade_table(ov_bind(rows),paste0(prefix,'_page_manifest_EN'),out)
  ov_bind(exports)
}

ade_forests <- function(r,ids,cfg,out,prefix,model='B1',generation='all') {
  if(!length(ids))return(data.frame())
  e<-r$effects[r$effects$feature_id%in%ids&r$effects$model==model,,drop=FALSE]
  e<-e[if(generation=='Overall')e$generation=='Overall'else e$generation!='Overall',,drop=FALSE]
  e<-e[tf_bool(e$analysis_valid)&is.finite(e$ratio_CI_low)&is.finite(e$ratio_CI_high)&e$ratio_CI_low>0,,drop=FALSE]
  stopifnot(nrow(e)>0)
  rr<-range(c(e$ratio_CI_low,e$ratio_CI_high));ticks<-ade_ratio_breaks(rr[1],rr[2]);limits<-2^(log2(rr)+c(-.12,.12))
  pages<-split(ids,ceiling(seq_along(ids)/28L));exports<-rows<-list()
  for(i in seq_along(pages)){
    keep<-pages[[i]];z<-e[e$feature_id%in%keep,,drop=FALSE]
    z$feature_id<-factor(z$feature_id,levels=rev(keep));z$generation<-factor(z$generation,levels=if(generation=='Overall')'Overall'else c('G0','G1','G2'))
    z$direction<-factor(ifelse(z$beta_log2<0,'Lower with steatosis','Higher with steatosis'),levels=c('Lower with steatosis','Higher with steatosis'))
    z$point_fill<-if(generation=='Overall')z$direction else z$generation
    p<-ggplot2::ggplot(z,ggplot2::aes(abundance_ratio,feature_id))+
      ggplot2::geom_vline(xintercept=1,linetype=2,colour='#888888',linewidth=.25)+
      ggplot2::geom_segment(ggplot2::aes(x=ratio_CI_low,xend=ratio_CI_high,yend=feature_id),colour='black',linewidth=.3)+
      ggplot2::geom_point(ggplot2::aes(fill=point_fill),shape=21,size=1.6,colour='black',stroke=.2)+
      ggplot2::scale_fill_manual(values=if(generation=='Overall')setNames(unname(cfg$binary_colors),levels(z$direction))else cfg$generation_colors)+
      ggplot2::scale_x_log10(limits=limits,breaks=ticks,labels=ade_ratio_text)+
      ggplot2::scale_y_discrete(limits=rev(keep),labels=function(v)ade_ylabels(r$annotation,v),drop=FALSE)+
      ggplot2::labs(title=sprintf('%s: page %d/%d',if(r$kind=='A')'Generation-specific steatosis associations'else if(generation=='Overall')'Steatosis-associated pathways'else'Generation-specific pathway associations',i,length(pages)),
        x='Adjusted abundance ratio (95% CI)',y=NULL)+ade_theme()+ggplot2::theme(legend.position='none')
    if(generation!='Overall')p<-p+ggplot2::facet_grid(.~generation,drop=FALSE)
    stem<-sprintf('%s_%02d',prefix,i)
    exports[[i]]<-ade_export(p,stem,out,cfg,width=183,height=max(55,32+length(keep)*4))
    rows[[i]]<-data.frame(figure=stem,page=i,display_order=seq_along(keep)+(i-1)*28,feature_id=keep,
      selection='All primary BH-supported features; no count cap',common_ratio_lower=limits[1],common_ratio_upper=limits[2])
  }
  ade_table(ov_bind(rows),paste0(prefix,'_page_manifest_EN'),out)
  ov_bind(exports)
}

ade_searchable_table <- function(r,out) {
  esc<-function(s){s<-gsub('&','&amp;',as.character(s),fixed=TRUE);s<-gsub('<','&lt;',s,fixed=TRUE);gsub('>','&gt;',s,fixed=TRUE)}
  t<-r$effects
  h<-r$heterogeneity
  if(nrow(h)){
    k<-match(paste(t$model,t$feature_id),paste(h$model,h$feature_id))
    t$heterogeneity_P<-h$P_raw[k];t$heterogeneity_q<-h$P_BH[k]
  }
  fields<-intersect(c('model','feature_id','feature_label','pathway_name','generation','N','cases','abundance_ratio',
    'ratio_CI_low','ratio_CI_high','P_raw','P_BH','heterogeneity_P','heterogeneity_q','status'),names(t))
  t<-t[,fields,drop=FALSE]
  for(v in names(t))if(is.numeric(t[[v]]))t[[v]]<-ifelse(is.na(t[[v]]),'',format(t[[v]],digits=5,trim=TRUE))
  rows<-apply(t,1,function(z)paste0('<tr>',paste0('<td>',esc(z),'</td>',collapse=''),'</tr>'))
  title<-if(r$kind=='A')'All SGB generation associations'else'All functional pathway associations'
  html<-paste0('<!doctype html><meta charset="utf-8"><title>',title,'</title>',
    '<style>body{font:14px Arial;margin:24px;color:#171717}table{border-collapse:collapse;font-size:12px}th,td{padding:7px;border-bottom:1px solid #ddd;white-space:nowrap}th{position:sticky;top:0;background:#eef4f8}input{width:500px;padding:10px}tr:hover{background:#fff2e5}</style>',
    '<h1>',title,'</h1><p>All attempted results are available in the accompanying CSV files. Search by feature, model, generation or status. Invalid estimates remain flagged. No top-feature display cap.</p>',
    '<input id="query" placeholder="Search SGB, species, pathway, generation or model" oninput="filterRows()">',
    '<p id="count"></p><table><thead><tr>',paste0('<th>',esc(fields),'</th>',collapse=''),'</tr></thead><tbody>',paste(rows,collapse='\n'),'</tbody></table>',
    '<script>function filterRows(){let q=document.getElementById("query").value.toLowerCase(),n=0;document.querySelectorAll("tbody tr").forEach(r=>{let ok=r.textContent.toLowerCase().includes(q);r.style.display=ok?"":"none";n+=ok});document.getElementById("count").textContent=n+" rows"}filterRows()</script>')
  writeLines(html,file.path(out,'tables/All_results_searchable_EN.html'),useBytes=TRUE)
}

ade_reports <- function(r,cfg,out) {
  for(nm in c('effects','heterogeneity')){
    z<-r[[nm]];if(nrow(z))ade_table(z,if(nm=='effects')'All_feature_effects_EN'else'All_joint_generation_tests_EN',out)
  }
  ade_table(r$annotation,'Feature_annotations_EN',out)
  diagnostics<-ov_bind(lapply(r$fits,function(f){z<-f$diagnostics;z$mode<-f$mode;z}))
  ade_table(diagnostics,'All_model_diagnostics_EN',out)
  ade_table(r$counts,'Analysis_sample_counts_EN',out);ade_table(r$summary,'Analysis_model_summary_EN',out)
  ade_searchable_table(r,out)
  methods<-c('# Methods and figure legends','',
    'Outcome: observed ultrasound hepatic steatosis, not MASLD. The adopted two-generation primary-label cohort is retained. Primary labels do not establish direct ancestry or current cohabitation.',
    'Model 3 includes generation, generation-specific age natural splines (3 df per generation), sex, library batch, sequencing pool, seasonal sine/cosine and recent antibiotic use, with a primary-family random intercept. Model 4 additionally includes smoking, alcohol, diet score, physical activity, log energy and continuous BMI; adult Model 3/4 have identical participants.',
    'Official MaAsLin2 1.24.1, normalization NONE after fixed-denominator relative abundance, transformation LOG, standardize FALSE. Zero values are replaced by half the smallest positive feature value in the actual model sample, then log2-transformed. Ratios refer to the zero-replaced relative-abundance scale, not absolute microbial load.',
    'Generation-specific effects are estimated in one model with three generation-specific disease indicators and the unchanged nuisance design. The two independent contrasts G0-G1 and G1-G2 are tested jointly with lmerTest::contestMD (Satterthwaite two-df F-test). The models have the same mean space as a disease-by-generation interaction. Individual coefficients and raw P are verified against the official package output. Singular fits remain flagged; nonconverged fits are not discoveries.',
    'Multiplicity: separate BH families for the joint generation tests and the three-generation feature effects within each model. Overall pathway disease associations form their own families. Attempted but invalid fits contribute P=1 to the correction denominator, with public q left missing. All primary thresholds are q<0.05; no selected-subset BH is used.',
    if(r$kind=='A')'Figure 41a is a volcano-style heterogeneity overview of every valid primary joint test. Its x coordinate is 2^(max generation beta - min generation beta), the largest divided by the smallest generation-specific disease abundance ratio; its y coordinate is -log10(joint BH q). The nonnegative x coordinate has no enrichment/depletion direction. The three-generation joint test is not a signed two-group comparison. Blue points have q>=0.05 and orange points q<0.05. All ten SGB IDs displayed in panel b are labelled at their corresponding panel-a points with collision-avoiding leader lines. X limits are fixed at 1 to 32 with ticks 1,2,4,8,16,32; the horizontal dashed line marks q=0.05. All original q values are unchanged.'else
      'Overview: every valid attempted primary overall pathway test is shown, including nonsignificant pathways; the axes show the adjusted disease abundance ratio and -log10(BH q).',
    if(r$kind=='A')'Figure 41b displays exactly the ten valid all-age Model-3 SGBs with the smallest raw joint-interaction P values, with ties resolved by SGB ID, regardless of whether q<0.05. Rows are a ranked exploratory display, not ten FDR discoveries. Stars beside row labels mark the original joint BH q<0.05. Points and bars show the three generation-specific adjusted disease abundance ratios and 95% CIs; G0/G1/G2 retain their fixed colours. The log axis is fixed at 0.01 to 10 with ticks 0.01,0.1,1,10, and a dashed null at 1. No interval is truncated. All 587 primary joint results and the adult supplements remain in the complete tables.'else
      'Detailed heatmaps/forests include every primary BH-supported feature, ordered by q, raw P and ID, with pagination only for legibility. Heatmaps use common symmetric log2 colour limits within a display family; keys show ratios. Grey cells denote filtered/unavailable estimates. Asterisks mark the secondary feature-effect BH q<0.05 and are distinct from row selection by primary heterogeneity q. Forest dots and bars show ratios and 95% intervals; physical axes are audited and all pages share limits.',
    'All unique, unambiguous species binomials are italicized and followed by SGB IDs; unresolved/conflicting names retain IDs. Full original taxonomy and pathway names are in the annotation table. All attempted numerical results are retained independently of display selection.',
    'This is an exploratory same-cohort extension, not independent validation. Generation and chronological age are closely related; G2 has relatively few cases. No direction of transmission or causal mediation is inferred.',
    'References: Wang et al. Nature Medicine 2024;30:1722-1731, doi:10.1038/s41591-024-03038-y (age/metabolic interaction; our SGB interaction is an adaptation). Mallick et al. PLOS Computational Biology 2021;17:e1009442, doi:10.1371/journal.pcbi.1009442. Kuznetsova et al. Journal of Statistical Software 2017;82(13), doi:10.18637/jss.v082.i13 (contrast tests).')
  if(r$kind=='D')methods<-c(methods,
    'Functional input is observed HUMAnN pathway abundance, not a prediction from SGB identities. Only unstratified community pathways are tested. Normalization uses all supplied community rows before filtering, including special rows when present, following the HUMAnN community relative-abundance convention. The current export supplies 620 biological community rows and no UNMAPPED/UNINTEGRATED/UNGROUPED rows: the denominator therefore comprises all reported biological community pathways, not all raw reads or genes. Original abundance units are not established by the assay name counts; the recorded source sums are variable and local relative-abundance normalization is explicit. Filter: >=10% prevalence and nonzero variance. Community totals and species-stratified pathway quantities need not be additive. Functional potential is not transcription or metabolite production.',
    'Functional references: Mei et al. Nature Medicine 2024;30:2265-2276, doi:10.1038/s41591-024-03067-7, Fig. 3; Kim et al. Nature Metabolism 2025;7:1476-1492, doi:10.1038/s42255-025-01318-6, Fig. 3. These source designs are adapted to ultrasound steatosis, primary-family dependence and the current adjustment sets; no claim of an exact multi-cohort/multi-omic replication is made.',
    'Specific differences from Mei et al.: the present single-cohort extension uses official MaAsLin2 family-intercept models rather than cross-cohort MMUPHin meta-analysis. It retains the project filter of nonzero detection in >=10% of samples, rather than the source minimum abundance of 0.0002, correlated-feature representative selection and dominance filtering. All eligible supplied pathways are retained. Source q<0.10 is replaced by the predeclared project q<0.05. No strain-specific pathway carriage model is inferred from these community-level data.')
  writeLines(methods,file.path(out,'reports/Methods_and_figure_legends_EN.md'))
  lines<-c('# 本次分析结果','',capture.output(print(r$summary,row.names=FALSE)),
    '',if(r$kind=='A')'主要结论依据全龄Model 3的两自由度代际异质性BH检验；每代显著性不能替代交互检验。'else
      '总体通路关联与代际交互分属不同检验族；结果代表宏基因组功能潜力。',
    '全部结果、失败诊断、注释和检出率过滤记录均保留。',
    if(r$kind=='A')'主图改为全部受检SGB的交互总览＋交互原始P值排名前10的三代效应森林图；未通过FDR者也展示。星号标记正式BH q<0.05发现，固定横轴不随数据自动变化。'else
      '显著多少展示多少，分页不改变筛选范围。')
  primary<-if(r$kind=='A')r$heterogeneity[r$heterogeneity$model=='B1',,drop=FALSE]else
    r$effects[r$effects$model=='B1'&r$effects$generation=='Overall',,drop=FALSE]
  discoveries<-primary[primary$FDR_lt_005,,drop=FALSE]
  sensitivity<-data.frame(feature_id=discoveries$feature_id,primary_P=discoveries$P_raw,primary_q=discoveries$P_BH)
  for(m in c('Adult_B1_same_sample','B2')){
    t<-if(r$kind=='A')r$heterogeneity else r$effects[r$effects$generation=='Overall',,drop=FALSE]
    t<-t[t$model==m,,drop=FALSE];j<-match(sensitivity$feature_id,t$feature_id)
    sensitivity[[paste0(m,'_valid')]]<-tf_bool(t$analysis_valid[j])
    sensitivity[[paste0(m,'_q')]]<-t$P_BH[j]
    sensitivity[[paste0(m,'_BH_supported')]]<-tf_bool(t$FDR_lt_005[j])
    if(r$kind=='D')sensitivity[[paste0(m,'_same_direction')]]<-sign(discoveries$beta_log2)==sign(t$beta_log2[j])
  }
  ade_table(sensitivity,'Primary_discovery_sensitivity_summary_EN',out)
  lines<-c(lines,'',sprintf('全龄主分析：%d个合格受检%s，%d个有效检验，%d个通过BH q<0.05。',
    nrow(primary),if(r$kind=='A')'SGB代际交互'else'功能通路',sum(tf_bool(primary$analysis_valid)),nrow(discoveries)),
    sprintf('上述主发现中，成人基础模型支持%d个，加入生活方式和BMI的成人模型支持%d个。',
      sum(sensitivity$Adult_B1_same_sample_BH_supported),sum(sensitivity$B2_BH_supported)),
    if(r$kind=='A')'这里的“支持”指正式代际交互检验通过相应模型自己的完整BH检验族。它不等于某一代单独P<0.05。'else
      '这里的“支持”指通路总体关联通过相应成人模型自己的完整BH检验族。')
  writeLines(lines,file.path(out,'reports/Results_summary_zh.md'))
}

ade_figures <- function(r,cfg,out) {
  if(r$kind=='A'){
    sys.source(file.path(cfg$project_root,'R/generation_A_top10_figure.R'),envir=.GlobalEnv)
    discoveries<-r$heterogeneity[r$heterogeneity$model=='B1'&r$heterogeneity$FDR_lt_005,,drop=FALSE]
    discoveries<-discoveries[match(ade_order(discoveries),discoveries$feature_id),,drop=FALSE]
    ade_table(discoveries,'Primary_heterogeneity_discoveries_EN',out)
    # Old display manifests belong to archived figures, not the new current set.
    obsolete<-list.files(file.path(out,'tables'),
      '^Figure41_All_(heterogeneous_SGB_(forest|heatmap)_page_manifest|SGB_generation_heterogeneity_source)_EN\\.(csv|html)$',full.names=TRUE)
    if(length(obsolete)){
      stopifnot(all(normalizePath(dirname(obsolete),winslash='/')==normalizePath(file.path(out,'tables'),winslash='/')))
      if(!all(file.remove(obsolete)))stop('Could not retire old A display manifests from the new run')
    }
    exports<-list(ade_A_top10_combined(r,cfg,out))
  }else{
    exports<-list(ade_overview(r,cfg,out))
    t<-r$effects[r$effects$model=='B1'&r$effects$generation=='Overall'&r$effects$FDR_lt_005,,drop=FALSE];ids<-ade_order(t)
    ade_table(t[match(ids,t$feature_id),,drop=FALSE],'Primary_pathway_discoveries_EN',out)
    exports[[2]]<-ade_heatmaps(r,ids,cfg,out,'Figure42_All_significant_pathway_heatmap',generation='Overall')
    exports[[3]]<-ade_forests(r,ids,cfg,out,'Figure42_All_significant_pathway_forest',generation='Overall')
    h<-r$heterogeneity[r$heterogeneity$model=='B1'&r$heterogeneity$FDR_lt_005,,drop=FALSE];hi<-ade_order(h)
    ade_table(h[match(hi,h$feature_id),,drop=FALSE],'Pathway_generation_heterogeneity_discoveries_EN',out)
    exports[[4]]<-ade_heatmaps(r,hi,cfg,out,'Figure42_All_heterogeneous_pathway_heatmap',models='B1')
    exports[[5]]<-ade_forests(r,hi,cfg,out,'Figure42_All_heterogeneous_pathway_forest')
    if(!length(ids))writeLines('No pathway passed primary overall disease BH q<0.05; no nominal substitute set was drawn.',file.path(out,'reports/No_primary_pathway_discoveries_EN.txt'))
  }
  manifest<-ov_bind(exports);ade_table(manifest,'Figure_export_manifest_EN',out)
  stopifnot(all(vapply(seq_len(nrow(manifest)),function(i)all(file.exists(file.path(out,'figures',paste0(manifest$figure[i],c('.pdf','.png','.tiff','.svg'))))),logical(1))))
  invisible(manifest)
}
