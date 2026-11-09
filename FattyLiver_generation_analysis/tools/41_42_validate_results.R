# Independent numerical and delivery checks; never fits or changes estimates.
root <- Sys.getenv('FATTYLIVER_ROOT', 'E:/YFS/workspace/FattyLiver_generation_analysis')
source(file.path(root, 'R/generation_ADE_extensions.R'), encoding='UTF-8')
cfg <- ade_setup(root)
x <- tf_load(tf_source(root))
args <- commandArgs(TRUE)
modules <- c(A='41_Generation_heterogeneity', D='42_Functional_pathways')
if (length(args)) modules <- modules[names(modules) %in% args]
checks <- list()
record <- function(kind, test, ok, detail='') {
  checks[[length(checks)+1L]] <<- data.frame(analysis=kind, test=test, pass=isTRUE(ok), detail=detail)
}
for (kind in names(modules)) {
  base <- file.path(root, 'outputs', modules[[kind]], 'current')
  paths <- list.files(base, '^analysis_result.rds$', recursive=TRUE, full.names=TRUE)
  if (!length(paths)) stop('No adopted result for ', kind)
  path <- paths[which.max(file.info(paths)$mtime)]
  out <- dirname(dirname(path)); r <- readRDS(path)
  record(kind, 'original inputs unchanged', all(unname(tools::md5sum(r$protected$path)) == r$protected$md5))
  for (name in names(r$fits)) {
    f <- r$fits[[name]]
    for (what in c('effects', 'heterogeneity')) {
      z <- f[[what]]
      if (!nrow(z)) next
      valid <- tf_bool(z$analysis_valid)
      p <- rep(1, nrow(z)); p[valid] <- z$P_raw[valid]
      q <- p.adjust(p, 'BH', n=nrow(z))
      record(kind, paste(name, what, 'complete-family BH'),
        all(is.na(z$P_BH[!valid])) && isTRUE(all.equal(z$P_BH[valid], q[valid], tolerance=1e-12)))
      record(kind, paste(name, what, 'discovery eligibility'),
        identical(z$FDR_lt_005, valid & !is.na(z$P_BH) & z$P_BH < .05))
      if (what == 'effects') {
        record(kind, paste(name, 'official coefficients and raw P'), all(z$official_result_match[valid]))
        record(kind, paste(name, 'ratio and CI conversion'),
          max(abs(z$abundance_ratio[valid] - 2^z$beta_log2[valid])) < 1e-12 &&
          max(abs(z$ratio_CI_low[valid] - 2^z$CI_low[valid])) < 1e-12 &&
          max(abs(z$ratio_CI_high[valid] - 2^z$CI_high[valid])) < 1e-12)
        ci <- qt(.975, z$df[valid]) * z$SE[valid]
        record(kind, paste(name, 'Satterthwaite t intervals'),
          max(abs((z$CI_high[valid] - z$beta_log2[valid]) - ci)) < 1e-10)
      } else record(kind, paste(name, 'joint two-df test'), all(z$numerator_df[valid] == 2))
    }
    record(kind, paste(name, 'filter coverage'), sum(f$filter$retained) == length(unique(f$effects$feature_id)))
    record(kind, paste(name, 'no repeated participant'), !anyDuplicated(f$sample_ids))
    if (f$mode == 'generation') {
      d <- v3_design(x, x$metadata$sample_id %in% f$sample_ids, if(f$model=='B2')v3_B2 else v3_B1)
      md <- d$metadata; baseX <- cbind(1, as.matrix(d$design[, d$nuisance, drop=FALSE]))
      disease <- md$fatty_liver
      separate <- sapply(c('G0','G1','G2'), function(g)disease*as.integer(md$generation==g))
      traditional <- cbind(disease, separate[,c('G1','G2'),drop=FALSE])
      X1 <- cbind(baseX, separate); X2 <- cbind(baseX, traditional)
      record(kind, paste(name, 'equivalent disease-by-generation mean space'),
        qr(X1)$rank==ncol(X1) && qr(X2)$rank==ncol(X2) &&
        max(abs(qr.resid(qr(X1),X2)))<1e-9 && max(abs(qr.resid(qr(X2),X1)))<1e-9)
    }
  }
  record(kind, 'adult models use identical participants',
    identical(r$fits$Adult_B1_same_sample$sample_ids, r$fits$B2$sample_ids))
  if (kind == 'D') {
    imp <- readRDS(file.path(out, 'private/functional_input.rds'))
    record(kind, 'full community denominator before filtering',
      nrow(imp$relative) == 620L && max(abs(colSums(imp$relative)-1)) < 1e-12)
    record(kind, 'community and generation use same participants',
      identical(r$fits$B1$sample_ids, r$fits$B1_generation$sample_ids))
  }
  manifest <- read.csv(file.path(out, 'tables/Figure_export_manifest_EN.csv'), check.names=FALSE,fileEncoding='UTF-8-BOM')
  record(kind, 'all four formats for every figure', all(vapply(manifest$figure, function(stem)
    all(file.exists(file.path(out, 'figures', paste0(stem, c('.pdf','.png','.tiff','.svg'))))), logical(1))))
  if(kind=='A' && file.exists(file.path(out,'tables/Top10_joint_generation_tests_EN.csv'))){
    top<-read.csv(file.path(out,'tables/Top10_joint_generation_tests_EN.csv'),fileEncoding='UTF-8-BOM')
    allh<-r$heterogeneity[r$heterogeneity$model=='B1'&tf_bool(r$heterogeneity$analysis_valid),,drop=FALSE]
    expected<-head(allh$feature_id[order(allh$P_raw,allh$feature_id)],10)
    record(kind,'top ten selected by full-family raw interaction P',identical(top$feature_id,expected))
    record(kind,'top ten includes q>=0.05 without recoding discovery status',
      nrow(top)==10 && sum(top$FDR_lt_005)==sum(allh$FDR_lt_005) && any(top$P_BH>=.05))
    ov<-read.csv(file.path(out,'tables/All_SGB_heterogeneity_overview_source_EN.csv'),fileEncoding='UTF-8-BOM')
    record(kind,'overview includes all valid primary tests',nrow(ov)==nrow(allh)&&setequal(ov$feature_id,allh$feature_id))
    te<-read.csv(file.path(out,'tables/Top10_generation_effects_EN.csv'),fileEncoding='UTF-8-BOM')
    record(kind,'three effects for every displayed SGB',nrow(te)==30&&all(table(te$feature_id)==3))
    record(kind,'fixed forest domain retains every confidence interval',all(te$ratio_CI_low>=.01)&all(te$ratio_CI_high<=10))
    ax<-read.csv(file.path(out,'tables/Figure41_Fixed_axis_contract_EN.csv'),fileEncoding='UTF-8-BOM')
    record(kind,'literal overview and forest ticks',identical(ax$x_ticks,c('1;2;4;8;16;32','0.01;0.1;1;10')))
    label_file<-file.path(out,'tables/Figure41_Overview_label_manifest_EN.csv')
    if(file.exists(label_file)){
      labs<-read.csv(label_file,fileEncoding='UTF-8-BOM')
      record(kind,'overview labels match all ten forest IDs',
        identical(labs$feature_id,expected)&&identical(labs$label,sub('^t__','',expected)))
      at<-match(labs$feature_id,ov$feature_id)
      record(kind,'ID leader lines anchor to unchanged point coordinates',
        max(abs(labs$point_x-ov$effect_spread_ratio[at]))<1e-10 &&
        max(abs(labs$point_y-ov$minus_log10_BH_q[at]))<1e-10)
      svg<-paste(readLines(file.path(out,'figures',paste0(manifest$figure[1],'.svg')),warn=FALSE),collapse='\n')
      record(kind,'each top-ten ID rendered in both vector panels',all(vapply(labs$label,function(id)
        length(regmatches(svg,gregexpr(id,svg,fixed=TRUE))[[1]])>=2L,logical(1))))
    }
  }
  for (stem in manifest$figure) {
    audit <- function(suffix)jsonlite::read_json(file.path(out,'audit',paste0(stem,suffix)),simplifyVector=TRUE)
    fonts <- audit('.effective-text-audit.json'); collisions <- audit('.collision-audit.json'); align <- audit('.alignment.json')
    # The approved exporter allows 8-pt panel tags; body text retains its
    # 7.01-pt ceiling and is checked independently by text_above_contract.
    record(kind, paste(stem,'physical text size'), fonts$verdict=='PASS' && fonts$minimum_pt>=5 &&
      length(fonts$text_above_contract)==0 && fonts$maximum_pt<=8.01)
    record(kind, paste(stem,'strict collision audit'), collisions$verdict=='PASS' && collisions$summary$fail==0 && collisions$summary$warn==0)
    record(kind, paste(stem,'actual facet alignment'), align$verdict %in% c('PASS','NOT_APPLICABLE_SINGLE_PANEL'))
  }
  pagefiles <- list.files(file.path(out, 'tables'), '_page_manifest_EN.csv$', full.names=TRUE)
  for (p in pagefiles) {
    pg <- read.csv(p,fileEncoding='UTF-8-BOM')
    target <- if (kind=='A') r$heterogeneity else if (grepl('heterogeneous', basename(p))) r$heterogeneity else r$effects
    target <- target[target$model=='B1' & target$FDR_lt_005, , drop=FALSE]
    if (kind=='D' && !grepl('heterogeneous', basename(p))) target <- target[target$generation=='Overall', , drop=FALSE]
    record(kind, paste('complete pagination', basename(p)),
      !anyDuplicated(pg$feature_id) && setequal(pg$feature_id, target$feature_id))
  }
  local <- ov_bind(checks); local <- local[local$analysis==kind, , drop=FALSE]
  fl_write_csv(local, file.path(out, 'audit/independent_validation.csv'))
}
allchecks <- ov_bind(checks)
print(aggregate(pass~analysis,allchecks,function(z)paste(sum(z),'/',length(z),'passed')),row.names=FALSE)
stopifnot(all(allchecks$pass))
cat('PASS: all numerical, source-preservation and delivery checks\n')
