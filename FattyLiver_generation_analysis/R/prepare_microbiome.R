# Counts-based preprocessing, under the user-directed input interpretation.
fl_processed_readme_en <- function(d, cfg, data_file) {
  template <- file.path(cfg$project_root, "templates", "processed_README_en.md")
  lines <- readLines(template, encoding = "UTF-8", warn = FALSE)
  values <- list(VERSION = d$version, DATA_FILE = data_file, PROJECT_ROOT = cfg$project_root,
    N_SAMPLES = nrow(d$metadata), N_CASES = sum(d$metadata$fatty_liver),
    N_MINORS = sum(d$metadata$age_group == "Minor"), N_FULL_FEATURES = nrow(d$full$counts),
    N_PRIMARY_FEATURES = nrow(d$matrices$counts), N_FLAGGED = sum(d$metadata$flag_review))
  for (key in names(values)) lines <- gsub(paste0("{{", key, "}}"), as.character(values[[key]]), lines, fixed = TRUE)
  if (any(grepl("{{", lines, fixed = TRUE))) stop("Unresolved English README template field")
  lines
}

fl_processed_set <- function(counts, metadata, annotation, log_pseudocount) {
  fl_qc_validate_matrix(counts)
  relative <- sweep(counts, 2L, metadata$matrix_total, "/")
  transformed <- fl_qc_transforms(counts)
  list(metadata = metadata, features = annotation,
    matrices = list(counts = counts, relative = relative,
      relative_closed = transformed$composition,
      log10 = log10(relative + log_pseudocount), clr = transformed$clr),
    transformation = list(log10_pseudocount = log_pseudocount,
      clr_zero_method = "multiplicative replacement",
      clr_zero_delta = transformed$delta,
      relative_denominator = "all nonzero SGB counts in each sample before prevalence filtering",
      relative_closed_denominator = "retained SGB counts in each sample",
      clr_reference = "geometric mean across this set of retained SGBs"))
}

fl_process_counts <- function(a, x, annotation, cfg) {
  fl_qc_validate_matrix(a)
  # Do not silently round fractional input to satisfy the adopted counts contract.
  if (any(a != floor(a))) stop("Fractional counts found: input counts interpretation must be revisited")
  scope <- fl_three_generation_scope(x)
  x <- scope$cohort
  if (anyDuplicated(x$sample_id) || anyNA(match(x$sample_id, colnames(a)))) stop("Metadata/assay IDs mismatch")
  a <- a[, match(x$sample_id, colnames(a)), drop = FALSE]
  ai <- match(rownames(a), annotation$feature_id)
  if (anyNA(ai)) stop("Feature annotation missing")
  annotation <- annotation[ai, , drop = FALSE]
  sample <- fl_qc_sample_metrics(a, x, cfg$qc_low_mad_multiplier)
  x <- sample$metadata
  filter <- fl_qc_feature_filters(a, x$generation, cfg$qc_sensitivity_prevalence, cfg$qc_prevalence)
  feature <- cbind(annotation, filter$features[, -1L])
  nonzero <- !feature$all_zero
  full_counts <- a[nonzero, , drop = FALSE]
  # Integer storage preserves exact values and reduces the saved package size.
  if (max(full_counts) <= .Machine$integer.max) storage.mode(full_counts) <- "integer"
  full_relative <- fl_qc_composition(full_counts)
  delta <- min(full_relative[full_relative > 0]) / 2
  x$shannon <- -colSums(ifelse(full_relative > 0,
    full_relative * log(pmax(full_relative, .Machine$double.xmin)), 0))
  x$observed_sgbs <- colSums(full_counts > 0)
  retain <- c("sample_id", "participant_id", "family_id", "generation", "age", "age_group", "sex",
    "fat_code", "fatty_liver", "BMI18", "ANTBV18", "dietscore_g18", "smoke18", "alkodpd18",
    "prep_batch", "seq_pool", "PCR_cycles", "DNA_extraction", "stool_ultrasound_gap_days",
    "stool_minus_ultrasound_days", "matrix_total", "detected_features", "flag_low_features",
    "flag_low_total", "flag_review", "automatic_exclusion", "shannon", "observed_sgbs")
  metadata <- x[intersect(retain, names(x))]; rownames(metadata) <- NULL
  full_features <- feature[nonzero, , drop = FALSE]
  main_ids <- feature$feature_id[feature$keep_primary_any_generation]
  selected <- full_counts[match(main_ids, rownames(full_counts)), , drop = FALSE]
  main <- fl_processed_set(selected, metadata, feature[feature$keep_primary_any_generation, ], delta)
  by_generation <- lapply(c("G0", "G1", "G2"), function(g) {
    j <- which(metadata$generation == g)
    ids <- feature$feature_id[feature[[paste0("keep_", g)]]]
    if (!length(j) || !length(ids)) return(NULL)
    counts <- full_counts[match(ids, rownames(full_counts)), j, drop = FALSE]
    fl_processed_set(counts, metadata[j, , drop = FALSE], feature[match(ids, feature$feature_id), ], delta)
  })
  names(by_generation) <- c("G0", "G1", "G2")
  by_generation <- Filter(Negate(is.null), by_generation)
  d <- c(list(schema_version = 1L, analysis_ready = TRUE, readiness_scope = "counts_based_analysis",
    cohort_contract = scope$contract, family_eligibility = scope$summary,
    input_contract = list(type = "nonnegative_integer_counts",
      adopted_from = "user-directed working interpretation, 2026-10-19",
      original_generation_provenance_verified = FALSE,
      metaphlan_native_relative_abundance_verified = FALSE,
      abundance_definition = "counts divided by the full sample SGB count total",
      batch_corrected = FALSE, rarefied = FALSE),
    full = list(counts = full_counts, relative = full_relative, features = full_features),
    by_generation = by_generation, feature_audit = feature,
    filter_summary = filter$summary, sample_thresholds = sample$thresholds,
    analysis_sets = list(all = rep(TRUE, nrow(metadata)), no_review_flags = !metadata$flag_review,
      no_recent_antibiotics = !is.na(metadata$ANTBV18) & metadata$ANTBV18 == 0,
      stool_within_90_days = !is.na(metadata$stool_ultrasound_gap_days) & metadata$stool_ultrasound_gap_days <= 90)), main)
  d
}

fl_processed_publish <- function(d, root, version, activate = TRUE) {
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", version)) stop("Invalid processed version")
  folder <- file.path(root, version)
  if (file.exists(folder)) stop("Processed version exists; never overwritten: ", folder)
  if (!dir.create(folder, recursive = TRUE)) stop("Cannot create processed version directory")
  writeLines("BUILDING", file.path(folder, "STATUS"))
  complete <- FALSE
  on.exit(if (!complete) writeLines("FAILED", file.path(folder, "STATUS")), add = TRUE)
  d$version <- version
  fl_validate_processed_microbiome(d, deep = TRUE)
  path <- file.path(folder, "microbiome_processed.rds")
  saveRDS(d, path, compress = "gzip")
  roundtrip <- readRDS(path)
  fl_validate_processed_microbiome(roundtrip, deep = TRUE)
  if (!identical(roundtrip$metadata$sample_id, d$metadata$sample_id)) stop("Saved package order changed")
  fl_write_csv(data.frame(file = "microbiome_processed.rds", bytes = file.info(path)$size,
    md5 = unname(tools::md5sum(path)), n_samples = nrow(d$metadata),
    n_features_full = nrow(d$full$counts), n_features_primary = nrow(d$matrices$counts)),
    file.path(folder, "dataset_manifest.csv"))
  fl_write_csv(d$feature_audit, file.path(folder, "feature_filter_audit.csv"))
  fl_write_csv(d$filter_summary, file.path(folder, "feature_filter_summary.csv"))
  writeLines("COMPLETED", file.path(folder, "STATUS"))
  # Only this small pointer is updated. Previous version data remain untouched.
  if (activate) writeLines(version, file.path(root, "CURRENT"))
  complete <- TRUE
  list(folder = folder, file = path)
}

fl_build_processed_microbiome <- function(prepared, out_dir, cfg) {
  cat("Reading source counts; processing without depending on QC output files...\n")
  a <- fl_qc_read_matrix(cfg$microbiome_file)$matrix
  x <- prepared$cohort
  md <- fl_read_tsv(cfg$manifest_file, c("sample_id", "prep_batch", "seq_pool", "PCR_cycles", "DNA_extraction"), cfg$manifest_encoding)
  mi <- match(x$sample_id, md$sample_id)
  for (nm in c("prep_batch", "seq_pool", "PCR_cycles", "DNA_extraction")) x[[nm]] <- fl_text(fl_column(md, nm))[mi]
  normalize_family <- function(v) {
    v <- fl_text(v); numeric <- !is.na(v) & grepl("^[0-9]+(\\.0+)?$", v)
    v[numeric] <- sub("\\.0+$", "", v[numeric])
    v[numeric] <- sub("^0+(?=[0-9])", "", v[numeric], perl = TRUE)
    v
  }
  x$family_id <- normalize_family(x$family_id); raw <- normalize_family(x$tutkno80perhe)
  both <- !is.na(x$family_id) & !is.na(raw)
  if (any(x$family_id[both] != raw[both])) stop("Primary family labels conflict")
  x$stool_minus_ultrasound_days <- as.numeric(x$stool_date - x$ultrasound_date)
  annotation <- fl_qc_annotate(rownames(a), cfg$taxonomy_map_file)
  d <- fl_process_counts(a, x, annotation, cfg)
  d$source_paths <- list(counts = cfg$microbiome_file, clinical = cfg$clinical_file,
    manifest = cfg$manifest_file, taxonomy = cfg$taxonomy_map_file)
  d$parameters <- list(prevalence = cfg$qc_prevalence, main_filter = "any generation reaches prevalence threshold",
    single_generation_filter = "that generation reaches prevalence threshold",
    alpha_input = "full nonzero SGB counts before prevalence filtering",
    matrix_orientation = "features x samples", primary_sample_exclusions = "baseline eligibility plus primary families represented in G0/G1/G2",
    low_information_flags = "review only; not automatic exclusions",
    comparison_scale = "use common full-denominator log10 abundances for cross-generation coefficient comparisons",
    clr_note = "generation-specific CLR uses its own feature reference; compare within that generation")
  summary <- do.call(rbind, c(list(data.frame(group = "Overall", n = nrow(d$metadata),
    cases = sum(d$metadata$fatty_liver), features = nrow(d$matrices$counts))),
    lapply(names(d$by_generation), function(g) {
      z <- d$by_generation[[g]]
      data.frame(group = g, n = nrow(z$metadata), cases = sum(z$metadata$fatty_liver), features = nrow(z$matrices$counts))
    })))
  publication <- fl_processed_publish(d, cfg$processed_microbiome_root, cfg$run_id, activate = FALSE)
  fl_write_csv(summary, file.path(out_dir, "processed_data_summary.csv"))
  fl_write_csv(data.frame(item = c("analysis_ready", "readiness_scope", "input_interpretation",
    "original_generation_provenance_verified", "native_MetaPhlAn_relative_abundance_verified", "cohort_rule"),
    value = c("TRUE", "counts_based_analysis", "user-directed counts assumption", "FALSE", "FALSE",
      d$cohort_contract$rule)),
    file.path(out_dir, "processing_contract.csv"))
  inputs <- c(cfg$microbiome_file, cfg$taxonomy_map_file, cfg$clinical_file, cfg$manifest_file, cfg$family_file)
  inputs <- inputs[file.exists(inputs)]
  fl_write_csv(data.frame(file = inputs, md5 = unname(tools::md5sum(inputs))),
    file.path(publication$folder, "source_fingerprints.csv"))
  code <- file.path(cfg$project_root, c("config.R", "R/common.R", "R/entry_helpers.R", "R/microbiome_qc.R",
    "R/prepare_microbiome.R", "R/load_processed_microbiome.R", "scripts/05_prepare_microbiome.R", "tests/processed_microbiome_test.R", "templates/processed_README_en.md"))
  fl_write_csv(data.frame(file = code, md5 = unname(tools::md5sum(code))), file.path(publication$folder, "code_fingerprints.csv"))
  cat("Published processed dataset: ", publication$file, "\n", sep = "")
  lines <- c("# 后续分析统一使用的菌群数据", "",
    "按本次指定的counts输入假设处理；数据已生成并通过验证，可供counts归一化方案的后续分析直接读取。", "",
    "仅纳入G0/G1/G2每代至少一名基线合格成员的主家系。物种检出率过滤、低信息量阈值及转换在此限定队列内计算；未成年人保留。",
    sprintf("- 人数：%d；脂肪肝病例：%d；未成年人：%d。", nrow(d$metadata), sum(d$metadata$fatty_liver), sum(d$metadata$age_group == "Minor")),
    sprintf("- 全部非零特征：%d；三代共用过滤特征：%d。", nrow(d$full$counts), nrow(d$matrices$counts)),
    sprintf("- %d名低信息量复核标记者全部保留。", sum(d$metadata$flag_review)),
    paste0("- 当前版本：", cfg$run_id), paste0("- 数据文件：", publication$file), "",
    "## 数据内容", "",
    "- metadata：临床信息、代际、主家系编号、批次、时间间隔和质控标记；与矩阵列顺序一致。",
    "- matrices$counts：过滤后的原始计数子集；计数值未更改。",
    "- matrices$relative：counts除以过滤前该样本全部SGB计数之和；范围0–1，过滤后列和可能小于1。",
    "- matrices$relative_closed：保留特征内重新归一化，列和为1，供群落距离等使用。",
    "- matrices$log10：对relative加入固定伪计数后取log10，所有代际使用相同伪计数。",
    "- matrices$clr：对保留特征进行乘法零值替换后作自然对数CLR转换。",
    "- full：检出率过滤前的非零counts及relative；metadata中的Shannon和检出特征数由它计算。",
    "- by_generation：分别预先处理好的G0/G1/G2矩阵，采用各代检出率≥10%的特征。",
    "- analysis_sets：全部样本、无复核标记、近期未使用抗生素、采便与超声相隔不超过90天的预先记录索引；主分析仍保留全部样本。", "",
    "## 后续脚本读取方式", "", "```r",
    paste0("source(\"", cfg$project_root, "/R/load_processed_microbiome.R\")"),
    "d <- fl_load_processed_microbiome()",
    "metadata <- d$metadata", "relative <- d$matrices$relative", "clr <- d$matrices$clr",
    "g1 <- fl_load_processed_microbiome(generation = \"G1\")", "```", "",
    "只读取已生成数据，不再重新读SAS或原始菌群RDS。默认读取CURRENT指定版本，并核对文件校验值。各版本保留，CURRENT是唯一会更新的指针。", "",
    "## 口径说明", "",
    "analysis_ready=TRUE仅指counts输入合同下处理完成、可供分析使用；不表示上游生成流程获得独立核验。原始来源未核实的记录仍保留，不能把计算得到的relative称为直接导入的MetaPhlAn原生相对丰度。",
    "不做稀释抽样或自动批次校正。未成年、抗生素使用和复核标记均不导致自动删人。只使用tutkno80perhe。单代CLR参考特征集合不同，不直接比较不同集合的CLR回归系数；代际效应对照使用共同尺度的log10丰度。",
    "数据包含研究编号，仅供本地研究；不要将其作为公开论文Source Data。")
  lines <- append(lines, "[English documentation](README_en.md)", after = 2L)
  for (path in c(file.path(out_dir, "processed_data_readme_zh.md"), file.path(publication$folder, "README_zh.md"))) {
    con <- file(path, "w", encoding = "UTF-8"); writeLines(lines, con); close(con)
  }
  doc <- d; doc$version <- cfg$run_id
  english <- fl_processed_readme_en(doc, cfg, publication$file)
  for (path in c(file.path(out_dir, "README_en.md"), file.path(publication$folder, "README_en.md"))) {
    text <- if (identical(dirname(path), out_dir))
      gsub("(README_zh.md)", "(processed_data_readme_zh.md)", english, fixed = TRUE) else english
    con <- file(path, "w", encoding = "UTF-8"); writeLines(text, con); close(con)
  }
  capture.output(sessionInfo(), file = file.path(publication$folder, "sessionInfo.txt"))
  writeLines("PASSED: in-memory and saved package alignment, counts, normalization, CLR and generation subsets.",
    file.path(publication$folder, "VALIDATION.txt"))
  writeLines(cfg$run_id, file.path(cfg$processed_microbiome_root, "CURRENT"))
  invisible(publication)
}
