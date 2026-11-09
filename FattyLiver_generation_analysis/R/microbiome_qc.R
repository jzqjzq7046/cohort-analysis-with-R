# Base-R QC helpers: no data are read by sourcing this file.
fl_qc_validate_matrix <- function(a) {
  if (!is.matrix(a) || !is.numeric(a) || !nrow(a) || !ncol(a)) stop("Expected numeric feature-by-sample matrix")
  fl_unique_id(rownames(a), "Taxonomy feature ID")
  fl_unique_id(colnames(a), "Taxonomy sample ID")
  if (any(!is.finite(a)) || any(a < 0)) stop("Nonfinite/negative abundance values")
  if (any(colSums(a) <= 0)) stop("Empty microbiome sample; resolve before preprocessing")
  invisible(a)
}

fl_qc_read_matrix <- function(path) {
  # The dense assay can be read without installing the full Bioconductor stack.
  # Explicitly require the named counts slot; never silently pick another assay.
  object <- suppressWarnings(readRDS(path))
  get <- function(z, name) attributes(z)[[name]]
  assays <- get(get(get(object, "assays"), "data"), "listData")
  if (is.null(assays) || !"counts" %in% names(assays)) stop("RDS lacks the expected dense counts assay")
  a <- assays[["counts"]]
  fl_qc_validate_matrix(a)
  list(matrix = a, assay_names = names(assays),
    metadata_keys = names(get(object, "metadata")), object_class = class(object)[1L])
}

fl_qc_annotate <- function(ids, path) {
  con <- bzfile(path, "rt"); on.exit(close(con))
  map <- utils::read.delim(con, header = FALSE, sep = "\t", quote = "", comment.char = "",
    col.names = c("official_id", "taxonomy"), stringsAsFactors = FALSE)
  if (anyDuplicated(map$official_id)) stop("Official map has duplicate keys")
  raw <- sub("^t__", "", ids)
  i <- match(raw, map$official_id); grouped <- is.na(i)
  i[grouped] <- match(paste0(raw[grouped], "_group"), map$official_id)
  data.frame(feature_id = ids, official_id = map$official_id[i], taxonomy = map$taxonomy[i],
    mapping = ifelse(is.na(i), "unmapped", ifelse(grouped, "group_suffix", "exact")),
    stringsAsFactors = FALSE)
}

fl_qc_source_inventory <- function(cfg) {
  files <- list.files(cfg$source_scan_root, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE, include.dirs = FALSE)
  n <- tolower(basename(files))
  raw <- grepl("\\.(fastq|fq)(\\.gz)?$|\\.(bam|sam)$|\\.bowtie2(out|\\.bz2)$", n)
  upstream <- grepl("fastqc|multiqc|kneaddata|knead_data|negative_control|blank_control|sequencing_qc|read_qc", n)
  candidates <- files[grepl("metaphlan.*relative_abundance.*\\.(tsv|txt)(\\.gz)?$", n)]
  header <- lapply(candidates, function(p) {
    con <- if (grepl("\\.gz$", p)) gzfile(p, "rt") else file(p, "rt")
    lines <- tryCatch(readLines(con, n = 4L, warn = FALSE), finally = close(con))
    version <- if (length(lines) && startsWith(lines[1L], "#mpa_")) sub("^#", "", lines[1L]) else "not recorded"
    data.frame(path = p, bytes = file.info(p)$size, database_header = version,
      matches_vJan25 = grepl("vJan25_CHOCOPhlAnSGB_202503", version, fixed = TRUE),
      used_in_analysis = FALSE, stringsAsFactors = FALSE)
  })
  profiles <- if (length(header)) do.call(rbind, header) else data.frame(path = character(),
    bytes = numeric(), database_header = character(), matches_vJan25 = logical(), used_in_analysis = logical())
  hits <- data.frame(path = files[raw | upstream], kind = ifelse(raw[raw | upstream],
    "raw_sequence_or_alignment", "candidate_upstream_report"), stringsAsFactors = FALSE)
  status <- data.frame(item = c("filename_scan", "raw_sequences_or_alignments", "upstream_qc_candidates",
    "matching_relative_profile", "negative_controls", "resequencing_map"),
    n = c(length(files), sum(raw), sum(upstream), sum(profiles$matches_vJan25), NA, NA),
    status = c("completed", if (any(raw)) "located_not_validated" else "not_located",
      if (any(upstream)) "located_not_validated" else "not_located",
      if (any(profiles$matches_vJan25)) "located_requires_reconciliation" else "not_located",
      "not_verified", "not_verified"), stringsAsFactors = FALSE)
  list(profiles = profiles, hits = hits, status = status)
}

fl_qc_lower_cut <- function(z, multiplier = 3) {
  l <- log10(z + 1)
  max(0, 10^(median(l) - multiplier * stats::mad(l)) - 1)
}

fl_qc_sample_metrics <- function(a, x, multiplier) {
  i <- match(x$sample_id, colnames(a)); if (anyNA(i)) stop("Cohort sample absent from assay")
  b <- a[, i, drop = FALSE]
  x$matrix_total <- colSums(b)
  x$detected_features <- colSums(b > 0)
  cut_f <- fl_qc_lower_cut(x$detected_features, multiplier)
  cut_t <- fl_qc_lower_cut(x$matrix_total, multiplier)
  x$flag_low_features <- x$detected_features < cut_f
  x$flag_low_total <- x$matrix_total < cut_t
  x$flag_review <- x$flag_low_features | x$flag_low_total
  x$automatic_exclusion <- FALSE
  list(metadata = x, thresholds = data.frame(metric = c("detected_features", "matrix_total"),
    rule = paste0("log10(x+1) median minus ", multiplier, " MAD; review only"),
    lower_cutoff = c(cut_f, cut_t)))
}

fl_qc_feature_filters <- function(a, generation, thresholds, primary) {
  nonzero <- a > 0
  d <- data.frame(feature_id = rownames(a), prevalence_overall = rowMeans(nonzero))
  for (g in c("G0", "G1", "G2")) {
    keep <- generation == g
    d[[paste0("prevalence_", g)]] <- if (any(keep)) rowMeans(nonzero[, keep, drop = FALSE]) else NA_real_
  }
  d$all_zero <- rowSums(nonzero) == 0
  by_g <- as.matrix(d[paste0("prevalence_", c("G0", "G1", "G2"))])
  d$keep_primary_any_generation <- rowSums(by_g >= primary, na.rm = TRUE) > 0 & !d$all_zero
  for (g in c("G0", "G1", "G2")) d[[paste0("keep_", g)]] <-
    !is.na(d[[paste0("prevalence_", g)]]) & d[[paste0("prevalence_", g)]] >= primary
  summary <- do.call(rbind, lapply(thresholds, function(t) {
    data.frame(threshold = t, overall = sum(d$prevalence_overall >= t),
      G0 = sum(by_g[, 1] >= t, na.rm = TRUE), G1 = sum(by_g[, 2] >= t, na.rm = TRUE),
      G2 = sum(by_g[, 3] >= t, na.rm = TRUE), any_generation = sum(rowSums(by_g >= t, na.rm = TRUE) > 0))
  }))
  list(features = d, summary = summary)
}

fl_qc_composition <- function(a) {
  fl_qc_validate_matrix(a)
  sweep(a, 2L, colSums(a), "/")
}

fl_qc_transforms <- function(composition) {
  # Input is closed on the retained SGB set. No rarefaction is performed.
  p <- fl_qc_composition(composition)
  delta <- min(p[p > 0]) / 2
  replaced <- p
  for (j in seq_len(ncol(p))) {
    zero <- p[, j] == 0; nzero <- sum(zero)
    if (nzero) {
      dj <- min(delta, 0.5 / nzero)
      replaced[zero, j] <- dj
      replaced[!zero, j] <- p[!zero, j] * (1 - nzero * dj)
    }
  }
  logged <- log(replaced)
  clr <- sweep(logged, 2L, colMeans(logged), "-")
  # A single outcome-blind pseudocount, fixed for this exported feature set.
  list(composition = p, log10 = log10(p + delta), clr = clr, delta = delta)
}

fl_qc_summaries <- function(x) {
  groups <- c(fl_groups(x), fl_groups(x, TRUE)[-1L])
  do.call(rbind, lapply(names(groups), function(g) {
    d <- groups[[g]]
    q <- function(v, prob) if (length(v)) unname(stats::quantile(v, prob)) else NA_real_
    data.frame(group = g, n = nrow(d), cases = sum(d$fatty_liver == 1),
      low_features_n = sum(d$flag_low_features), low_total_n = sum(d$flag_low_total),
      review_flag_n = sum(d$flag_review), automatic_exclusions_n = sum(d$automatic_exclusion),
      detected_min = q(d$detected_features, 0), detected_q25 = q(d$detected_features, .25),
      detected_median = q(d$detected_features, .5), detected_q75 = q(d$detected_features, .75),
      detected_max = q(d$detected_features, 1), matrix_total_median = q(d$matrix_total, .5),
      antibiotics_yes = sum(d$ANTBV18 == 1, na.rm = TRUE), antibiotics_missing = sum(is.na(d$ANTBV18)),
      flagged_cases = sum(d$flag_review & d$fatty_liver == 1),
      flagged_antibiotics_available = sum(d$flag_review & !is.na(d$ANTBV18)),
      flagged_antibiotics_yes = sum(d$flag_review & d$ANTBV18 == 1, na.rm = TRUE))
  }))
}

fl_qc_batch_tables <- function(x, field) {
  values <- fl_text(x[[field]]); values[is.na(values)] <- "Missing"
  rows <- lapply(sort(unique(values)), function(v) {
    d <- x[values == v, , drop = FALSE]
    data.frame(field = field, batch = v, n = nrow(d), G0 = sum(d$generation == "G0"),
      G1 = sum(d$generation == "G1"), G2 = sum(d$generation == "G2"),
      cases = sum(d$fatty_liver == 1), controls = sum(d$fatty_liver == 0),
      minors = sum(d$age_group == "Minor"), flagged = sum(d$flag_review),
      detected_median = median(d$detected_features), age_median = median(d$age, na.rm = TRUE))
  })
  overall <- do.call(rbind, rows)
  within <- do.call(rbind, lapply(c("G0", "G1", "G2"), function(g) {
    d <- x[x$generation == g, , drop = FALSE]
    v <- values[x$generation == g]
    do.call(rbind, lapply(sort(unique(v)), function(b) data.frame(field = field, generation = g,
      batch = b, n = sum(v == b), cases = sum(v == b & d$fatty_liver == 1),
      controls = sum(v == b & d$fatty_liver == 0))))
  }))
  list(overall = overall, within = within)
}

fl_qc_cramers <- function(a, b) {
  tab <- table(a, b)
  if (min(dim(tab)) < 2L) return(NA_real_)
  e <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  sqrt(sum((tab - e)^2 / e) / (sum(tab) * min(dim(tab) - 1)))
}

fl_qc_rank <- function(x) {
  groups <- fl_groups(x)
  do.call(rbind, lapply(names(groups), function(g) {
    d <- groups[[g]]
    d <- d[complete.cases(d[c("age", "sex", "prep_batch", "seq_pool")]), , drop = FALSE]
    d$sex <- factor(d$sex); d$prep_batch <- factor(d$prep_batch); d$seq_pool <- factor(d$seq_pool)
    terms <- c(if (length(unique(d$age)) >= 4L) "splines::ns(age, df=3)" else "age",
      if (nlevels(d$sex) > 1L) "sex")
    rows <- lapply(c("none", "prep_batch", "seq_pool", "both"), function(batch) {
      added <- switch(batch, none = character(), prep_batch = "prep_batch", seq_pool = "seq_pool",
        both = c("prep_batch", "seq_pool"))
      added <- added[vapply(added, function(k) nlevels(d[[k]]) > 1L, logical(1))]
      base <- stats::model.matrix(stats::reformulate(c(terms, added)), d)
      with_case <- cbind(base, fatty_liver = d$fatty_liver)
      r <- qr(base)$rank; rr <- qr(with_case)$rank
      data.frame(group = g, batch_terms = batch, n = nrow(d), design_columns = ncol(with_case),
        design_rank = rr, rank_deficient = rr < ncol(with_case),
        fatty_liver_identifiable = rr > r, residual_df = nrow(d) - rr)
    })
    do.call(rbind, rows)
  }))
}

fl_qc_family_batch <- function(x) {
  family <- fl_text(x$family_id); good <- !is.na(family)
  do.call(rbind, lapply(c("prep_batch", "seq_pool"), function(field) {
    f <- family[good]; b <- fl_text(x[[field]])[good]
    keep <- !is.na(b); f <- f[keep]; b <- b[keep]
    by_f <- split(b, f); sizes <- lengths(by_f); multi <- by_f[sizes > 1L]
    pairs <- sum(sizes * (sizes - 1) / 2)
    same <- sum(vapply(by_f, function(v) sum(table(v) * (table(v) - 1) / 2), numeric(1)))
    data.frame(field = field, families_with_two_or_more = length(multi),
      multi_member_families_single_batch = sum(vapply(multi, function(v) length(unique(v)) == 1L, logical(1))),
      within_family_pairs = pairs, same_batch_pairs = same,
      same_batch_pair_pct = if (pairs) 100 * same / pairs else NA_real_)
  }))
}

fl_qc_pca <- function(clr, k = 20L, seed = 20261018L) {
  # Seeded randomized SVD, used only for descriptive technical diagnostics.
  centered <- sweep(clr, 1L, rowMeans(clr), "-")
  k <- min(k, nrow(centered), ncol(centered) - 1L)
  ell <- min(k + 10L, nrow(centered), ncol(centered))
  set.seed(seed)
  q <- qr.Q(qr(centered %*% matrix(rnorm(ncol(centered) * ell), ncol(centered), ell)))
  for (i in 1:2) q <- qr.Q(qr(centered %*% crossprod(centered, q)))
  small <- crossprod(q, centered)
  s <- svd(small, nu = k, nv = 0L)
  scores <- crossprod(small, s$u[, seq_len(k), drop = FALSE])
  colnames(scores) <- paste0("PC", seq_len(k)); rownames(scores) <- colnames(clr)
  list(scores = scores, fraction = sum(s$d[seq_len(k)]^2) / sum(centered^2),
    k = k, seed = seed)
}

fl_qc_batch_pca <- function(scores, x, fraction) {
  do.call(rbind, lapply(c("Overall", "G0", "G1", "G2"), function(g) {
    keep <- if (g == "Overall") rep(TRUE, nrow(x)) else x$generation == g
    d <- x[keep, , drop = FALSE]; y <- scores[keep, , drop = FALSE]
    complete <- complete.cases(d[c("age", "sex", "prep_batch", "seq_pool")])
    d <- d[complete, , drop = FALSE]; y <- y[complete, , drop = FALSE]
    base <- c("splines::ns(age, df=3)", if (length(unique(d$sex)) > 1L) "sex",
      if (g == "Overall") "generation")
    residual_ss <- function(terms) {
      m <- stats::model.matrix(stats::reformulate(terms), d)
      sum(qr.resid(qr(m), y)^2)
    }
    ss <- residual_ss(base)
    out <- lapply(c("prep_batch", "seq_pool", "both"), function(field) {
      add <- if (field == "both") c("prep_batch", "seq_pool") else field
      add <- add[vapply(add, function(k) length(unique(d[[k]])) > 1L, logical(1))]
      fit_ss <- residual_ss(c(base, add))
      data.frame(group = g, batch_terms = field, n = nrow(d),
        partial_R2_top_PCs = max(0, 1 - fit_ss / ss), PC_variance_fraction = fraction,
        inference = "descriptive_only_no_P_values_no_family_independence_assumption")
    })
    do.call(rbind, out)
  }))
}

fl_qc_html <- function(sections, path) {
  esc <- function(z) {
    z <- gsub("&", "&amp;", as.character(z), fixed = TRUE)
    z <- gsub("<", "&lt;", z, fixed = TRUE); gsub(">", "&gt;", z, fixed = TRUE)
  }
  table_html <- function(d) {
    if (!nrow(d)) return("<p>No records.</p>")
    cells <- function(row, tag) paste0("<tr>", paste0("<", tag, ">", esc(row), "</", tag, ">", collapse = ""), "</tr>")
    paste0("<div class='scroll'><table>", cells(names(d), "th"),
      paste(apply(d, 1, cells, tag = "td"), collapse = ""), "</table></div>")
  }
  parts <- c("<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><title>脂肪肝菌群质控</title>",
    "<style>body{font:15px Arial,'Microsoft YaHei',sans-serif;max-width:1400px;margin:30px auto;padding:0 20px;line-height:1.6;color:#263b46}table{border-collapse:collapse;font-size:12px}td,th{padding:7px 10px;border:1px solid #d8e1e6;text-align:left}th{background:#edf3f5}.scroll{overflow:auto}.note{padding:16px;background:#fff2d8;border-left:4px solid #b7832b}h2{margin-top:32px}</style></head><body>",
    "<h1>脂肪肝分析前菌群质控</h1>",
    "<p class='note'><b>矩阵技术检查和预处理已完成；文献同口径的物种相对丰度分析尚未就绪。</b>当前RDS只有counts，来源及与原始相对丰度的关系未核实。导出的读数构成比、对数矩阵、CLR均为临时探索输入，不能声称是已核实的MetaPhlAn物种相对丰度。未因低信息量、抗生素使用或年龄自动删除样本。</p>",
    "<p>仅使用tutkno80perhe，队列仅纳入G0/G1/G2每代至少一名基线合格成员的家系。技术诊断的R²为描述性估计，不是PERMANOVA检验，也不提供家系独立假设下的P值。原始FASTQ、去宿主和阴性对照报告未核验；扫描未找到文件不代表上游没有做质控。</p>")
  for (name in names(sections)) parts <- c(parts, paste0("<h2>", esc(name), "</h2>"), table_html(sections[[name]]))
  con <- file(path, "w", encoding = "UTF-8"); on.exit(close(con))
  writeLines(c(parts, "</body></html>"), con)
}

fl_build_microbiome_qc <- function(prepared, out_dir, cfg) {
  cat("Reading the 2025 taxonomy matrix and checking provenance...\n")
  source <- fl_qc_read_matrix(cfg$microbiome_file); a <- source$matrix
  unit_audit <- data.frame(item = c("object_class", "assay_names", "unit_provenance_metadata_keys",
    "feature_rows", "sample_columns", "missing_values", "nonfinite_values", "negative_values",
    "noninteger_values", "minimum_value", "maximum_value", "sample_total_min",
    "sample_total_median", "sample_total_max", "unit_verified"),
    value = c(source$object_class, paste(source$assay_names, collapse = ";"),
      if (length(source$metadata_keys)) paste(source$metadata_keys, collapse = ";") else "none",
      nrow(a), ncol(a), sum(is.na(a)), sum(!is.finite(a)), sum(a < 0),
      sum(a != floor(a)), min(a), max(a), min(colSums(a)), median(colSums(a)), max(colSums(a)), "NO"))
  if (!setequal(colnames(a), fl_read_tsv(cfg$manifest_file, "sample_id", cfg$manifest_encoding)$sample_id))
    stop("Taxonomy sample set differs from the full manifest")
  inventory <- fl_qc_source_inventory(cfg)
  annotation <- fl_qc_annotate(rownames(a), cfg$taxonomy_map_file)
  x <- prepared$cohort
  normalize_family <- function(v) {
    v <- fl_text(v)
    numeric <- !is.na(v) & grepl("^[0-9]+(\\.0+)?$", v)
    v[numeric] <- sub("\\.0+$", "", v[numeric])
    v[numeric] <- sub("^0+(?=[0-9])", "", v[numeric], perl = TRUE)
    v
  }
  raw_family <- normalize_family(x$tutkno80perhe)
  x$family_id <- normalize_family(x$family_id)
  both <- !is.na(raw_family) & !is.na(x$family_id)
  if (any(raw_family[both] != x$family_id[both])) stop("Primary family map disagrees with clinical labels")
  x$stool_minus_ultrasound_days <- as.numeric(x$stool_date - x$ultrasound_date)
  manifest <- fl_read_tsv(cfg$manifest_file,
    c("sample_id", "prep_batch", "seq_pool", "PCR_cycles", "DNA_extraction"), cfg$manifest_encoding)
  mi <- match(x$sample_id, manifest$sample_id)
  for (nm in c("prep_batch", "seq_pool", "PCR_cycles", "DNA_extraction")) x[[nm]] <- fl_text(fl_column(manifest, nm))[mi]
  s <- fl_qc_sample_metrics(a, x, cfg$qc_low_mad_multiplier); x <- s$metadata
  cohort_matrix <- a[, match(x$sample_id, colnames(a)), drop = FALSE]
  filters <- fl_qc_feature_filters(cohort_matrix, x$generation, cfg$qc_sensitivity_prevalence, cfg$qc_prevalence)
  feature <- cbind(annotation, filters$features[, -1L])
  full_keep <- !feature$all_zero; primary <- feature$keep_primary_any_generation
  composition_full <- fl_qc_composition(cohort_matrix[full_keep, , drop = FALSE])
  transformed <- fl_qc_transforms(cohort_matrix[primary, , drop = FALSE])
  x$shannon_provisional <- -colSums(ifelse(composition_full > 0,
    composition_full * log(pmax(composition_full, .Machine$double.xmin)), 0))
  summaries <- fl_qc_summaries(x)
  batches <- lapply(c("prep_batch", "seq_pool"), function(f) fl_qc_batch_tables(x, f))
  batch_table <- do.call(rbind, lapply(batches, `[[`, "overall"))
  batch_within <- do.call(rbind, lapply(batches, `[[`, "within"))
  design <- fl_qc_rank(x)
  family_batch <- fl_qc_family_batch(x)
  batch_balance <- do.call(rbind, lapply(c("prep_batch", "seq_pool"), function(f) {
    do.call(rbind, lapply(c("Overall", "G0", "G1", "G2"), function(g) {
      d <- if (g == "Overall") x else x[x$generation == g, , drop = FALSE]
      data.frame(field = f, group = g, n = nrow(d),
        cramers_V_generation = if (g == "Overall") fl_qc_cramers(d[[f]], d$generation) else NA_real_,
        cramers_V_case_status = fl_qc_cramers(d[[f]], d$fatty_liver))
    }))
  }))
  batch_variables <- do.call(rbind, lapply(c("prep_batch", "seq_pool", "PCR_cycles", "DNA_extraction"), function(f)
    data.frame(field = f, available = sum(!is.na(x[[f]])), missing = sum(is.na(x[[f]])),
      n_levels = length(unique(x[[f]][!is.na(x[[f]])])))))
  cat("Computing descriptive batch diagnostics (no disease association tests)...\n")
  pca <- fl_qc_pca(transformed$clr, seed = cfg$qc_seed)
  batch_pca <- fl_qc_batch_pca(pca$scores, x, pca$fraction)
  availability <- do.call(rbind, lapply(c("age", "sex", "BMI18", "ANTBV18", "dietscore_g18",
    "smoke18", "alkodpd18", "prep_batch", "seq_pool"), function(f) {
    groups <- c(fl_groups(x), fl_groups(x, TRUE)[-1L])
    do.call(rbind, lapply(names(groups), function(g) {
      d <- groups[[g]]; v <- d[[f]]
      data.frame(field = f, group = g, n = nrow(d), available = sum(!is.na(v)), missing = sum(is.na(v)))
    }))
  }))
  time_groups <- c(fl_groups(x), fl_groups(x, TRUE)[-1L])
  time <- do.call(rbind, lapply(names(time_groups), function(g) {
    d <- time_groups[[g]]; v <- d$stool_ultrasound_gap_days
    data.frame(group = g, n = nrow(d), available = sum(!is.na(v)),
      median_absolute_days = median(v, na.rm = TRUE), within_30_days = sum(v <= 30, na.rm = TRUE),
      within_90_days = sum(v <= 90, na.rm = TRUE), over_90_days = sum(v > 90, na.rm = TRUE),
      over_365_days = sum(v > 365, na.rm = TRUE),
      stool_before_ultrasound = sum(d$stool_minus_ultrasound_days < 0, na.rm = TRUE),
      stool_after_ultrasound = sum(d$stool_minus_ultrasound_days > 0, na.rm = TRUE),
      same_day = sum(d$stool_minus_ultrasound_days == 0, na.rm = TRUE))
  }))
  readiness <- data.frame(item = c("sample_alignment", "primary_family_alignment", "matrix_integrity", "official_taxonomy_map",
    "feature_filtering", "low_information_review", "batch_diagnostics", "count_origin",
    "validated_species_relative_abundance", "upstream_sequence_qc", "literature_matched_association_ready"),
    status = c("passed", "passed", "passed", if (any(feature$mapping == "unmapped")) "partial" else "passed",
      "completed", "flagged_not_excluded", "completed_descriptive", "unverified", "unavailable",
      "not_verified", "NO"))
  tables <- list(matrix_unit_audit = unit_audit, sample_qc_summary = summaries, low_information_thresholds = s$thresholds,
    feature_filter_summary = filters$summary, feature_filter_audit = feature,
    batch_summary = batch_table, batch_by_generation_case = batch_within,
    batch_variable_inventory = batch_variables, batch_balance = batch_balance,
    model_design_rank = design, family_batch_overlap = family_batch,
    batch_pca_descriptive = batch_pca, covariate_availability = availability,
    sample_time_alignment = time, source_relative_profile_inventory = inventory$profiles,
    upstream_file_candidates = inventory$hits, upstream_qc_status = inventory$status,
    analysis_readiness = readiness)
  for (nm in names(tables)) fl_write_csv(tables[[nm]], file.path(out_dir, paste0(nm, ".csv")))
  private <- file.path(dirname(out_dir), "private")
  # Strip birth dates and raw clinical fields from the exported analysis package.
  retain <- c("sample_id", "participant_id", "family_id", "generation", "age", "age_group", "sex",
    "fat_code", "fatty_liver", "BMI18", "ANTBV18", "dietscore_g18", "smoke18", "alkodpd18",
    "prep_batch", "seq_pool", "PCR_cycles", "DNA_extraction", "stool_ultrasound_gap_days",
    "stool_minus_ultrasound_days", "matrix_total",
    "detected_features", "flag_low_features", "flag_low_total", "flag_review", "automatic_exclusion", "shannon_provisional")
  metadata <- x[retain]
  parameters <- list(prevalence = cfg$qc_prevalence, primary_rule = "prevalent in at least one generation",
    per_generation_testing_requires = "use keep_G0 / keep_G1 / keep_G2",
    alpha_feature_set = "all nonzero SGBs before prevalence filtering",
    primary_feature_set = "union of generation-specific prevalence filters",
    unit = "provisional_total_sum_scaled_count_composition_not_verified_MetaPhlAn_relative_abundance",
    multiplicative_zero_delta = transformed$delta, log10_pseudocount = transformed$delta,
    automatic_low_information_exclusion = FALSE, PCA = "randomized SVD, 20 components, two power iterations",
    PC_variance_fraction = pca$fraction, seed = cfg$qc_seed)
  saveRDS(list(analysis_ready = FALSE, readiness = readiness, parameters = parameters,
    metadata = metadata, features = feature, validated_relative_abundance = NULL,
    provisional_count_composition_full = composition_full,
    provisional_count_composition_filtered = transformed$composition,
    provisional_log10_filtered = transformed$log10,
    provisional_clr_filtered = transformed$clr), file.path(private, "microbiome_qc_bundle.rds"))
  saveRDS(list(metadata = metadata, pc_scores = pca$scores), file.path(private, "sample_qc_details.rds"))
  fl_write_csv(metadata[metadata$flag_review, ], file.path(private, "samples_for_review.csv"))
  writeLines(fl_json(parameters), file.path(out_dir, "preprocessing_parameters.json"), useBytes = TRUE)
  fingerprints <- data.frame(source = c(cfg$microbiome_file, cfg$taxonomy_map_file),
    md5 = unname(tools::md5sum(c(cfg$microbiome_file, cfg$taxonomy_map_file))))
  fl_write_csv(fingerprints, file.path(dirname(out_dir), "audit", "microbiome_input_fingerprints.csv"))
  report_sections <- setNames(list(readiness, unit_audit, summaries, s$thresholds, filters$summary,
    batch_variables, batch_table, batch_balance, design, family_batch, batch_pca, time,
    inventory$profiles, inventory$status), c("分析就绪状态", "矩阵数值及单位核查", "样本质控汇总",
    "低信息量复核阈值", "特征过滤数量", "批次变量可用性", "批次与代际及病例分布",
    "批次分布关联强度（描述性）", "协变量设计矩阵秩", "家系与批次重叠",
    "前20主成分的批次解释比例（描述性）", "样本采集与超声时间间隔",
    "原始相对丰度表版本核查", "上游文件核验状态"))
  fl_qc_html(report_sections, file.path(out_dir, "microbiome_qc_report.html"))
  overall <- summaries[summaries$group == "Overall", ]
  lines <- c("# 脂肪肝分析前菌群质控结果", "",
    sprintf("纳入%d人，其中脂肪肝%d例；全部年龄保留。", nrow(x), sum(x$fatty_liver)), "",
    sprintf("- 原始2025版矩阵：%d个特征 × %d份样本；最终队列全零特征%d个。", nrow(a), ncol(a), sum(feature$all_zero)),
    sprintf("- 按每代检出率至少%.0f%%取并集：保留%d个SGB；各代检验使用各自过滤标记。", 100 * cfg$qc_prevalence, sum(primary)),
    sprintf("- 低信息量复核标记%d人；自动删除0人。具体编号只保存在private/samples_for_review.csv。", overall$review_flag_n),
    sprintf("- G2未成年人%d人、脂肪肝%d例，均保留。", sum(x$age_group == "Minor"), sum(x$age_group == "Minor" & x$fatty_liver == 1)),
    sprintf("- 复核样本中脂肪肝%d例；抗生素信息可用%d人，其中%d人报告近期使用。", overall$flagged_cases,
      overall$flagged_antibiotics_available, overall$flagged_antibiotics_yes),
    sprintf("- 批次诊断前%d个主成分覆盖CLR变异的%.1f%%；结果仅描述技术关联，不是PERMANOVA或因果效应。", pca$k, 100 * pca$fraction), "",
    sprintf("- 采便与超声日期均可用%d人，时间差中位数%.0f天；%d人相隔超过90天。", time$available[time$group == "Overall"],
      time$median_absolute_days[time$group == "Overall"], time$over_90_days[time$group == "Overall"]), "",
    "## 必须保留的限制", "",
    "1. 当前RDS仅有counts且没有单位/生成流程元数据。不能凭整数值或列和确认其为原始reads、绝对菌量或标准物种相对丰度。",
    "2. 检出的历史相对丰度表按文件头核对版本，未作为2025版替代。mpa_vJan25...species.txt.bz2是分类注释字典，并非样本丰度表。",
    "3. 已导出总和归一化的读数构成比、log10和CLR临时矩阵，但analysis_ready=FALSE、validated_relative_abundance=NULL；正式文献同口径分析前须获得2025版原始relative_abundance/profile，或数据提供者的转换脚本及说明。",
    "4. 对E:/YFS的文件名扫描未核验上游去接头、低质量过滤、去宿主、阴性/阳性对照及复测流程；未找到报告不能证明上游未做质控。",
    "5. 不将matrix_total当作已核实测序深度；不因低多样性、抗生素使用或未成年身份自动排除。",
    "6. 只使用tutkno80perhe。没有脂肪肝菌群关联检验、预测或中介模型。", "",
    "## 后续处理规则", "",
    "- 保留全部通过矩阵完整性检查的样本，低信息量标记仅用于复核及敏感性分析，不作为当前排除依据。",
    "- 主过滤为任一代检出率至少10%；单代物种检验另用该代的10%标记，未知分类SGB以稳定ID保留。",
    "- 当前prep_batch与seq_pool同时直接作为分类变量时存在秩亏，不能原样照搬两个批次因子的模型。建议主模型使用prep_batch，seq_pool作替代敏感性调整；PCR_cycles非缺失值只有一个水平时不作为协变量。",
    "- 这些批次处理建议以描述性诊断为依据；不对矩阵实施自动批次校正，不把技术关联解释成因果效应。",
    "- 采便与超声时间差超过90天仅作对齐标记，当前不删人；后续关联分析需要预先规定时间窗敏感性分析，不能把不同日期当成完全同期。",
    "- 优先补齐2025版原始物种相对丰度/profile或明确的转换脚本，完成对齐比对并更新单位核验后，再进入文献同口径的菌群关联分析。仅把文件放进目录不会自动解除analysis_ready=FALSE。", "",
    "## 文件", "", "- microbiome_qc_report.html：完整可浏览报告。",
    "- analysis_readiness.csv：各检查及就绪状态。",
    "- feature_filter_audit.csv：所有SGB的官方注释、各代检出率和保留标记。",
    "- batch_by_generation_case.csv / model_design_rank.csv：批次分布及模型可识别性。",
    "- private/microbiome_qc_bundle.rds：含研究编号的本地临时数据包，不作为论文Source Data。", "",
    "方法参考：https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4 ; https://github.com/biobakery/kneaddata",
    "脂肪肝基础分析模板：Kim et al. Nature Metabolism 2025. DOI:10.1038/s42255-025-01318-6")
  con <- file(file.path(out_dir, "microbiome_qc_findings_zh.md"), "w", encoding = "UTF-8")
  writeLines(lines, con); close(con)
  code_files <- file.path(cfg$project_root, c("config.R", "R/common.R", "R/entry_helpers.R",
    "R/microbiome_qc.R", "scripts/04_microbiome_qc.R", "tests/microbiome_qc_test.R"))
  fl_write_csv(data.frame(file = code_files, md5 = unname(tools::md5sum(code_files))),
    file.path(dirname(out_dir), "audit", "code_fingerprints.csv"))
  cat(sprintf("QC complete: %d people, %d primary SGBs, %d review flags. Provisional matrices only.\n",
    nrow(x), sum(primary), overall$review_flag_n))
  invisible(tables)
}
