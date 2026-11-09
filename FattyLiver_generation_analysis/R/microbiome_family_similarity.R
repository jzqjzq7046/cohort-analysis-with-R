# Beghini et al., Nature 2025, Methods: Microbiome null permutations.
# See the module README for every departure from the source studies.
fl_ms_require_packages <- function() {
  if (!requireNamespace("vegan", quietly = TRUE)) stop("vegan is required; no automatic installation")
  invisible(TRUE)
}

fl_ms_blocks <- function(m, scheme, width) {
  required <- c("generation", "sex", "age", "family_id")
  if (grepl("prep_batch$", scheme)) required <- c(required, "prep_batch")
  if (!all(required %in% names(m))) stop("Required permutation metadata missing")
  if (length(width) != 1L || !is.finite(width) || width <= 0) stop("Invalid age-band width")
  text <- function(x) { x <- as.character(x); x[is.na(x) | !nzchar(trimws(x))] <- "<missing>"; x }
  age_band <- ifelse(is.na(m$age), "<missing>", as.character(floor(m$age / width)))
  columns <- list(text(m$generation), text(m$sex), age_band)
  if (grepl("prep_batch$", scheme)) columns <- c(columns, list(text(m$prep_batch)))
  key <- do.call(interaction, c(columns, list(drop = TRUE, lex.order = TRUE)))
  groups <- unname(split(seq_len(nrow(m)), key))
  movable <- cross_family_movable <- rep(FALSE, nrow(m))
  for (z in groups) {
    movable[z] <- length(z) > 1L
    cross_family_movable[z] <- length(unique(m$family_id[z])) > 1L
  }
  list(groups = groups, movable = movable, cross_family_movable = cross_family_movable,
    summary = data.frame(scheme = scheme, n = nrow(m), blocks = length(groups),
      singleton_blocks = sum(lengths(groups) == 1L), movable_people = sum(movable),
      cross_family_movable_people = sum(cross_family_movable),
      age_missing_n = sum(is.na(m$age)), sex_missing_n = sum(is.na(m$sex)),
      prep_batch_missing_n = if ("prep_batch" %in% required) sum(is.na(m$prep_batch)) else NA_integer_))
}

fl_ms_permute_nodes <- function(groups, n) {
  p <- seq_len(n)
  for (z in groups) if (length(z) > 1L) p[z] <- z[sample.int(length(z))]
  p
}

fl_ms_plan_pairs <- function(m, D, generations) {
  a <- which(m$generation == generations[1L]); b <- which(m$generation == generations[2L])
  if (!length(a) || !length(b)) stop("A generation is absent")
  same <- outer(m$family_id[a], m$family_id[b], "==")
  ij <- which(same, arr.ind = TRUE)
  e1 <- a[ij[, 1L]]; e2 <- b[ij[, 2L]]
  if (!length(e1) || all(same)) stop("Need both within- and between-family pairs")
  values <- D[a, b, drop = FALSE]
  summarize <- function(v, label) data.frame(comparison = paste(generations, collapse = "-"),
    pair_type = label, pairs = length(v), mean_distance = mean(v),
    median_distance = median(v), Q1 = unname(quantile(v, .25)), Q3 = unname(quantile(v, .75)))
  within <- values[same]; between <- values[!same]
  list(comparison = paste(generations, collapse = "-"), a = a, b = b, e1 = e1, e2 = e2,
    n_within = length(within), n_between = length(between), total_distance = sum(values),
    within_mean = mean(within), between_mean = mean(between),
    summary = rbind(summarize(within, "within_primary_family"), summarize(between, "between_primary_families")))
}

fl_ms_compute <- function(m, relative, cfg) {
  if (!is.numeric(cfg$ms_permutations) || cfg$ms_permutations < 99L ||
      cfg$ms_permutations != floor(cfg$ms_permutations)) stop("Use at least 99 integer permutations")
  scope <- fl_three_generation_scope(m)
  if (!all(scope$keep)) stop("Input includes incomplete families; regenerate restricted data with 05")
  m <- scope$cohort
  if (anyDuplicated(m$sample_id) || !identical(colnames(relative), m$sample_id) ||
      !is.matrix(relative) || any(!is.finite(relative)) || any(relative < 0) ||
      max(abs(colSums(relative) - 1)) > 1e-10) stop("Invalid full-denominator relative matrix or sample alignment")
  D <- as.matrix(vegan::vegdist(t(relative), method = "bray"))
  if (any(!is.finite(D)) || max(abs(D - t(D))) > 1e-12 || max(abs(diag(D))) > 1e-12 ||
      any(D < -1e-12 | D > 1 + 1e-12)) stop("Invalid Bray-Curtis matrix")
  plans <- lapply(cfg$ms_pairs, function(g) fl_ms_plan_pairs(m, D, g))
  names(plans) <- vapply(plans, `[[`, character(1), "comparison")
  results <- nulls <- block_audits <- list()
  B <- as.integer(cfg$ms_permutations)
  for (k in seq_along(cfg$ms_schemes)) {
    scheme <- cfg$ms_schemes[k]
    blocks <- fl_ms_blocks(m, scheme, cfg$ms_age_band_years)
    block_audits[[k]] <- blocks$summary
    set.seed(cfg$ms_seed + k - 1L)
    null <- matrix(NA_real_, B, length(plans), dimnames = list(NULL, names(plans)))
    cat("Node permutation scheme: ", scheme, "; B = ", B,
      "; cross-family movable people = ", sum(blocks$cross_family_movable), "\n", sep = "")
    for (r in seq_len(B)) {
      p <- fl_ms_permute_nodes(blocks$groups, nrow(m))
      for (j in seq_along(plans)) {
        z <- plans[[j]]
        null[r, j] <- mean(D[cbind(p[z$e1], p[z$e2])])
      }
      if (r %% 2000L == 0L) cat("  ", r, "/", B, " permutations completed\n", sep = "")
    }
    for (j in seq_along(plans)) {
      z <- plans[[j]]; v <- null[, j]; observed <- z$within_mean
      valid <- sd(v) > 1e-12
      # Inclusive tails and +1 avoid zero Monte Carlo P values. No P adjustment.
      lower <- (1 + sum(v <= observed + 1e-12)) / (B + 1)
      upper <- (1 + sum(v >= observed - 1e-12)) / (B + 1)
      null_id <- paste(scheme, z$comparison, sep = ":")
      results[[null_id]] <- data.frame(comparison = z$comparison, scheme = scheme,
        n_people = length(unique(c(z$a, z$b))), families = length(unique(m$family_id)),
        within_pairs = z$n_within, between_pairs = z$n_between,
        observed_within_mean = observed, observed_between_mean = z$between_mean,
        observed_within_minus_between = observed - z$between_mean,
        null_mean = mean(v), observed_minus_null_mean = observed - mean(v),
        null_reference_low = unname(quantile(v, .025)), null_reference_high = unname(quantile(v, .975)),
        null_SD = sd(v), p_raw_two_sided = if (valid) min(1, 2 * min(lower, upper)) else NA_real_,
        permutations = B, distinct_null_statistics = length(unique(round(v, 12))),
        within_edges_both_endpoints_not_cross_family_movable = sum(!blocks$cross_family_movable[z$e1] &
          !blocks$cross_family_movable[z$e2]),
        status = if (valid) "ESTIMATED" else "NO_EFFECTIVE_PERMUTATION_VARIATION",
        p_adjustment = "none; raw two-sided", interval_type = "95% null reference, NOT effect confidence interval")
      nulls[[null_id]] <- data.frame(comparison = z$comparison, scheme = scheme,
        permutation = seq_len(B), within_mean_distance = v)
    }
  }
  out <- list(results = do.call(rbind, results), pair_summary = do.call(rbind, lapply(plans, `[[`, "summary")),
    exchangeability = do.call(rbind, block_audits), null_distribution = do.call(rbind, nulls),
    cohort = data.frame(n = nrow(m), families = scope$contract$families,
      cases = sum(m$fatty_liver), minors = sum(m$age < 18, na.rm = TRUE),
      features = nrow(relative), age_band_years = cfg$ms_age_band_years))
  for (name in names(out)) rownames(out[[name]]) <- NULL
  out
}

fl_ms_report <- function(r, out_dir, version) {
  main <- r$results[r$results$scheme == "generation_sex_age", , drop = FALSE]
  batch <- r$results[r$results$scheme == "generation_sex_age_prep_batch", , drop = FALSE]
  lines <- c("# 07：三代家庭菌群相似性结果", "", "## Material Passport", "",
    "- Origin Skill: academic-research-suite / experiment-agent; nature-academic-search; statistical-analysis.",
    "- Verification Status: ANALYZED — completed computation; see audit for independent checks.",
    paste0("- Input version: ", version), "- Date: 2026-10-19.", "",
    sprintf("纳入%d个家庭、%d人（未成年人%d人）；使用全部%d个非零SGB的完整分母相对丰度计算Bray–Curtis。",
      r$cohort$families, r$cohort$n, r$cohort$minors, r$cohort$features), "",
    "距离越小表示物种组成越相似。以下配对数不是独立样本量；P值来自保留节点重复关系的个体置换，非将配对行视为独立样本的检验。", "",
    "## 原始总体相似性参照（不作为充分调整结论）", "",
    "| 比较 | 家系内配对 | 家系外配对 | 家系内平均距离 | 家系外平均距离 | 主置换零分布均值 | 原始双侧P | 加建库批次限制P |",
    "|---|---:|---:|---:|---:|---:|---:|---:|")
  for (i in seq_len(nrow(main))) {
    z <- main[i, ]; b <- batch[match(z$comparison, batch$comparison), ]
    lines <- c(lines, sprintf("| %s | %d | %d | %.4f | %.4f | %.4f | %.4g | %.4g |",
      z$comparison, z$within_pairs, z$between_pairs, z$observed_within_mean,
      z$observed_between_mean, z$null_mean, z$p_raw_two_sided, b$p_raw_two_sided))
  }
  lines <- c(lines, "", "主交换范围为代际×性别×5岁年龄层；技术敏感性再限制建库批次。均为9999次（自检时使用较少次数）节点置换，报告未经多重校正的原始双侧P。", "",
    "主结果及敏感性结果应同时阅读，不能按较小P值挑选方案。限制置换依赖零假设下的层内可交换性；5岁分层并非连续年龄完全调整，未测共居、饮食及技术因素仍可能影响观察关联。", "",
    "`null_reference_low/high`是零假设置换参考范围，不是效应置信区间。`observed_minus_null_mean`是相对该条件零假设的偏离，不是回归系数。批次限制后的可移动人数和冻结边数见exchangeability.csv及results.csv。", "",
    "家系标签不证明直接亲子/祖孙关系；相似性不能解释为菌株共享、垂直传播、遗传或脂肪肝的菌群中介。本模块没有按脂肪肝或物种显著性筛选，也不代表已经进行了脂肪肝菌群关联分析。", "",
    "方法依据及全部本数据适配见scripts/07_cross_generation_microbiome_similarity_README.md。",
    "核心来源：[Beghini et al., Nature 2025](https://www.nature.com/articles/s41586-024-08222-1)、[Valles-Colomer et al., Nature Microbiology 2022](https://www.nature.com/articles/s41564-021-01021-8)。")
  if (!is.null(r$adjusted_results)) {
    ar <- r$adjusted_results
    review <- c("## 全面协变量审阅：优先阅读", "",
      "新增结果是双端协变量回归调整的家系内外Bray–Curtis差异，使用Freedman–Lane残差MRQAP和学生化t。负系数表示控制所列协变量后家系内距离较小。",
      "全龄core_base/core_adjusted严格同样本；成人四模型严格同样本；营养两模型严格同样本。原319三代家庭先确定，完整资料和成人子集不重新定义家系资格。",
      "本轮9种方案在结果产生前固定。全部报告，不按P选择；BMI、糖尿病与药物模型是条件敏感性，不宣称因果直接效应。",
      "P=(1+sum(abs(t置换)>=abs(t观察)))/(9999+1)，9999次时最小0.0001；与原参照双尾加倍算法不同。无独立配对OLS置信区间。", "",
      "| 方案 | 代际 | 人数 | IF配对 | 调整距离差 | 原始双侧P | 状态 |",
      "|---|---|---:|---:|---:|---:|---|")
    for (i in seq_len(nrow(ar))) {
      z <- ar[i, ]
      review <- c(review, sprintf("| %s | %s | %d | %s | %s | %s | %s |", z$scenario,
        z$comparison, z$n_people, format(z$within_pairs, trim = TRUE),
        format(signif(z$adjusted_distance_difference, 5), trim = TRUE),
        format(signif(z$p_raw_two_sided, 5), trim = TRUE), z$status))
    }
    review <- c(review, "", "双端必调：年龄自然样条3df、性别、建库批次、log10输入矩阵总量及同批次指标。core新增调查中心、城乡、近期抗生素、便秘、稀便、乳酸菌制剂、季节周期、采集至记录冷冻间隔及粪便与超声时间差；控制调查中心时另加同中心指标。",
      "成人生活方式进一步含吸烟、饮酒、饮食评分、体力活动及能量；BMI和糖尿病/药物分层加入；营养替代用能量/纤维/糖代替综合饮食分数。无近期扰动方案排除近期抗生素或胃肠炎。",
      "仅在角色内性别×5岁年龄层进行节点残差置换，仍依赖条件可交换性；不能保证消除未测共居、精细饮食及地理混杂。矩形二模网络是本研究适配。",
      "方法：[Dekker et al. 2007 §3.2.1](https://pmc.ncbi.nlm.nih.gov/articles/PMC2798974/)、[Beghini et al. Nature 2025](https://www.nature.com/articles/s41586-024-08222-1)。", "")
    lines <- append(lines, review, after = 2L)
  }
  con <- file(file.path(out_dir, "results_zh.md"), "w", encoding = "UTF-8"); writeLines(lines, con); close(con)
  escape <- function(x) { x <- gsub("&", "&amp;", as.character(x), fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }
  html <- c("<!doctype html><html><head><meta charset='utf-8'><title>Familial microbiome similarity</title>",
    "<style>body{font:14px Arial,sans-serif;margin:25px;line-height:1.6}table{border-collapse:collapse}th,td{padding:6px;border-bottom:1px solid #ddd;white-space:nowrap}.scroll{overflow:auto}</style></head><body>",
    "<h1>三代家庭菌群相似性 / Familial microbiome similarity</h1>",
    "<p>Bray–Curtis + constrained node permutations; raw two-sided P. Pair counts are not independent sample sizes. Null reference intervals are not effect confidence intervals.</p>",
    "<p>物种相似性不等于菌株传播。严格文献框架与本研究适配详见方法README；结果解释见results_zh.md。</p>")
  if (!is.null(r$adjusted_results)) html <- c(html,
    "<h2>全面协变量审阅 / Covariate review: primary interpretation</h2>",
    "<p>双端协变量调整的距离差；Freedman–Lane残差MRQAP，学生化t，原始双侧P。9999次时最小P=0.0001。原始总体置换仅作参照；完整方法与中文解释见results_zh.md。</p>")
  display <- c(if (!is.null(r$adjusted_results)) c("adjusted_results", "adjusted_pair_summary", "adjusted_exchangeability"),
    "cohort", "results", "pair_summary", "exchangeability")
  for (name in display) {
    d <- r[[name]]
    html <- c(html, paste0("<h2>", name, "</h2><div class='scroll'><table><tr>",
      paste0("<th>", escape(names(d)), "</th>", collapse = ""), "</tr>"))
    for (i in seq_len(nrow(d))) html <- c(html, paste0("<tr>", paste0("<td>",
      escape(vapply(d, function(x) if (is.numeric(x)) format(signif(x[i], 5), trim = TRUE) else as.character(x[i]), character(1))), "</td>", collapse = ""), "</tr>"))
    html <- c(html, "</table></div>")
  }
  con <- file(file.path(out_dir, "microbiome_family_similarity.html"), "w", encoding = "UTF-8")
  writeLines(c(html, "</body></html>"), con); close(con)
}

fl_ms_code_paths <- function(cfg) file.path(cfg$project_root, c("R/approved_design.R", "config/approved_design_20261020.md", "config.R", "config/microbiome_family_similarity.R",
  "R/common.R", "R/entry_helpers.R", "R/load_processed_microbiome.R", "R/microbiome_family_similarity.R",
  "R/confounders.R", "R/microbiome_similarity_adjusted.R", "scripts/mrqap_permutation_kernel.py",
  "config/confounder_review_plan_20261019.md", "tests/microbiome_similarity_adjusted_test.R",
  "scripts/07_cross_generation_microbiome_similarity.R", "scripts/07_cross_generation_microbiome_similarity_README.md",
  "tests/microbiome_family_similarity_test.R"))

fl_ms_run <- function(cfg) {
  context <- fl_output_start_run(cfg); cfg <- context$cfg
  audit <- file.path(context$run_dir, "audit"); out <- file.path(context$run_dir, "similarity")
  private <- file.path(context$run_dir, "private")
  dir.create(audit); dir.create(out); dir.create(private); writeLines("RUNNING", file.path(context$run_dir, "run_status.txt"))
  tryCatch({
    d <- fl_load_processed_microbiome(cfg$project_root, version = cfg$ms_data_version)
    version <- d$version; m <- d$metadata; relative <- d$full$relative
    rm(d); gc(FALSE)
    cat("07 input: ", nrow(m), " people; ", nrow(relative), " SGBs; version ", version, "\n", sep = "")
    results <- fl_ms_compute(m, relative, cfg)
    results$results$result_role <- "legacy_overall_similarity_reference; not fully confounder adjusted"
    reviewed <- fl_cv_load(cfg, m)
    adjusted <- fl_ms_adjusted_compute(reviewed$metadata, relative, cfg, private)
    results <- c(results, adjusted)
    fl_write_csv(data.frame(covariate_version = reviewed$version), file.path(audit, "covariate_version.csv"))
    for (name in names(results)) {
      destination <- if (name %in% c("null_distribution", "adjusted_null_distribution")) audit else out
      fl_write_csv(results[[name]], file.path(destination, paste0(name, ".csv")))
    }
    fl_ms_report(results, out, version)
    input <- file.path(cfg$processed_microbiome_root, version, "microbiome_processed.rds")
    paths <- unique(c(input, reviewed$source_files, fl_ms_code_paths(cfg)))
    fl_write_csv(data.frame(path = paths, md5 = unname(tools::md5sum(paths))), file.path(audit, "input_code_fingerprints.csv"))
    fl_write_csv(data.frame(data_version = version, family_rule = "G0_G1_G2_each_at_least_one",
      distance = "Bray-Curtis; full-denominator relative; all nonzero SGBs",
      permutations = cfg$ms_permutations, seed = cfg$ms_seed, age_band_years = cfg$ms_age_band_years,
      p_value_policy = "raw two-sided; no multiplicity adjustment", privacy = "aggregate tables only"),
      file.path(audit, "analysis_contract.csv"))
    fl_write_csv(data.frame(scenarios = paste(cfg$ms_review_scenarios, collapse = ";"),
      estimator = "OLS distance difference; dual endpoint covariates; no independent-dyad P/CI",
      inference = "Freedman-Lane reduced residual MRQAP; absolute studentized t; add-one Monte Carlo P",
      exchangeability = "separate generation roles; sex x 5-year age band",
      covariate_version = reviewed$version, permutations = cfg$ms_permutations,
      first5_lm_max_error = max(results$adjusted_results$first5_lm_statistic_max_error, na.rm = TRUE)),
      file.path(audit, "adjusted_analysis_contract.csv"))
    capture.output(sessionInfo(), file = file.path(audit, "sessionInfo.txt"))
    writeLines(c("COMPLETED", "Read exchangeability diagnostics and methodological limitations before interpretation."),
      file.path(context$run_dir, "run_status.txt"))
    final <- fl_output_complete_run(context)
    cat("Completed 07: ", file.path(final, "similarity", "microbiome_family_similarity.html"), "\n", sep = "")
    invisible(results)
  }, error = function(e) { fl_output_fail_run(context, conditionMessage(e)); stop(e) })
}
