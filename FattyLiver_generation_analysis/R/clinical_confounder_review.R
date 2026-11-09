# Expanded clinical covariate review, not a replacement for the legacy analyses.
# Long et al. 2019: https://pmc.ncbi.nlm.nih.gov/articles/PMC6758911/
# VanderWeele 2019: https://pmc.ncbi.nlm.nih.gov/articles/PMC6447501/
# Covariates are selected by the prespecified causal/measurement role, never P.
# Study adaptations: observed older age/sex composition, contemporary survey
# centre/urbanicity, all-age core plus adult lifestyle and conditional models.

fl_cg_review_older <- function(m, generations) {
  built <- fl_cg_make_comparison(m, generations)
  d <- built$data
  older <- m[m$generation == generations[1L] & !is.na(m$family_id), , drop = FALSE]
  blocks <- split(seq_len(nrow(older)), older$family_id)
  summaries <- if (length(blocks)) do.call(rbind, lapply(names(blocks), function(f) {
    z <- older[blocks[[f]], , drop = FALSE]
    data.frame(family_id = f,
      # All observed older members define exposure; do not silently use a
      # different subset of them to characterize observation opportunity.
      older_mean_age = if (anyNA(z$age)) NA_real_ else mean(z$age),
      older_female_fraction = if (anyNA(z$sex)) NA_real_ else mean(z$sex == "Female"),
      older_age_missing_n = sum(is.na(z$age)),
      older_sex_missing_n = sum(is.na(z$sex)), stringsAsFactors = FALSE)
  })) else data.frame(family_id = character(), older_mean_age = numeric(),
    older_female_fraction = numeric(), older_age_missing_n = integer(),
    older_sex_missing_n = integer())
  k <- match(d$family_id, summaries$family_id)
  for (v in setdiff(names(summaries), "family_id")) d[[v]] <- summaries[[v]][k]
  stopifnot(!anyDuplicated(d$participant_id), !anyNA(k))
  d
}

fl_cg_review_specs <- function(nutrition = TRUE) {
  legacy <- c("age", "sex", "n_older_observed")
  core <- c(legacy, "older_mean_age", "older_female_fraction", "survey_center", "urbanicity")
  life <- c("alcohol", "smoking", "diet_score", "physical_activity", "energy_log")
  drugs <- c("diabetes_med", "lipid_med", "bp_med", "rheum_med")
  adult_complete <- c(core, life, "BMI18", drugs)
  models <- list(
    list("R0_legacy_core_same_sample", "review_core", legacy, FALSE, "core_reassessment", "same_sample_reference", core),
    list("R1_review_primary", "review_core", core, TRUE, "core_reassessment", "all_age_clinical_core", core),
    list("R2_adult_core_same_sample", "adult_review", core, FALSE, "adult_nested", "adult_same_sample_reference", adult_complete),
    list("R3_adult_lifestyle", "adult_review", c(core, life), FALSE, "adult_nested", "adult_lifestyle_extension", adult_complete),
    list("R4_adult_lifestyle_BMI", "adult_review", c(core, life, "BMI18"), FALSE, "adult_nested", "adiposity_conditioned_sensitivity", adult_complete),
    list("R5_adult_lifestyle_BMI_medications", "adult_review", c(core, life, "BMI18", drugs), FALSE, "adult_nested", "adiposity_medication_conditioned_sensitivity", adult_complete))
  if (isTRUE(nutrition)) {
    nutrient <- c("alcohol", "smoking", "physical_activity", "energy_log", "fiber", "sugar")
    models <- c(models, list(
      list("R6_adult_nutrition_core_same_sample", "adult_nutrition", core, FALSE,
        "nutrition_alternative", "nutrition_same_sample_reference", c(core, nutrient)),
      list("R7_adult_nutrition_alternative", "adult_nutrition", c(core, nutrient), FALSE,
        "nutrition_alternative", "nutrient_alternative_sensitivity", c(core, nutrient))))
  }
  do.call(rbind, lapply(models, function(x) data.frame(model = x[[1]], sample_set = x[[2]],
    covariates = paste(x[[3]], collapse = "+"), primary = x[[4]],
    comparison_pair = x[[5]], analysis_role = x[[6]],
    complete_case_fields = paste(x[[7]], collapse = "+"),
    references = "Long2019;VanderWeele2019;LiangZeger1986;Geepack2006;Konis2007",
    study_adaptation = paste("Primary family-label any-case exposure; observed older count/age/sex composition;",
      "survey centre and urbanicity; linear ages to limit complexity; adult complete-case nested comparisons;",
      "BMI/medication conditions may remove shared pathways and are not the total familial association"),
    stringsAsFactors = FALSE)))
}

fl_cg_review_compute <- function(m, cfg) {
  # Input is already scoped and validated by fl_cv_load. Calling the legacy
  # validator would deliberately discard the extended covariates.
  specs <- fl_cg_review_specs(if (is.null(cfg$cg_review_nutrition)) TRUE else cfg$cg_review_nutrition)
  fields <- unique(unlist(strsplit(specs$complete_case_fields, "+", fixed = TRUE)))
  derived <- c("n_older_observed", "older_mean_age", "older_female_fraction")
  required <- unique(c(setdiff(fields, derived), "family_id", "participant_id", "sample_id", "generation",
    "fatty_liver", "age_group"))
  if (!all(required %in% names(m))) stop("Clinical review missing fields: ",
    paste(setdiff(required, names(m)), collapse = ", "))
  if (anyDuplicated(m$participant_id) || anyNA(m$family_id)) stop("Clinical review identity failure")
  presence <- table(m$family_id, factor(m$generation, levels = c("G0", "G1", "G2")))
  if (!nrow(m) || any(presence == 0)) stop("Clinical review requires the fixed three-generation family cohort")
  plans <- results <- cells <- audits <- missingness <- list()
  for (g in cfg$cg_pairs) {
    comparison <- paste(g, collapse = "-")
    d <- fl_cg_review_older(m, g)
    sets <- list()
    for (sample_set in unique(specs$sample_set)) {
      s <- specs[specs$sample_set == sample_set, , drop = FALSE]
      if (length(unique(s$complete_case_fields)) != 1L) stop("Same-sample contract differs")
      need <- strsplit(s$complete_case_fields[1L], "+", fixed = TRUE)[[1L]]
      eligible <- if (grepl("^adult_", sample_set)) !is.na(d$age) & d$age >= 18 else rep(TRUE, nrow(d))
      sets[[sample_set]] <- eligible & complete.cases(d[, need, drop = FALSE])
    }
    for (v in fields) for (stratum in c("All", "Adult", "Minor")) {
      keep <- if (stratum == "All") rep(TRUE, nrow(d)) else !is.na(d$age_group) & d$age_group == stratum
      z <- d[keep, , drop = FALSE]
      missingness[[paste(comparison, v, stratum)]] <- data.frame(comparison = comparison,
        variable = v, age_stratum = stratum, n = nrow(z), missing_n = sum(is.na(z[[v]])),
        observed_n = sum(!is.na(z[[v]])), observed_cases = sum(z$fatty_liver[!is.na(z[[v]])]),
        stringsAsFactors = FALSE)
    }
    for (i in seq_len(nrow(specs))) {
      s <- specs[i, , drop = FALSE]
      keep <- sets[[s$sample_set]]
      z <- droplevels(d[keep, , drop = FALSE])
      fit <- fl_cg_fit(z, s, comparison, cfg)
      r <- fit$result
      r$review_primary <- s$primary
      r$analysis_role <- s$analysis_role
      r$comparison_pair <- s$comparison_pair
      r$complete_case_fields <- s$complete_case_fields
      r$conditioned_fields <- s$covariates
      r$cohort_n <- nrow(m)
      r$eligible_comparison_n <- nrow(d)
      r$excluded_n <- nrow(d) - nrow(z)
      r$excluded_cases <- sum(d$fatty_liver[!keep])
      r$excluded_minors_n <- sum(d$age_group[!keep] == "Minor", na.rm = TRUE)
      r$family_scope <- "Fixed families observed in all three generations before covariate exclusion"
      key <- paste(comparison, s$model)
      results[[key]] <- r
      cells[[key]] <- fit$cells
      plans[[key]] <- cbind(data.frame(comparison = comparison, stringsAsFactors = FALSE), s)
    }
    for (sample_set in unique(specs$sample_set)) {
      names_in_set <- specs$model[specs$sample_set == sample_set]
      reference_ids <- sort(as.character(d$participant_id[sets[[sample_set]]]))
      # Verify from each model's actual complete-case requirements, not only a
      # repeated comparison of the same stored logical vector.
      sample_matches <- vapply(names_in_set, function(model) {
        s <- specs[specs$model == model, , drop = FALSE]
        cov <- strsplit(s$covariates, "+", fixed = TRUE)[[1L]]
        ids <- sort(as.character(d$participant_id[sets[[sample_set]] &
          (if (length(cov)) complete.cases(d[, cov, drop = FALSE]) else rep(TRUE, nrow(d)))]))
        identical(ids, reference_ids)
      }, logical(1))
      audits[[paste(comparison, sample_set)]] <- data.frame(comparison = comparison,
        sample_set = sample_set, models = paste(names_in_set, collapse = ";"),
        same_participants = all(sample_matches), n = length(reference_ids),
        cases = sum(d$fatty_liver[sets[[sample_set]]]),
        families = length(unique(d$family_id[sets[[sample_set]]])),
        complete_case_fields = specs$complete_case_fields[match(sample_set, specs$sample_set)],
        stringsAsFactors = FALSE)
      if (!all(sample_matches)) stop("Clinical review same-sample contract failed")
    }
  }
  out <- list(review_models = do.call(rbind, results), review_model_cells = do.call(rbind, cells),
    review_model_plan = do.call(rbind, plans), review_same_sample_audit = do.call(rbind, audits),
    review_missingness = do.call(rbind, missingness))
  out <- lapply(out, function(x) { rownames(x) <- NULL; x })
  if (sum(out$review_models$review_primary) != length(cfg$cg_pairs)) stop("Review primary contract failed")
  out
}

fl_cg_review_note <- function(results, path) {
  r <- results$review_models
  if (is.null(r)) stop("Missing clinical review results")
  fmt <- function(x, digits = 3L) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
  rows <- function(d) vapply(seq_len(nrow(d)), function(i) paste0("| ", d$comparison[i], " | ",
    d$model[i], " | ", d$n[i], "/", d$cases[i], " | ", fmt(d$OR[i]), " (",
    fmt(d$CI_low[i]), ", ", fmt(d$CI_high[i]), ") | ",
    if (is.na(d$p_value[i])) "NA" else format(d$p_value[i], digits = 4L), " | ", d$status[i], " |"), character(1))
  lines <- c("# 脂肪肝跨代关联：混杂因素审阅", "",
    paste0("研究范围仍为基线已观察到 G0/G1/G2 的", results$cohort_summary$families,
      "个主家系；后续协变量缺失或成人限制不重新定义家庭资格。只使用主家系标签，不推定已核验的直接亲子或祖孙关系。"), "",
    "每项比较中，年轻成员每人一行；暴露是已观察到的对应年长一代成员中是否至少一人有超声脂肪肝。Logistic GEE、独立工作相关、按家系聚类的稳健标准误；所有P为原始双侧P，不作多重比较校正。", "",
    "## 预设的全年龄主要模型", "",
    "调整年轻成员年龄、性别、年长成员观察人数、年长成员平均年龄和女性比例、调查地区和城乡。人数等变量不变时仅按预设规则移除不变项，并在结果中记录。年龄采用线性项，以限制稀少事件下的复杂度；这并不保证完全控制年龄非线性。", "",
    "年长成员年龄和性别构成用于补充 any-case 暴露的观察机会差异，属于本研究适配；并非原始文献逐字复现的模型。", "",
    "| 比较 | 模型 | 人数/病例 | OR (95%CI) | 原始P | 状态 |",
    "|---|---|---:|---|---:|---|", rows(r[r$review_primary, ]), "",
    "## 成人扩展与同样本比较", "",
    "成人 core、生活方式、生活方式+BMI、生活方式+BMI+药物四层使用同一批资料完整者。生活方式项为饮酒（份/天）、吸烟、饮食评分、身体活动、对数能量摄入。进一步条件于BMI及糖尿病/降脂/降压/风湿病药物时，可能控制共同代谢途径或疾病后的治疗，因此这些结果是条件敏感性分析，不作为家庭总体关联的替代。药物仅按现有大类解释，不能改称具体药物如二甲双胍或PPI。", "",
    "| 比较 | 模型 | 人数/病例 | OR (95%CI) | 原始P | 状态 |",
    "|---|---|---:|---|---:|---|", rows(r[r$sample_set == "adult_review", ]), "",
    "营养替代分析在其自身同样本中，以膳食纤维、糖及能量替代饮食评分，同时调整饮酒、吸烟和身体活动；不能把它与不同样本模型的变化直接归因于某个营养素。未成年人问卷缺失不填为零或从全年龄主要模型中剔除。", "",
    "近期抗生素主要影响同期菌群测量，未将其默认作为脂肪肝家庭聚集的核心共同原因。菌群模块另作抗生素调整/排除分析。肝酶、血脂、血糖、胰岛素等并非因存在就全部加入；其中一些可能是结局相关代谢表型、共同途径或后果，全部条件化会改变所估计的关系。", "",
    "低事件/参数、分离、秩亏及不收敛均保留并报告，不按显著性删变量或自动换用忽略家系结构的Firth模型。G2原始仅17例，成人完整资料事件更少；即使得出数值，置信区间和小样本限制仍须优先解释。", "",
    "## 文献与适配边界", "",
    "- Long et al., Liver International 2019：跨代脂肪肝关联及年龄、性别、饮酒、BMI分层调整框架。https://pmc.ncbi.nlm.nih.gov/articles/PMC6758911/",
    "- VanderWeele 2019：依据因果/共同原因角色选择混杂因素，不按单变量P选择；区别混杂因素、途径变量与碰撞点。https://pmc.ncbi.nlm.nih.gov/articles/PMC6447501/",
    "- Liang & Zeger 1986：GEE框架；本研究使用家系聚类稳健方差。上述完整协变量清单、老一代观察构成和同样本比较是预先记录的研究适配，不能称为上述论文完全相同的分析。")
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
  invisible(path)
}
