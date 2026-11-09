# Scientific covariate scenarios fixed in config/confounder_review_plan_20261019.md.
# Family eligibility is fixed before any scenario-specific complete-case selection.
fl_dv_review_compute <- function(m, D, cfg, private, out, audit) {
  scenarios <- fl_cv_scenarios(m)
  stopifnot(identical(scenarios$core_base$mask, scenarios$core_adjusted$mask),
    identical(scenarios$adult_core$mask, scenarios$adult_BMI$mask),
    identical(scenarios$bowel_base$mask, scenarios$bowel_adjusted$mask),
    identical(scenarios$DNA_base$mask, scenarios$core_DNA$mask))
  alpha <- hetero <- diagnostics <- beta <- exchange <- flow <- list()
  main <- NULL; fullD <- as.matrix(D)
  for (s in names(scenarios)) {
    z <- scenarios[[s]]; j <- which(z$mask); ms <- m[j, , drop = FALSE]
    cat("08 confounder review: ", s, "; n=", length(j), "; cases=", sum(ms$fatty_liver), "\n", sep = ""); flush.console()
    cf <- cfg; cf$dv_extra_covariates <- z$covariates; cf$dv_review_mode <- TRUE
    folder <- file.path(private, paste0("review_", s)); dir.create(folder)
    a <- fl_dv_alpha_fit(ms, cf, s)
    a$results$review_primary <- s == "core_adjusted"
    a$results$covariates_added <- paste(z$covariates, collapse = ";")
    a$results$status <- "ESTIMATED"
    a$interaction$review_primary <- s == "core_adjusted"
    alpha[[s]] <- a$results; hetero[[s]] <- a$interaction; diagnostics[[s]] <- a$diagnostic
    saveRDS(a, file.path(folder, "alpha_model.rds"), version = 3)
    # Independently reconstruct the Gaussian independence-GEE sandwich from X.
    X <- model.matrix(a$fit); err <- as.numeric(a$metadata$shannon - X %*% coef(a$fit))
    bread <- solve(crossprod(X)); scores <- rowsum(X * err, a$metadata$family_id)
    V <- bread %*% crossprod(scores) %*% bread
    v_error <- max(abs(V - vcov(a$fit)))
    if (!is.finite(v_error) || v_error > 1e-7) stop("Independent GEE sandwich mismatch: ", s)
    diagnostics[[s]]$independent_sandwich_max_error <- v_error
    if (identical(s, "core_adjusted")) main <- a
    rm(a, X, bread, scores, V); gc(FALSE)
    b <- fl_dv_beta_fit(ms, as.dist(fullD[j, j, drop = FALSE]), cf, folder)
    b$results$scenario <- s; b$results$review_primary <- s == "core_adjusted"
    b$results$covariates_added <- paste(z$covariates, collapse = ";")
    b$results$status <- "ESTIMATED"
    b$exchangeability$scenario <- s
    beta[[s]] <- b$results; exchange[[s]] <- b$exchangeability
    fl_write_csv(b$null_distribution, file.path(audit, paste0("review_beta_null_", s, ".csv")))
    for (g in c("G0", "G1", "G2")) {
      k <- ms$generation == g
      flow[[paste(s, g)]] <- data.frame(scenario = s, generation = g, n = sum(k),
        cases = sum(ms$fatty_liver[k]), minors = sum(ms$age[k] < 18),
        families = length(unique(ms$family_id[k])), excluded_from_fixed_generation = sum(m$generation == g) - sum(k),
        label = z$label, covariates_added = paste(z$covariates, collapse = ";"))
    }
    rm(b); gc(FALSE)
  }
  ans <- list(Shannon_models = do.call(rbind, alpha), generation_heterogeneity = do.call(rbind, hetero),
    composition_association = do.call(rbind, beta), exchangeability = do.call(rbind, exchange),
    sample_counts = do.call(rbind, flow), model_diagnostics = do.call(rbind, diagnostics))
  for (nm in names(ans)) fl_write_csv(ans[[nm]], file.path(if (nm == "model_diagnostics") audit else out,
    paste0("review_", nm, ".csv")))
  fl_dv_figures(main$metadata, main$results, main$fit, out, audit)
  ans
}

fl_dv_review_report <- function(review, out, covariate_version) {
  alpha <- review$Shannon_models; beta <- review$composition_association
  main <- alpha[alpha$review_primary, ]; bm <- beta[beta$review_primary, ]
  sample_total <- function(s, v) sum(review$sample_counts[review$sample_counts$scenario == s, v])
  scope_note <- paste0("基线仍为319个既定三代家庭的1275人。核心完整样本", sample_total("core_adjusted", "n"),
    "人（", sample_total("core_adjusted", "cases"), "例），保留", sample_total("core_adjusted", "minors"),
    "名未成年人；", 1275 - sample_total("core_adjusted", "n"), "人因核心协变量缺失退出，不重新筛三代家系。")
  note <- c("# 08 全面协变量审阅后的结果", "", paste0("协变量数据版本：", covariate_version), "",
    scope_note,
    "review_*文件为本次审阅后结果；Table2/Table3旧模型仅供历史基础参照。图2/图3已改用本次核心完整样本及其调整结果。",
    "核心：各代3df年龄样条、性别、建库批次、log10输入矩阵总量；另加入调查地区、城乡、近期抗生素、便秘、稀便、乳酸菌制剂、采便季节正余弦、采集至记录冷冻间隔、粪便与超声时间差。",
    "核心模型纳入排便/采样因素是条件关联敏感性选择：排便也可能位于疾病作用路径上；不能把它们一概当作已证实的混杂因素，也不能把所有条件效应解释为总因果效应。",
    paste0("core_base与core_adjusted同样本；成人四模型同样本", sample_total("adult_core", "n"),
      "人；营养替代两模型同样本", sample_total("nutrient_base", "n"), "人。变化同时报告效应量和区间，不按P选择模型。"),
    "成人生活方式：吸烟、饮酒、饮食评分、成人PAI活动指数及能量；再加入成人BMI；最后加入糖尿病及降糖/降脂/降压/抗风湿药物。药物类别不能等同特异metformin/statin/PPI。",
    "营养替代以纤维、糖和能量替代综合饮食评分；DNA浓度单独敏感性；排除近期抗生素/胃肠炎为限制样本敏感性。BMI、糖尿病和药物可能涉及共同病理或疾病后果，因此这些是条件敏感性，不强称消除全部混杂。",
    "儿童未被排除出核心分析；儿童与成人活动量表不同，未拼接。未记录的吸烟/饮酒不填0。无可靠教育收入、特异PPI/metformin或标准Bristol评分，未臆造代理。",
    "Shannon：Gaussian identity GEE，家系sandwich；同一模型给出3代系数及2df异质性。独立重建每个模型的sandwich核验数值。",
    "组成：完整SGB Bray–Curtis、Lingoes校正partial dbRDA，Condition家系和同方案协变量；9999次家系内及同家系同代残差置换。每个方案前10次结果与vegan逐项核验。",
    "beta R2_total_corrected相对于Lingoes修正总惯量；R2_partial相对于条件化后的惯量。它们不是同一解释比例。家系内条件组成关联与GEE人群边际均值差为不同目标。",
    "所有P均为原始P：GEE双侧，dbRDA无方向F上尾；无多重比较校正。可交换性假设不能由调整变量或小P证明，须结合exchangeability表；尤其同家系同代G1无法交换。",
    "家系标签未核验直接亲子或同住关系；相对丰度不能证明菌株传播；观察性横断面分析不支持因果/中介结论。",
    "方法依据及适配见config/confounder_review_plan_20261019.md与scripts/08_microbiome_steatosis_diversity_README.md。")
  writeLines(enc2utf8(note), file.path(out, "results_zh.md"), useBytes = TRUE)
  page <- paste0("<!doctype html><html lang='zh-CN'><meta charset='utf-8'><title>08 全面协变量审阅</title>",
    "<style>body{font:16px/1.7 Arial,'Microsoft YaHei',sans-serif;max-width:1300px;margin:25px auto;padding:15px;color:#24384a}table{border-collapse:collapse;font-size:12px}th,td{border:1px solid #cbd5df;padding:6px;white-space:nowrap}th{background:#edf3f7}.scroll{overflow:auto}img{max-width:100%}</style>",
    "<h1>08 脂肪肝与菌群多样性：全面协变量审阅</h1><p>", paste(note[-c(1:4)], collapse = "</p><p>"), "</p>",
    "<h2>核心调整：Shannon均值差及95%CI</h2>", fl_dv_html_table(main),
    "<h2>核心调整：家庭内条件组成关联</h2>", fl_dv_html_table(bm),
    "<h2>全部同样本/敏感性Shannon结果</h2>", fl_dv_html_table(alpha),
    "<h2>代际异质性</h2>", fl_dv_html_table(review$generation_heterogeneity),
    "<h2>全部条件组成关联</h2>", fl_dv_html_table(beta),
    "<h2>各方案分代样本与病例数</h2>", fl_dv_html_table(review$sample_counts),
    "<h2>置换信息</h2>", fl_dv_html_table(review$exchangeability),
    "<img src='Figure2_Shannon_distribution.png'><img src='Figure3_Shannon_adjusted_difference.png'></html>")
  writeLines(enc2utf8(page), file.path(out, "microbiome_steatosis_diversity.html"), useBytes = TRUE)
}
