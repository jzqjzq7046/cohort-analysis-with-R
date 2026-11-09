# Clinical model template: Long et al. Liver Int 2019;39:740-747, section 2.5.
# GEE: Liang & Zeger 1986. Raw two-sided P values, per user decision 2026-10-19.
# See the accompanying README for explicit departures from the published design.

fl_cg_require_packages <- function() {
  packages <- c("geepack", "detectseparation")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Required R packages missing: ", paste(absent, collapse = ", "),
    ". Install explicitly before running; this script never installs packages automatically.")
  invisible(packages)
}

fl_cg_code_paths <- function(cfg) file.path(cfg$project_root, c(
  "scripts/06_cross_generation_clinical.R", "R/cross_generation_clinical.R",
  "config/cross_generation_clinical.R", "config/cross_generation_clinical_references.csv",
  "scripts/06_cross_generation_clinical_README.md", "R/load_processed_microbiome.R",
  "config.R", "R/common.R", "R/entry_helpers.R", "tests/cross_generation_clinical_test.R",
  "R/confounders.R", "R/clinical_confounder_review.R", "tests/clinical_confounder_review_test.R"))

fl_cg_family_key <- function(x) {
  x <- trimws(as.character(x)); x[is.na(x) | x == ""] <- NA_character_
  numeric_label <- !is.na(x) & grepl("^[0-9]+(\\.0+)?$", x)
  x[numeric_label] <- sub("\\.0+$", "", x[numeric_label])
  x[numeric_label] <- sub("^0+(?=[0-9])", "", x[numeric_label], perl = TRUE)
  x
}

fl_cg_validate_metadata <- function(m) {
  # family_id in the processed dataset was cross-checked against tutkno80perhe.
  fields <- c("sample_id", "participant_id", "family_id", "generation", "age",
    "age_group", "sex", "fat_code", "fatty_liver", "BMI18", "alkodpd18")
  if (!all(fields %in% names(m))) stop("Required clinical fields missing: ",
    paste(setdiff(fields, names(m)), collapse = ", "))
  m <- m[, fields, drop = FALSE]
  for (field in c("sample_id", "participant_id")) {
    v <- as.character(m[[field]])
    if (anyNA(v) || any(!nzchar(trimws(v))) || anyDuplicated(v)) stop("Invalid/duplicate ", field)
    m[[field]] <- v
  }
  if (!nrow(m) || anyNA(m$generation) || any(!m$generation %in% c("G0", "G1", "G2")))
    stop("Invalid generation")
  for (field in c("age", "fat_code", "fatty_liver", "BMI18", "alkodpd18")) {
    if (!is.numeric(m[[field]]) || any(!is.na(m[[field]]) & !is.finite(m[[field]])))
      stop("Non-numeric or nonfinite clinical field: ", field)
  }
  if (anyNA(m$fat_code) || any(!m$fat_code %in% 1:3) || anyNA(m$fatty_liver) ||
      any(m$fatty_liver != as.integer(m$fat_code %in% c(2, 3)))) stop("Outcome coding mismatch")
  if (any(!is.na(m$age) & (m$age < 0 | m$age > 120)) ||
      any(!is.na(m$BMI18) & m$BMI18 <= 0) ||
      any(!is.na(m$alkodpd18) & m$alkodpd18 < 0)) stop("Invalid age, BMI or alcohol; inspect source coding")
  m$family_id <- fl_cg_family_key(m$family_id)
  m$sex <- as.character(m$sex)
  if (any(!is.na(m$sex) & !m$sex %in% c("Female", "Male"))) stop("Invalid sex labels")
  # Exact age at ultrasound, not calendar-year age; missing age is not adulthood.
  m$age_group <- ifelse(is.na(m$age), "Age unavailable", ifelse(m$age < 18, "Minor", "Adult"))
  rownames(m) <- NULL
  m
}

fl_cg_make_comparison <- function(m, generations) {
  older <- m[m$generation == generations[1L] & !is.na(m$family_id), , drop = FALSE]
  younger <- m[m$generation == generations[2L], , drop = FALSE]
  blocks <- split(seq_len(nrow(older)), older$family_id)
  senior <- if (length(blocks)) do.call(rbind, lapply(names(blocks), function(f) {
    z <- older[blocks[[f]], , drop = FALSE]
    data.frame(family_id = f, n_older_observed = nrow(z), n_older_cases = sum(z$fatty_liver),
      older_any_case = as.integer(any(z$fatty_liver == 1)), stringsAsFactors = FALSE)
  })) else data.frame(family_id = character(), n_older_observed = integer(),
    n_older_cases = integer(), older_any_case = integer())
  rownames(senior) <- NULL
  k <- match(younger$family_id, senior$family_id)
  known_family <- !is.na(younger$family_id)
  available <- known_family & !is.na(k)
  d <- younger[available, , drop = FALSE]
  k <- k[available]
  for (field in setdiff(names(senior), "family_id")) d[[field]] <- senior[[field]][k]
  d$comparison <- rep(paste(generations, collapse = "-"), nrow(d))
  if (anyDuplicated(d$participant_id) || any(d$n_older_observed < 1L) ||
      any(d$older_any_case != as.integer(d$n_older_cases > 0L))) stop("Comparison construction failure")
  d$sex <- factor(d$sex, levels = c("Female", "Male"))
  rownames(d) <- NULL
  flow <- data.frame(comparison = paste(generations, collapse = "-"),
    stage = c("All younger-generation cohort members", "Younger members without a primary family label",
      "Known label but no eligible observed older-generation member", "Eligible younger individuals"),
    n = c(nrow(younger), sum(!known_family), sum(known_family & is.na(match(younger$family_id,
      senior$family_id))), nrow(d)), stringsAsFactors = FALSE)
  stopifnot(flow$n[1L] == sum(flow$n[2:4]))
  list(data = d, flow = flow)
}

fl_cg_describe <- function(d, comparison) {
  out <- list(); j <- 0L
  for (stratum in c("All", "Adult", "Minor", "Age unavailable")) for (exposure in 0:1) {
    z <- d[d$older_any_case == exposure & (stratum == "All" | d$age_group == stratum), , drop = FALSE]
    j <- j + 1L; a <- z$age[!is.na(z$age)]; b <- z$BMI18[!is.na(z$BMI18)]
    out[[j]] <- data.frame(comparison = comparison, age_stratum = stratum,
      older_any_case = exposure,
      exposure_label = if (exposure) "At least one observed older case" else "No observed older case",
      n = nrow(z), cases = sum(z$fatty_liver), controls = sum(z$fatty_liver == 0),
      prevalence_pct = if (nrow(z)) 100 * mean(z$fatty_liver) else NA_real_,
      families = length(unique(z$family_id)), females_n = sum(z$sex == "Female", na.rm = TRUE),
      sex_missing_n = sum(is.na(z$sex)), age_missing_n = sum(is.na(z$age)),
      age_mean = if (length(a)) mean(a) else NA_real_,
      age_sd = if (length(a) > 1L) sd(a) else NA_real_,
      BMI_nonmissing_n = length(b), BMI_mean = if (length(b)) mean(b) else NA_real_,
      BMI_sd = if (length(b) > 1L) sd(b) else NA_real_, alcohol_missing_n = sum(is.na(z$alkodpd18)),
      mean_observed_older_n = if (nrow(z)) mean(z$n_older_observed) else NA_real_)
  }
  do.call(rbind, out)
}

fl_cg_model_specs <- function() {
  # Same-sample pairs distinguish adjustment from complete-case sample changes.
  # Models beyond M1 are supportive; they do not replace the primary model by P value.
  data.frame(model = c("M0_crude", "M1_primary", "M1_adult_BMI_complete",
      "M2_adult_BMI", "M2_adult_alcohol_complete", "M3_adult_BMI_alcohol"),
    sample_set = c("all", "core", "adult_BMI", "adult_BMI", "adult_alcohol", "adult_alcohol"),
    covariates = c("", "age+sex+n_older_observed", "age+sex+n_older_observed",
      "age+sex+n_older_observed+BMI18", "age+sex+n_older_observed+BMI18",
      "age+sex+n_older_observed+BMI18+alkodpd18"),
    primary = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
    comparison_pair = c("", "", "BMI_adjustment", "BMI_adjustment", "Alcohol_adjustment", "Alcohol_adjustment"),
    references = c("Long2019;LiangZeger1986;Geepack2006;Konis2007",
      "Long2019;LiangZeger1986;Geepack2006;Konis2007",
      rep("Long2019;LiangZeger1986;Geepack2006;Konis2007", 4L)),
    study_adaptation = c("Family-label exposure; ultrasound outcome; saved-cohort observation scope",
      "All ages; observed-older-count adjustment; alcohol deferred to adult extension",
      "Adult complete-BMI comparison sample; observed-older-count adjustment",
      "Adult complete-BMI comparison sample; observed-older-count adjustment",
      "Adult complete-BMI/alcohol comparison sample; observed-older-count adjustment",
      "Adult complete-BMI/alcohol comparison sample; observed-older-count adjustment"), stringsAsFactors = FALSE)
}

fl_cg_select_sets <- function(d) {
  core <- complete.cases(d[, c("age", "sex"), drop = FALSE])
  adult_BMI <- core & d$age_group == "Adult" & !is.na(d$BMI18)
  adult_alcohol <- adult_BMI & !is.na(d$alkodpd18)
  list(all = rep(TRUE, nrow(d)), core = core, adult_BMI = adult_BMI, adult_alcohol = adult_alcohol)
}

fl_cg_fit <- function(d, spec, comparison, cfg) {
  n <- nrow(d); families <- length(unique(d$family_id)); sizes <- table(d$family_id)
  covariates <- if (nzchar(spec$covariates)) strsplit(spec$covariates, "+", fixed = TRUE)[[1L]] else character()
  constant <- covariates[vapply(covariates, function(v) length(unique(d[[v]])) < 2L, logical(1))]
  terms <- c("older_any_case", setdiff(covariates, constant))
  form <- reformulate(terms, response = "fatty_liver")
  row <- data.frame(comparison = comparison, model = spec$model, primary = spec$primary,
    sample_set = spec$sample_set, n = n, cases = sum(d$fatty_liver), controls = sum(d$fatty_liver == 0),
    families = families, singleton_families = sum(sizes == 1L),
    exposed_families = length(unique(d$family_id[d$older_any_case == 1])),
    unexposed_families = length(unique(d$family_id[d$older_any_case == 0])),
    exposed_n = sum(d$older_any_case == 1), unexposed_n = sum(d$older_any_case == 0),
    minors_n = sum(d$age_group == "Minor"),
    formula = paste(deparse(form), collapse = ""), dropped_constant_terms = paste(constant, collapse = ";"),
    design_columns = NA_integer_, design_rank = NA_integer_, events_per_parameter = NA_real_,
    status = "NOT_ESTIMATED", reason = "", warnings = "", gee_error = NA_integer_,
    initial_glm_converged = NA,
    working_correlation = cfg$cg_working_correlation,
    log_OR = NA_real_, robust_SE = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
    p_value = NA_real_, p_value_adjustment = "none (raw two-sided P)",
    confidence_level = 1 - cfg$cg_alpha, CI_multiplicity_adjusted = FALSE,
    references = spec$references, study_adaptation = spec$study_adaptation, stringsAsFactors = FALSE)
  cells <- expand.grid(older_any_case = 0:1, younger_case = 0:1, KEEP.OUT.ATTRS = FALSE)
  cells$n <- vapply(seq_len(nrow(cells)), function(i) sum(d$older_any_case == cells$older_any_case[i] &
    d$fatty_liver == cells$younger_case[i]), integer(1))
  cells$comparison <- comparison; cells$model <- spec$model
  finish <- function(status, reason = "", warnings = character()) {
    row$status <- status; row$reason <- reason; row$warnings <- paste(unique(warnings), collapse = " | ")
    list(result = row, cells = cells)
  }
  if (!n || families < 2L) return(finish("INSUFFICIENT_DATA", "Fewer than two observed families"))
  if (length(unique(d$fatty_liver)) < 2L || length(unique(d$older_any_case)) < 2L)
    return(finish("NO_VARIATION", "Outcome or exposure has one level"))
  d <- droplevels(d[order(d$family_id, d$participant_id), , drop = FALSE])
  d$cluster_id <- match(d$family_id, unique(d$family_id))
  # Sorting is mandatory for geeglm: each family's records must be contiguous.
  mm <- model.matrix(form, d)
  row$design_columns <- ncol(mm); row$design_rank <- qr(mm)$rank
  row$events_per_parameter <- min(row$cases, row$controls) / max(1L, ncol(mm) - 1L)
  if (nrow(mm) <= ncol(mm) || qr(mm)$rank != ncol(mm))
    return(finish("RANK_DEFICIENT", "Insufficient residual degrees of freedom or aliased covariates"))
  warnings <- character()
  collect <- function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
  # Separation is a structural diagnostic. No independent Firth model is silently
  # substituted for a failed family-correlated GEE.
  sep <- tryCatch(withCallingHandlers(glm(form, data = d, family = binomial("logit"),
    method = detectseparation::detect_separation), warning = collect), error = function(e) e)
  if (inherits(sep, "error")) return(finish("SEPARATION_CHECK_FAILED", conditionMessage(sep), warnings))
  if (any(is.infinite(coef(sep))))
    return(finish("SEPARATION", "Complete/quasi-complete separation: no ordinary logistic GEE inference", warnings))
  if (anyNA(coef(sep))) return(finish("SEPARATION_CHECK_FAILED", "Undetermined separation coefficients", warnings))
  initial <- tryCatch(withCallingHandlers(glm(form, data = d, family = binomial("logit"),
    control = glm.control(maxit = cfg$cg_max_iterations)), warning = collect), error = function(e) e)
  if (inherits(initial, "error")) return(finish("FIT_FAILED", conditionMessage(initial), warnings))
  row$initial_glm_converged <- isTRUE(initial$converged)
  if (!row$initial_glm_converged || any(!is.finite(coef(initial))))
    return(finish("NOT_CONVERGED", "Initial logistic fit did not converge", warnings))
  fit <- tryCatch(withCallingHandlers(geepack::geeglm(form, data = d,
    id = cluster_id, family = binomial("logit"), corstr = cfg$cg_working_correlation,
    # Default scale.value is 1. Avoid forwarding an explicit scale.value/start
    # through model.frame in geepack releases that do not remove these extras.
    std.err = "san.se", scale.fix = TRUE,
    control = geepack::geese.control(maxit = cfg$cg_max_iterations)), warning = collect), error = function(e) e)
  if (inherits(fit, "error")) return(finish("FIT_FAILED", conditionMessage(fit), warnings))
  gee_error <- fit$geese$error
  if (length(gee_error) != 1L || is.na(gee_error))
    return(finish("NOT_CONVERGED", "GEE convergence status unavailable", warnings))
  row$gee_error <- as.integer(gee_error)
  # geeglm intentionally removes the glm $converged field; use geese$error here.
  if (row$gee_error != 0L) return(finish("NOT_CONVERGED", "GEE did not converge", warnings))
  if (any(grepl("did not converge", warnings, fixed = TRUE)))
    return(finish("NOT_CONVERGED", "A fitting stage reported non-convergence", warnings))
  beta <- coef(fit); v <- fit$geese$vbeta
  if (!is.matrix(v) || nrow(v) != length(beta) || any(!is.finite(beta)) || any(!is.finite(v)))
    return(finish("INVALID_COVARIANCE", "Nonfinite coefficients or robust sandwich covariance", warnings))
  i <- match("older_any_case", names(beta))
  if (is.na(i) || v[i, i] <= 0)
    return(finish("INVALID_COVARIANCE", "Non-positive or unavailable exposure variance", warnings))
  se <- sqrt(v[i, i])
  if (!is.finite(se) || se <= 0) return(finish("INVALID_COVARIANCE", "Non-positive exposure variance", warnings))
  z <- qnorm(1 - cfg$cg_alpha / 2)
  limits <- exp(c(beta[i], beta[i] - z * se, beta[i] + z * se))
  if (any(!is.finite(limits)) || any(limits <= 0))
    return(finish("NUMERIC_OVERFLOW", "Exponentiated effect or CI not finite and positive", warnings))
  row$log_OR <- unname(beta[i]); row$robust_SE <- se
  row$OR <- limits[1L]; row$CI_low <- limits[2L]; row$CI_high <- limits[3L]
  row$p_value <- 2 * pnorm(-abs(beta[i] / se))
  if (families < 50L) warnings <- c(warnings, "Few clusters: asymptotic sandwich inference may be unreliable")
  if (min(row$exposed_families, row$unexposed_families) < 10L)
    warnings <- c(warnings, "Sparse family-level exposure groups; interpret asymptotic inference cautiously")
  if (any(cells$n == 0L)) warnings <- c(warnings, "Zero exposure/outcome cell")
  if (row$events_per_parameter < 10) warnings <- c(warnings,
    "Low events per parameter (diagnostic flag only; not an automatic variable-selection rule)")
  finish(if (length(warnings)) "ESTIMATED_WITH_WARNINGS" else "ESTIMATED", warnings = warnings)
}

fl_cg_prepare_metadata <- function(m) {
  fl_three_generation_scope(fl_cg_validate_metadata(m))
}

fl_cg_compute <- function(m, cfg) {
  scope <- fl_cg_prepare_metadata(m)
  m <- scope$cohort; specs <- fl_cg_model_specs()
  flows <- descriptions <- coverage <- missing <- results <- cells <- pair_audits <- list()
  private <- list(); counter <- 0L
  for (g in cfg$cg_pairs) {
    comparison <- paste(g, collapse = "-"); built <- fl_cg_make_comparison(m, g); d <- built$data
    flows[[comparison]] <- built$flow; descriptions[[comparison]] <- fl_cg_describe(d, comparison)
    nvalues <- sort(unique(d$n_older_observed))
    coverage[[comparison]] <- if (length(nvalues)) do.call(rbind, lapply(nvalues, function(k) {
      z <- d[d$n_older_observed == k, , drop = FALSE]
      data.frame(comparison = comparison, n_older_observed = k, n = nrow(z),
        families = length(unique(z$family_id)), older_any_case_n = sum(z$older_any_case),
        younger_cases = sum(z$fatty_liver))
    })) else data.frame(comparison = character(), n_older_observed = integer(), n = integer(),
      families = integer(), older_any_case_n = integer(), younger_cases = integer())
    sets <- fl_cg_select_sets(d)
    for (v in c("age", "sex", "BMI18", "alkodpd18")) {
      for (stratum in c("All", "Adult", "Minor", "Age unavailable")) {
        z <- d[stratum == "All" | d$age_group == stratum, , drop = FALSE]
        missing[[paste(comparison, v, stratum)]] <- data.frame(comparison = comparison,
          age_stratum = stratum, field = v, n = nrow(z), missing_n = sum(is.na(z[[v]])),
          cases_among_missing = sum(z$fatty_liver[is.na(z[[v]])]))
      }
    }
    for (s in seq_len(nrow(specs))) {
      spec <- specs[s, , drop = FALSE]; z <- d[sets[[spec$sample_set]], , drop = FALSE]
      f <- fl_cg_fit(z, spec, comparison, cfg); counter <- counter + 1L
      f$result$eligible_comparison_n <- nrow(d)
      f$result$excluded_from_model_n <- nrow(d) - nrow(z)
      f$result$excluded_cases_n <- sum(d$fatty_liver) - sum(z$fatty_liver)
      results[[counter]] <- f$result; cells[[counter]] <- f$cells
    }
    for (pair in unique(specs$comparison_pair[nzchar(specs$comparison_pair)])) {
      ids <- which(specs$comparison_pair == pair)
      rows_a <- d$participant_id[sets[[specs$sample_set[ids[1L]]]]]
      rows_b <- d$participant_id[sets[[specs$sample_set[ids[2L]]]]]
      stopifnot(identical(rows_a, rows_b))
      pair_audits[[paste(comparison, pair)]] <- data.frame(comparison = comparison,
        adjustment = pair, first_model = specs$model[ids[1L]], second_model = specs$model[ids[2L]],
        n = length(rows_a), same_participants = TRUE)
    }
    if (isTRUE(cfg$cg_export_private)) private[[comparison]] <- d
  }
  result <- do.call(rbind, results)
  p <- which(result$primary)
  if (length(p) != 3L) stop("Expected exactly three prespecified primary comparisons")
  list(family_eligibility = scope$summary,
    cohort_summary = data.frame(n = nrow(m), families = scope$contract$families,
      cases = sum(m$fatty_liver), rule = scope$contract$rule),
    flow = do.call(rbind, flows), descriptive = do.call(rbind, descriptions),
    older_coverage = do.call(rbind, coverage), missingness = do.call(rbind, missing),
    models = result, model_cells = do.call(rbind, cells),
    same_sample_audit = do.call(rbind, pair_audits), model_plan = specs, private = private)
}

fl_cg_html <- function(results, path) {
  escape <- function(x) {
    x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE)
  }
  table_html <- function(d) {
    cols <- lapply(d, function(x) if (is.numeric(x)) ifelse(is.na(x), "", format(signif(x, 4),
      trim = TRUE)) else ifelse(is.na(x), "", as.character(x)))
    rows <- vapply(seq_len(nrow(d)), function(i) paste0("<tr>", paste0("<td>",
      escape(vapply(cols, `[`, character(1), i)), "</td>", collapse = ""), "</tr>"), character(1))
    paste0("<div class='scroll'><table><tr>", paste0("<th>", escape(names(d)), "</th>", collapse = ""),
      "</tr>", paste(rows, collapse = "\n"), "</table></div>")
  }
  notes <- c("<h1>Cross-generation hepatic steatosis / 跨代肝脂肪变关联</h1>",
    "<p>Each younger individual contributes once per generation comparison. All eligible ages are retained in the primary analysis.</p>",
    "<p>所有比较仅限G0/G1/G2每代至少一名基线合格成员的主家系；家系集合在模型协变量限制之前确定。All analyses are restricted to families represented in all three generations.</p>",
    "<p>仅使用 tutkno80perhe；直接亲子或祖孙关系尚未核验。未观察到年长成员者不归为未暴露；对照仅表示已观察成员中无病例。</p>",
    "<p>Legacy reference: logistic GEE with family-cluster sandwich SE; age, sex and number of observed older members. Revised core models additionally include observed older age/sex composition and survey centre/urbanicity. ORs describe cross-sectional odds, not incidence risk or causal transmission.</p>",
    "<p>按2026-10-19用户决定，所有比较报告原始双侧P值，不进行多重比较校正；95%CI为逐项区间。成人BMI和饮酒扩展使用配对的同一样本模型；未成年人仅不进入成人扩展。</p>",
    "<p>All results, including failed or non-estimable models, are retained below. No clinical P-value gates later microbiome work.</p>")
  sections <- list("Family eligibility / 三代家系筛选" = results$family_eligibility,
    "Restricted cohort / 限定队列" = results$cohort_summary,
    "Sample accounting / 样本流程" = results$flow,
    "Case distribution / 病例分布" = results$descriptive,
    "Legacy base GEE / 原基础调整参照" = results$models[results$models$primary, , drop = FALSE],
    "All planned models / 全部预设模型" = results$models,
    "Missingness / 缺失情况" = results$missingness,
    "Observed older-member coverage / 年长成员观察人数" = results$older_coverage,
    "Same-sample checks / 同样本核对" = results$same_sample_audit)
  if (!is.null(results$review_models)) {
    sections <- c(list("Revised primary clinical association / 审阅后的全年龄主要关联" =
      results$review_models[results$review_models$review_primary, , drop = FALSE],
      "All prespecified covariate review models / 全部调整扩展" = results$review_models,
      "Review same-sample checks / 扩展同样本核对" = results$review_same_sample_audit), sections)
    notes <- c(notes, "<p><b>本次审阅优先阅读review_models：全年龄模型增加年长成员年龄和女性比例及调查地区/城乡；成人扩展增加饮酒、吸烟、饮食、运动及能量，再加入BMI和药物。原18模型只作为参照。所有未估计状态和完整样本损失均保留。</b></p>")
  }
  lines <- c("<!doctype html><html><head><meta charset='utf-8'><title>Cross-generation clinical association</title>",
    "<style>body{font:14px Arial,sans-serif;margin:28px;line-height:1.6;color:#243746}table{border-collapse:collapse}th,td{border-bottom:1px solid #ddd;padding:6px;text-align:right;white-space:nowrap}th{background:#eef2f5}.scroll{overflow-x:auto;margin-bottom:24px}</style></head><body>", notes)
  for (label in names(sections)) lines <- c(lines, paste0("<h2>", escape(label), "</h2>"), table_html(sections[[label]]))
  con <- file(path, "w", encoding = "UTF-8"); on.exit(close(con))
  writeLines(c(lines, "</body></html>"), con)
}

fl_cg_run <- function(cfg) {
  # No data are loaded until this explicit entry point is called.
  source_paths <- fl_cg_code_paths(cfg)
  if (any(!file.exists(source_paths))) stop("Required code/documentation file missing")
  context <- fl_output_start_run(cfg); cfg <- context$cfg; run_dir <- context$run_dir
  output <- file.path(run_dir, "clinical"); audit <- file.path(run_dir, "audit")
  dir.create(output); dir.create(audit)
  status <- file.path(run_dir, "run_status.txt"); writeLines("RUNNING", status)
  cat("Starting cross-generation clinical analysis\nOutput: ", run_dir, "\n", sep = "")
  tryCatch({
    package <- fl_load_processed_microbiome(project_root = cfg$project_root,
      version = cfg$cg_data_version, verify_checksum = TRUE, require_three_generation = FALSE)
    metadata <- package$metadata; input_version <- package$version
    rm(package) # Abundance matrices are never used in this clinical module.
    results <- fl_cg_compute(metadata, cfg)
    reviewed <- fl_cv_load(cfg, metadata)
    expanded <- fl_cg_review_compute(reviewed$metadata, cfg)
    results <- c(results, expanded)
    fl_write_csv(data.frame(covariate_version = reviewed$version), file.path(audit, "covariate_version.csv"))
    for (name in setdiff(names(results), "private"))
      fl_write_csv(results[[name]], file.path(output, paste0(name, ".csv")))
    fl_cg_review_note(results, file.path(output, "confounder_review_zh.md"))
    fl_cg_html(results, file.path(output, "cross_generation_clinical.html"))
    if (isTRUE(cfg$cg_export_private)) {
      dir.create(file.path(run_dir, "private"))
      saveRDS(results$private, file.path(run_dir, "private", "younger_individual_datasets.rds"))
    }
    input_file <- file.path(cfg$processed_microbiome_root, input_version, "microbiome_processed.rds")
    files <- c(input_file, reviewed$source_files, source_paths)
    fl_write_csv(data.frame(path = files, md5 = unname(tools::md5sum(files))),
      file.path(audit, "input_code_fingerprints.csv"))
    fl_write_csv(data.frame(data_version = input_version, input_n = nrow(metadata),
      cohort_n = results$cohort_summary$n, families = results$cohort_summary$families,
      family_eligibility_rule = "G0_G1_G2_each_at_least_one; before covariate restrictions",
      input_scope = "Three-generation primary families within the saved ultrasound-eligible stool cohort",
      family_definition = "tutkno80perhe only; direct pedigree unverified",
      exposure_definition = "At least one eligible observed older member is a case",
      prevalence_CI = "Not calculated; no independent-binomial assumption",
      GEE_working_correlation = cfg$cg_working_correlation,
      variance = "Family-cluster robust sandwich; asymptotic normal Wald intervals",
      multiplicity = "None; raw two-sided P values for all models, per user decision 2026-10-19; pointwise 95% CI",
      optional_private_export = isTRUE(cfg$cg_export_private)), file.path(audit, "analysis_contract.csv"))
    if (!file.copy(file.path(cfg$project_root, "config", "cross_generation_clinical_references.csv"),
      file.path(audit, "references.csv"))) stop("Cannot copy references")
    capture.output(sessionInfo(), file = file.path(audit, "sessionInfo.txt"))
    issues <- sum(!grepl("^ESTIMATED", results$models$status)) +
      sum(!grepl("^ESTIMATED", results$review_models$status))
    writeLines(c("COMPLETED", paste("Models without an estimate:", issues),
      "Read model status/warnings before interpreting any estimate."), status)
    final_dir <- fl_output_complete_run(context)
    cat("Completed. Models without an estimate: ", issues, "\nReport: ",
      file.path(final_dir, "clinical", "cross_generation_clinical.html"), "\n", sep = "")
    invisible(results[setdiff(names(results), "private")])
  }, error = function(e) {
    fl_output_fail_run(context, conditionMessage(e)); stop(e)
  })
}
