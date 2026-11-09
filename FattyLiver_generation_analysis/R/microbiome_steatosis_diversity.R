# 08: prespecified community-level steatosis associations; no taxa screening.
fl_dv_code_paths <- function(cfg) file.path(cfg$project_root, c("config.R", "config/microbiome_steatosis_diversity.R",
  "R/common.R", "R/entry_helpers.R", "R/load_processed_microbiome.R", "R/microbiome_steatosis_diversity.R",
  "scripts/08_microbiome_steatosis_diversity.R", "scripts/08_microbiome_steatosis_diversity_README.md",
  "config/microbiome_steatosis_diversity_references.csv", "scripts/dbrda_permutation_kernel.py",
  "tests/microbiome_steatosis_diversity_test.R", "R/confounders.R", "R/diversity_confounder_review.R",
  "config/confounder_review_plan_20261019.md", "tests/diversity_confounder_review_test.R"))

fl_dv_prepare <- function(m, relative, counts) {
  required <- c("sample_id", "family_id", "generation", "age", "sex", "fatty_liver", "prep_batch", "matrix_total")
  if (!all(required %in% names(m)) || anyDuplicated(m$sample_id) ||
      !identical(colnames(relative), m$sample_id) || !identical(dimnames(relative), dimnames(counts)) ||
      any(!is.finite(relative)) || any(relative < 0) || max(abs(colSums(relative) - 1)) > 1e-10 ||
      any(!is.finite(counts)) || any(counts < 0) || any(colSums(counts) <= 0)) stop("Invalid input or alignment")
  scope <- fl_three_generation_scope(m)
  if (nrow(scope$cohort) != nrow(m)) stop("Input contains ineligible families")
  if (anyNA(m[required]) || any(!m$fatty_liver %in% c(0, 1)) || any(m$matrix_total <= 0) ||
      max(abs(colSums(counts) - m$matrix_total)) > 1e-8) stop("Invalid core metadata")
  shannon <- -colSums(ifelse(relative > 0, relative * log(pmax(relative, .Machine$double.xmin)), 0))
  if ("shannon" %in% names(m) && max(abs(m$shannon - shannon)) > 1e-10) stop("Stored Shannon mismatch")
  m$shannon <- shannon; m$observed_sgbs <- colSums(counts > 0)
  m$generation <- factor(m$generation, levels = c("G0", "G1", "G2"))
  m$sex <- factor(m$sex); m$prep_batch <- factor(m$prep_batch)
  m$log_total <- log10(m$matrix_total)
  m
}

fl_dv_design <- function(m, cfg) {
  m <- droplevels(m)
  if (any(table(m$generation) < cfg$dv_age_df + 2L)) stop("Insufficient ages for prespecified spline")
  age_names <- character(); knots <- list()
  for (g in levels(m$generation)) {
    j <- which(m$generation == g)
    basis <- splines::ns(m$age[j], df = cfg$dv_age_df)
    for (k in seq_len(ncol(basis))) {
      nm <- paste0("age_", g, "_", k); m[[nm]] <- 0; m[j, nm] <- basis[, k]
      age_names <- c(age_names, nm)
    }
    m[[paste0("case_", g)]] <- as.numeric(m$generation == g) * m$fatty_liver
    knots[[g]] <- data.frame(generation = g, interior_knots = paste(signif(attr(basis, "knots"), 12), collapse = ";"),
      boundary_knots = paste(signif(attr(basis, "Boundary.knots"), 12), collapse = ";"), df = cfg$dv_age_df)
  }
  nuisance <- unique(c("generation", age_names, "sex", "prep_batch", "log_total", cfg$dv_extra_covariates))
  # A restricted sensitivity sample can make e.g. antibiotic use constant.
  constant <- nuisance[vapply(m[nuisance], function(v) length(unique(v)) < 2L, logical(1))]
  list(metadata = m, nuisance = setdiff(nuisance, constant), constant_nuisance = constant,
       age_knots = do.call(rbind, knots))
}

fl_dv_alpha_fit <- function(m, cfg, scenario, bmi = FALSE) {
  m <- m[order(m$family_id, m$sample_id), , drop = FALSE]
  z <- fl_dv_design(m, cfg); m <- z$metadata
  case_terms <- paste0("case_", c("G0", "G1", "G2"))
  terms <- c(z$nuisance, case_terms, if (bmi) "BMI18")
  form <- stats::reformulate(terms, response = "shannon")
  requested_formula <- paste(deparse(form), collapse = " "); aliases <- character()
  if (isTRUE(cfg$dv_review_mode)) {
    # Retain the entire nuisance column space without forcing redundant dummies.
    # Keep all three outcome contrasts and never prune them by significance.
    Z <- model.matrix(reformulate(c(z$nuisance, if (bmi) "BMI18")), m)
    qZ <- qr(Z); keep <- qZ$pivot[seq_len(qZ$rank)]
    aliases <- setdiff(colnames(Z), colnames(Z)[keep])
    Z <- Z[, keep, drop = FALSE]
    nm <- sprintf("nuisance_%03d", seq_len(ncol(Z)))
    for (i in seq_along(nm)) m[[nm[i]]] <- Z[, i]
    form <- reformulate(c(nm, case_terms), response = "shannon", intercept = FALSE)
  }
  X <- stats::model.matrix(form, m); rank <- qr(X)$rank
  if (rank != ncol(X) || nrow(X) <= rank) stop("Alpha design is rank deficient: ", scenario)
  m$cluster_id <- as.integer(factor(m$family_id))
  captured <- character()
  fit <- withCallingHandlers(geepack::geeglm(form, id = cluster_id, data = m,
    family = gaussian("identity"), corstr = "independence", std.err = "san.se"),
    warning = function(w) { captured <<- c(captured, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (fit$geese$error != 0L || any(!is.finite(coef(fit)))) stop("GEE failed: ", scenario)
  b <- coef(fit); V <- stats::vcov(fit); ix <- match(case_terms, names(b))
  se <- sqrt(diag(V)[ix]); est <- b[ix]
  if (any(!is.finite(se)) || any(se <= 0)) stop("Invalid robust variance")
  results <- do.call(rbind, lapply(seq_along(ix), function(i) {
    g <- c("G0", "G1", "G2")[i]; j <- m$generation == g
    data.frame(scenario = scenario, generation = g, n = sum(j), cases = sum(m$fatty_liver[j]),
      families = length(unique(m$family_id[j])), estimate = unname(est[i]), robust_SE = se[i],
      CI_low = unname(est[i] - qnorm(.975) * se[i]), CI_high = unname(est[i] + qnorm(.975) * se[i]),
      z = unname(est[i] / se[i]), p_raw = 2 * pnorm(-abs(unname(est[i] / se[i]))),
      effect = "adjusted mean Shannon difference: steatosis minus no steatosis", p_adjustment = "none")
  }))
  H <- rbind(c(1, -1, 0), c(0, 1, -1)); delta <- H %*% b[ix]; W <- H %*% V[ix, ix] %*% t(H)
  chi <- as.numeric(crossprod(delta, solve(W, delta)))
  interaction <- data.frame(scenario = scenario, hypothesis = "G0=G1=G2 steatosis-Shannon coefficients",
    Wald_chisq = chi, df = 2L, p_raw = pchisq(chi, 2, lower.tail = FALSE), p_adjustment = "none")
  diagnostic <- data.frame(scenario = scenario, n = nrow(m), families = length(unique(m$family_id)),
    cases = sum(m$fatty_liver), minors = sum(m$age < 18), design_columns = ncol(X), design_rank = rank,
    residual_df = nrow(m) - rank, working_correlation = "independence", variance = "family sandwich",
    gee_error = fit$geese$error, max_abs_standardized_residual = max(abs(residuals(fit) / sd(residuals(fit)))),
    formula = requested_formula, constant_nuisance = paste(z$constant_nuisance, collapse = ";"),
    redundant_nuisance_columns = paste(aliases, collapse = ";"), warnings = paste(unique(captured), collapse = "; "))
  residual_summary <- do.call(rbind, lapply(split(seq_len(nrow(m)), interaction(m$generation, m$fatty_liver, drop = TRUE)), function(j)
    data.frame(scenario = scenario, generation = as.character(m$generation[j[1]]), fatty_liver = m$fatty_liver[j[1]],
      n = length(j), residual_mean = mean(residuals(fit)[j]), residual_SD = sd(residuals(fit)[j]))))
  list(results = results, interaction = interaction, diagnostic = diagnostic, fit = fit, metadata = m,
       residual_summary = residual_summary, age_knots = cbind(scenario = scenario, z$age_knots))
}

fl_dv_permutations <- function(block, B, seed) {
  set.seed(seed)
  P <- permute::shuffleSet(length(block), nset = B,
    control = permute::how(nperm = B, blocks = factor(block)))
  if (nrow(P) != B || anyNA(P)) stop("Unexpected permutation count")
  # Check every row without allocating another B x n character matrix.
  for (i in seq_len(nrow(P))) if (any(block[P[i, ]] != block)) stop("Permutation crossed a block")
  P
}

fl_dv_fast_setup <- function(fit, m, Z) {
  # Same reduced residual Gower matrix used by vegan::permutest.cca.
  E <- getFromNamespace("ordiYbar", "vegan")(fit, "partial")
  fam <- as.character(m$family_id)
  C <- Z[, !grepl("^family_id", colnames(Z)) & colnames(Z) != "(Intercept)", drop = FALSE]
  group_sum <- rowsum(C, fam); size <- table(fam)
  group_mean <- group_sum / as.numeric(size[rownames(group_sum)])
  C <- C - group_mean[match(fam, rownames(group_mean)), , drop = FALSE]
  C <- C[, colSums(C^2) > 1e-20, drop = FALSE]
  qC <- qr(C); Qc <- qr.Q(qC)[, seq_len(qC$rank), drop = FALSE]
  x <- qr.resid(qr(Z), m$fatty_liver); qx <- x / sqrt(sum(x^2))
  Q <- cbind(qx, Qc)
  if (max(abs(crossprod(Q) - diag(ncol(Q)))) > 1e-7 ||
      max(abs(rowsum(E, fam))) > 1e-7) stop("Fast kernel projection identities failed")
  list(E = E, Q = Q, nuisance_rank = ncol(Qc),
    residual_df = nrow(m) - fit$CCA$QR$rank - 1L)
}

fl_dv_fast_null <- function(setup, P, cfg, private_dir, scheme) {
  prefix <- file.path(private_dir, paste0("work_kernel_", scheme))
  binary <- function(value, suffix) {
    con <- file(paste0(prefix, suffix), "wb"); on.exit(close(con))
    writeBin(value, con, endian = "little")
  }
  binary(as.double(setup$E), "_E.bin"); binary(as.double(setup$Q), "_Q.bin")
  binary(as.integer(P), "_P.bin")
  args <- c(shQuote(file.path(cfg$project_root, "scripts/dbrda_permutation_kernel.py")),
    "--prefix", shQuote(prefix), "--n", nrow(setup$E), "--k", ncol(setup$Q),
    "--permutations", nrow(P), "--df", setup$residual_df)
  log <- paste0(prefix, ".log")
  status <- system2(cfg$python, args = args, stdout = log, stderr = log)
  if (status != 0L) stop("Optimized kernel failed; see ", log)
  values <- utils::read.csv(paste0(prefix, "_F.csv"))$F_null
  if (length(values) != nrow(P) || any(!is.finite(values))) stop("Incomplete kernel output")
  cat(paste(readLines(log, warn = FALSE), collapse = "\n"), "\n"); flush.console()
  values
}

fl_dv_beta_fit <- function(m, D, cfg, private_dir) {
  z <- fl_dv_design(m, cfg); m <- z$metadata
  nuisance <- c("family_id", z$nuisance)
  m$family_id <- factor(m$family_id)
  Z <- stats::model.matrix(stats::reformulate(nuisance), m)
  x_residual <- qr.resid(qr(Z), m$fatty_liver)
  if (sum(x_residual^2) < 1e-8) stop("Steatosis has no residual variation after conditioning")
  form <- as.formula(paste("D ~ fatty_liver + Condition(", paste(nuisance, collapse = " + "), ")"))
  fit <- vegan::dbrda(form, data = m, add = "lingoes")
  if (is.null(fit$CCA) || fit$CCA$qrank != 1L || fit$CA$rank < 1L) stop("Beta model not estimable")
  saveRDS(list(fit = fit, metadata = m, condition_matrix = Z, exposure_residual = x_residual),
    file.path(private_dir, "beta_model.rds"), version = 3)
  setup <- fl_dv_fast_setup(fit, m, Z)
  results <- nulls <- blocks_audit <- list()
  schemes <- list(family = as.character(m$family_id),
    family_generation = interaction(m$family_id, m$generation, drop = TRUE))
  for (s in seq_along(schemes)) {
    name <- names(schemes)[s]; block <- as.character(schemes[[s]])
    P <- fl_dv_permutations(block, cfg$dv_permutations, cfg$dv_seed + s)
    saveRDS(P[seq_len(min(10, nrow(P))), , drop = FALSE], file.path(private_dir, paste0("permutations_first10_", name, ".rds")))
    cat("08 beta optimized kernel: ", name, "\n", sep = ""); flush.console()
    Fnull <- fl_dv_fast_null(setup, P, cfg, private_dir, name)
    test <- vegan::permutest(fit, permutations = P[seq_len(min(10, nrow(P))), , drop = FALSE], model = "reduced", parallel = 1L)
    equivalence_error <- max(abs(head(Fnull, length(test$F.perm)) - as.numeric(test$F.perm)))
    if (equivalence_error > 1e-8) stop("Optimized permutation kernel differs from vegan: ", equivalence_error)
    observed <- as.numeric(test$F.0); residual_df <- test$df[2]
    if (any(!is.finite(Fnull))) stop("Nonfinite permutation statistic")
    # vegan's F statistic is nondirectional: its upper tail is the omnibus P, never doubled.
    p <- (1 + sum(Fnull >= observed - sqrt(.Machine$double.eps))) / (length(Fnull) + 1)
    results[[name]] <- data.frame(scheme = name, n = nrow(m), families = nlevels(m$family_id), cases = sum(m$fatty_liver),
      pseudo_F = observed, df_exposure = 1L, df_residual = residual_df,
      R2_total_corrected = fit$CCA$tot.chi / fit$tot.chi,
      R2_partial = fit$CCA$tot.chi / (fit$CCA$tot.chi + fit$CA$tot.chi), p_raw = p,
      permutations = length(Fnull), p_resolution = 1 / (length(Fnull) + 1), p_adjustment = "none",
      test = "nondirectional omnibus F; upper-tail permutation P", distance = "Bray-Curtis + Lingoes correction",
      estimand = "within-family composition association conditional on covariates",
      seed = cfg$dv_seed + s)
    nulls[[name]] <- data.frame(scheme = name, iteration = seq_along(Fnull), F_null = Fnull)
    sizes <- table(block)
    discordant <- vapply(split(m$fatty_liver, block), function(v) length(unique(v)) > 1L, logical(1))
    blocks_audit[[name]] <- data.frame(scheme = name, blocks = length(sizes), singleton_blocks = sum(sizes == 1L),
      movable_people = sum(sizes[sizes > 1L]), discordant_blocks = sum(discordant),
      people_in_discordant_blocks = sum(sizes[names(discordant)[discordant]]),
      condition_columns = ncol(Z), condition_rank = qr(Z)$rank,
      exposure_residual_SS = sum(x_residual^2), lingoes_constant = if (is.null(fit$ac)) 0 else fit$ac,
      optimized_kernel_vs_vegan_first10_max_error = equivalence_error,
      permutation_block_check = "all rows verified")
    # Remove only this function's exact temporary binary/CSV/log paths after checking.
    temporary <- paste0(file.path(private_dir, paste0("work_kernel_", name)), c("_E.bin", "_Q.bin", "_P.bin", "_F.csv", ".log"))
    if (any(normalizePath(dirname(temporary), winslash = "/", mustWork = TRUE) !=
        normalizePath(private_dir, winslash = "/", mustWork = TRUE))) stop("Unsafe temporary cleanup")
    unlink(temporary, recursive = FALSE)
    rm(P); gc(FALSE)
  }
  list(results = do.call(rbind, results), null_distribution = do.call(rbind, nulls),
       exchangeability = do.call(rbind, blocks_audit), fit = fit)
}

fl_dv_descriptive <- function(m) {
  do.call(rbind, lapply(split(seq_len(nrow(m)), interaction(m$generation, m$fatty_liver, drop = TRUE)), function(j) {
    data.frame(generation = as.character(m$generation[j[1]]), fatty_liver = m$fatty_liver[j[1]], n = length(j),
      families = length(unique(m$family_id[j])),
      Shannon_mean = mean(m$shannon[j]), Shannon_SD = sd(m$shannon[j]), Shannon_median = median(m$shannon[j]),
      Shannon_Q1 = unname(quantile(m$shannon[j], .25)), Shannon_Q3 = unname(quantile(m$shannon[j], .75)),
      observed_SGB_median = median(m$observed_sgbs[j]), observed_SGB_Q1 = unname(quantile(m$observed_sgbs[j], .25)),
      observed_SGB_Q3 = unname(quantile(m$observed_sgbs[j], .75)))
  }))
}

fl_dv_dispersion <- function(m, D) {
  all <- list()
  for (g in levels(m$generation)) {
    j <- which(m$generation == g)
    bd <- vegan::betadisper(as.dist(as.matrix(D)[j, j]), factor(m$fatty_liver[j]),
      type = "median", bias.adjust = TRUE, add = "lingoes")
    for (status in c(0L, 1L)) {
      v <- bd$distances[m$fatty_liver[j] == status]
      all[[paste(g, status)]] <- data.frame(generation = g, fatty_liver = status, n = length(v),
        mean_distance = mean(v), median_distance = median(v), Q1 = unname(quantile(v, .25)),
        Q3 = unname(quantile(v, .75)), analysis = "unadjusted descriptive dispersion; no inferential P")
    }
  }
  do.call(rbind, all)
}

fl_dv_save_plot <- function(p, path, width = 183, height = 115) {
  ggplot2::ggsave(paste0(path, ".png"), p, width = width, height = height, units = "mm", dpi = 600, bg = "white")
  ggplot2::ggsave(paste0(path, ".pdf"), p, width = width, height = height, units = "mm", device = grDevices::cairo_pdf)
  if (requireNamespace("svglite", quietly = TRUE))
    ggplot2::ggsave(paste0(path, ".svg"), p, width = width, height = height, units = "mm", device = svglite::svglite)
}

fl_dv_figures <- function(m, primary, main_fit, out, audit) {
  theme <- ggplot2::theme_classic(base_size = 10, base_family = "Arial") +
    ggplot2::theme(legend.position = "top", plot.title = ggplot2::element_text(size = 12),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0), plot.margin = ggplot2::margin(8, 12, 8, 8))
  m$status <- factor(m$fatty_liver, levels = c(0, 1), labels = c("No steatosis", "Steatosis"))
  m$display_group <- factor(paste(m$generation, m$status, sep = "\n"),
    levels = unlist(lapply(c("G0", "G1", "G2"), function(g) paste(g, levels(m$status), sep = "\n"))))
  count <- as.data.frame(table(m$display_group)); names(count) <- c("display_group", "n")
  p <- ggplot2::ggplot(m, ggplot2::aes(display_group, shannon, fill = status)) +
    ggplot2::geom_boxplot(width = .6, outlier.shape = NA, linewidth = .4, alpha = .7) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = .15, height = 0, seed = 82),
      size = .7, alpha = .3, colour = "#24384A") +
    ggplot2::geom_text(data = count, ggplot2::aes(x = display_group, y = Inf, label = paste0("n=", n)),
      inherit.aes = FALSE, vjust = 1.6, size = 3) +
    ggplot2::scale_fill_manual(values = c("#87CEEB", "#DDA66C")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.04, .12))) +
    ggplot2::labs(x = NULL, y = "Shannon diversity (natural log)", fill = NULL,
      title = "Microbiome diversity by generation and steatosis",
      caption = paste0("n=", nrow(m), ". Boxes: median and IQR; whiskers: 1.5 IQR.\nPoints: individual observations; descriptive, unadjusted distributions.")) + theme
  fl_dv_save_plot(p, file.path(out, "Figure2_Shannon_distribution"))
  primary$generation <- factor(primary$generation, levels = c("G2", "G1", "G0"))
  p <- ggplot2::ggplot(primary, ggplot2::aes(estimate, generation)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "#777777", linewidth = .4) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high), height = .12, linewidth = .7) +
    ggplot2::geom_point(size = 2.5, colour = "#CC247C") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .15)) +
    ggplot2::labs(x = "Adjusted Shannon difference (steatosis minus no steatosis)", y = NULL,
      title = "Steatosis–diversity association within each generation",
      caption = if ("review_primary" %in% names(primary))
        "Points: Gaussian GEE; bars: 95% family-sandwich confidence intervals.\nCore review: age, sex, batch, input total, geography, bowel symptoms and sampling factors."
      else "Points: Gaussian GEE estimates; bars: 95% family-sandwich confidence intervals.\nAdjusted for generation-specific age splines, sex, preparation batch and input matrix total.") + theme
  fl_dv_save_plot(p, file.path(out, "Figure3_Shannon_adjusted_difference"), height = 100)
  r <- data.frame(fitted = fitted(main_fit), residual = residuals(main_fit), generation = main_fit$data$generation)
  p <- ggplot2::ggplot(r, ggplot2::aes(fitted, residual, colour = generation)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = .4) +
    ggplot2::geom_point(size = .7, alpha = .4) +
    ggplot2::scale_colour_manual(values = c(G0 = "#CC247C", G1 = "#5292F7", G2 = "#4EA660")) +
    ggplot2::labs(x = "Fitted Shannon", y = "Residual", title = "Primary GEE: residuals versus fitted values") + theme
  fl_dv_save_plot(p, file.path(audit, "GEE_residual_diagnostic"))
}

fl_dv_html_table <- function(x) {
  esc <- function(v) { v <- gsub("&", "&amp;", as.character(v), fixed = TRUE); v <- gsub("<", "&lt;", v, fixed = TRUE); v }
  numeric <- vapply(x, is.numeric, logical(1))
  x[numeric] <- lapply(x[numeric], function(v) ifelse(is.na(v), "", format(signif(v, 5), trim = TRUE)))
  paste0("<div class='scroll'><table><thead><tr>", paste0("<th>", esc(names(x)), "</th>", collapse = ""),
    "</tr></thead><tbody>", paste(apply(x, 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>")), collapse = ""), "</tbody></table></div>")
}

fl_dv_report <- function(alpha, interaction, beta, out, version) {
  text <- c("# 08：脂肪肝与菌群整体多样性", "", paste0("输入版本：", version), "",
    "仅纳入既定三代家庭；主分析包含未成年人。Shannon 为个体内多样性主指标，检出 SGB 数仅描述。",
    "Shannon 以 Gaussian identity GEE 估计脂肪肝组减非脂肪肝组的均值差，家系稳健标准误；同一联合模型给出各代系数及2自由度代际异质性检验。",
    "beta 使用完整SGB相对丰度的Bray–Curtis距离与Lingoes校正；partial dbRDA以家系及预先指定协变量为Condition。其检验目标为家庭内条件组成关联，与alpha边际均值差不同。",
    "主置换允许家系内跨代残差交换；敏感性仅允许同家系同代残差交换。后者G1不能交换，置换范围及信息量须结合exchangeability表阅读。",
    "GEE为原始双侧P；组成检验为无方向F检验，按上尾计算原始P，不将F检验P加倍。均未作多重比较校正。",
    "R2_total_corrected为脂肪肝解释的Lingoes修正总惯量比例；R2_partial为控制协变量后剩余惯量中的比例，两者不可混用。",
    "分代dispersion表为未调整离散程度描述，不提供P，也不能证明残差交换性。年龄与代际高度相关，即使调整也不构成独立代际因果效应。",
    "adult_same_sample_base与adult_BMI使用完全相同的成人BMI完整样本；不把其与全年龄结果的差异归因于BMI，不作中介解释。",
    "原始counts按此前约定解释，上游生成来源未独立核实；matrix_total称输入矩阵总量，不称已核实的测序reads深度。",
    "家庭标签不是核验的亲子/祖孙关系；关联不是菌株传播或因果机制。本步未筛具体菌种、未预测、未中介。", "",
    "具体方法文献及本队列适配见 scripts/08_microbiome_steatosis_diversity_README.md；来源和代码指纹在audit。")
  writeLines(enc2utf8(text), file.path(out, "results_zh.md"), useBytes = TRUE)
  primary <- alpha[alpha$scenario == "primary_all_ages", ]
  page <- paste0("<!doctype html><html lang='zh-CN'><meta charset='utf-8'><title>08 脂肪肝与菌群多样性</title>",
    "<style>body{font:16px/1.7 Arial,'Microsoft YaHei',sans-serif;max-width:1200px;margin:30px auto;padding:15px;color:#24384a}table{border-collapse:collapse;font-size:13px}th,td{padding:7px;border:1px solid #d5dce2;white-space:nowrap}th{background:#edf3f7}.scroll{overflow:auto}img{max-width:100%}</style>",
    "<h1>08 脂肪肝与菌群整体多样性</h1><p>主分析包含319个家庭、1275人、233例超声脂肪肝；保留264名未成年人。</p>",
    "<h2>主模型：各代Shannon均值差</h2>", fl_dv_html_table(primary),
    "<h2>代际异质性检验</h2>", fl_dv_html_table(interaction),
    "<h2>家庭内条件组成关联</h2>", fl_dv_html_table(beta),
    "<h2>全部Shannon模型（包括敏感性和成人同样本扩展）</h2>", fl_dv_html_table(alpha),
    "<h2>结果图</h2><img src='Figure2_Shannon_distribution.png'><img src='Figure3_Shannon_adjusted_difference.png'>",
    "<h2>解释边界</h2><p>", paste(text[-c(1:4)], collapse = "</p><p>"), "</p></html>")
  writeLines(enc2utf8(page), file.path(out, "microbiome_steatosis_diversity.html"), useBytes = TRUE)
}

fl_dv_run <- function(cfg) {
  context <- fl_output_start_run(cfg); cfg <- context$cfg
  audit <- file.path(context$run_dir, "audit"); out <- file.path(context$run_dir, "diversity")
  private <- file.path(context$run_dir, "private")
  dir.create(audit); dir.create(out); dir.create(private)
  writeLines("RUNNING", file.path(context$run_dir, "run_status.txt"))
  tryCatch({
    d <- fl_load_processed_microbiome(cfg$project_root, version = cfg$dv_data_version)
    version <- d$version; m <- fl_dv_prepare(d$metadata, d$full$relative, d$full$counts)
    D <- vegan::vegdist(t(d$full$relative), method = "bray")
    cat("08 input: ", nrow(m), " people; ", length(unique(m$family_id)), " families; ", nrow(d$full$relative), " SGBs\n", sep = "")
    masks <- list(primary_all_ages = rep(TRUE, nrow(m)),
      no_recent_antibiotics = !is.na(m$ANTBV18) & m$ANTBV18 == 0,
      adult_same_sample_base = m$age >= 18 & !is.na(m$BMI18), adult_BMI = m$age >= 18 & !is.na(m$BMI18))
    stopifnot(identical(masks$adult_same_sample_base, masks$adult_BMI))
    fits <- lapply(names(masks), function(s) {
      cat("08 alpha: ", s, "\n", sep = ""); flush.console()
      fl_dv_alpha_fit(m[masks[[s]], , drop = FALSE], cfg, s, bmi = identical(s, "adult_BMI"))
    }); names(fits) <- names(masks)
    alpha <- do.call(rbind, lapply(fits, `[[`, "results"))
    heterogeneity <- do.call(rbind, lapply(fits, `[[`, "interaction"))
    saveRDS(fits, file.path(private, "alpha_models.rds"), version = 3)
    fl_write_csv(alpha, file.path(out, "Table2_Shannon_models.csv"))
    fl_write_csv(heterogeneity, file.path(out, "Table2_generation_heterogeneity.csv"))
    fl_write_csv(fl_dv_descriptive(m), file.path(out, "Table2_diversity_descriptive.csv"))
    fl_write_csv(do.call(rbind, lapply(fits, `[[`, "diagnostic")), file.path(audit, "alpha_model_diagnostics.csv"))
    fl_write_csv(do.call(rbind, lapply(fits, `[[`, "residual_summary")), file.path(audit, "alpha_residual_by_generation_status.csv"))
    fl_write_csv(do.call(rbind, lapply(fits, `[[`, "age_knots")), file.path(audit, "age_spline_knots.csv"))
    flow <- do.call(rbind, lapply(names(masks), function(s) data.frame(scenario = s,
      n = sum(masks[[s]]), cases = sum(m$fatty_liver[masks[[s]]]),
      families = length(unique(m$family_id[masks[[s]]])), excluded_from_baseline = sum(!masks[[s]]))))
    fl_write_csv(flow, file.path(out, "analysis_samples.csv"))
    fl_dv_figures(m, alpha[alpha$scenario == "primary_all_ages", ], fits[[1]]$fit, out, audit)
    beta <- fl_dv_beta_fit(m, D, cfg, private)
    fl_write_csv(beta$results, file.path(out, "Table3_composition_association.csv"))
    fl_write_csv(beta$exchangeability, file.path(out, "Table3_exchangeability.csv"))
    fl_write_csv(beta$null_distribution, file.path(audit, "beta_null_distribution.csv"))
    fl_write_csv(fl_dv_dispersion(m, D), file.path(out, "Table3_dispersion_descriptive.csv"))
    fl_dv_report(alpha, heterogeneity, beta$results, out, version)
    reviewed <- fl_cv_load(cfg, m)
    mr <- fl_dv_prepare(reviewed$metadata, d$full$relative, d$full$counts)
    review <- fl_dv_review_compute(mr, D, cfg, private, out, audit)
    fl_dv_review_report(review, out, reviewed$version)
    fl_write_csv(data.frame(covariate_version = reviewed$version), file.path(audit, "covariate_version.csv"))
    paths <- c(file.path(cfg$processed_microbiome_root, version, "microbiome_processed.rds"), reviewed$source_files, fl_dv_code_paths(cfg))
    fl_write_csv(data.frame(path = paths, md5 = unname(tools::md5sum(paths))), file.path(audit, "input_code_fingerprints.csv"))
    file.copy(file.path(cfg$project_root, "config/microbiome_steatosis_diversity_references.csv"), file.path(audit, "references.csv"))
    fl_write_csv(data.frame(data_version = version, n = nrow(m), families = length(unique(m$family_id)),
      cases = sum(m$fatty_liver), minors = sum(m$age < 18), features = nrow(d$full$relative),
      family_rule = "G0_G1_G2_each_at_least_one", primary_alpha = "Shannon GEE with case-generation interaction",
      beta = "partial dbRDA; family and family-generation reduced residual permutations", permutations = cfg$dv_permutations,
      p_policy = "raw: GEE two-sided; F permutation upper-tail omnibus; no multiplicity adjustment"),
      file.path(audit, "analysis_contract.csv"))
    capture.output(sessionInfo(), file = file.path(audit, "sessionInfo.txt"))
    writeLines(c("COMPLETED", "Exploratory community-level associations; inspect model and exchangeability diagnostics."),
      file.path(context$run_dir, "run_status.txt"))
    final <- fl_output_complete_run(context)
    cat("08 completed: ", file.path(final, "diversity/microbiome_steatosis_diversity.html"), "\n", sep = "")
    invisible(list(alpha = alpha, heterogeneity = heterogeneity, beta = beta$results, review = review, final = final))
  }, error = function(e) { fl_output_fail_run(context, conditionMessage(e)); stop(e) })
}
