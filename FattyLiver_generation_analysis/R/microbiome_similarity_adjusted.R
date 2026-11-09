# Confounder review for 07: Freedman--Lane residual MRQAP, studentized t.
# Dekker et al. 2007, section 3.2.1; Beghini et al. Nature 2025.
# Two disjoint generation roles and sex/age exchangeability restrictions are
# study adaptations. No independent-dyad standard errors or P values are used.

fl_msa_person_basis <- function(m, covariates) {
  covariates <- unique(c("age", "sex", "prep_batch", "log_total", covariates))
  if (!all(covariates %in% names(m))) stop("Missing adjustment fields: ",
    paste(setdiff(covariates, names(m)), collapse = ", "))
  pieces <- list(); constants <- character()
  for (v in covariates) {
    x <- m[[v]]
    if (anyNA(x) || (is.numeric(x) && any(!is.finite(x)))) stop("Incomplete covariate: ", v)
    if (length(unique(x)) < 2L) { constants <- c(constants, v); next }
    if (v == "age") {
      # Exactly three spline df where estimable; otherwise explicitly recorded
      # lower-dimensional fallback for small or synthetic complete-case sets.
      df <- min(3L, length(unique(x)) - 1L)
      z <- if (df == 1L) matrix(as.numeric(x), ncol = 1L) else splines::ns(x, df = df)
      colnames(z) <- paste0("age_ns", seq_len(ncol(z)))
    } else if (is.factor(x) || is.character(x) || is.logical(x)) {
      f <- droplevels(factor(x))
      z <- stats::model.matrix(~ f)[, -1L, drop = FALSE]
      colnames(z) <- paste0(v, "=", levels(f)[-1L])
    } else {
      z <- matrix(as.numeric(x), ncol = 1L, dimnames = list(NULL, v))
    }
    pieces[[v]] <- z
  }
  basis <- if (length(pieces)) do.call(cbind, pieces) else matrix(numeric(), nrow(m), 0L)
  list(basis = basis, constant_covariates = constants, covariates = covariates)
}

fl_msa_design <- function(m, D, generations, covariates) {
  a <- which(as.character(m$generation) == generations[1L])
  b <- which(as.character(m$generation) == generations[2L])
  if (length(a) < 3L || length(b) < 3L) stop("Fewer than three participants in a generation role")
  same <- outer(as.character(m$family_id[a]), as.character(m$family_id[b]), "==")
  if (!any(same) || all(same)) stop("Both within-family and between-family dyads are required")
  old <- fl_msa_person_basis(m[a, , drop = FALSE], covariates)
  young <- fl_msa_person_basis(m[b, , drop = FALSE], covariates)
  Z <- cbind(intercept = 1, old$basis[rep(seq_along(a), times = length(b)), , drop = FALSE],
    young$basis[rep(seq_along(b), each = length(a)), , drop = FALSE])
  colnames(Z) <- c("intercept", paste0("older:", colnames(old$basis)),
    paste0("younger:", colnames(young$basis)))
  if (length(unique(m$prep_batch[c(a, b)])) > 1L)
    Z <- cbind(Z, same_prep_batch = as.numeric(outer(as.character(m$prep_batch[a]),
      as.character(m$prep_batch[b]), "==")))
  if ("survey_center" %in% covariates && length(unique(m$survey_center[c(a, b)])) > 1L)
    Z <- cbind(Z, same_survey_center = as.numeric(outer(as.character(m$survey_center[a]),
      as.character(m$survey_center[b]), "==")))
  if ("seq_pool" %in% covariates && length(unique(m$seq_pool[c(a, b)])) > 1L)
    Z <- cbind(Z, same_seq_pool = as.numeric(outer(as.character(m$seq_pool[a]), as.character(m$seq_pool[b]), "==")))
  # Scale before QR to prevent unrelated units (energy, age, DNA) changing rank.
  if (ncol(Z) > 1L) for (j in 2:ncol(Z)) {
    s <- stats::sd(Z[, j]); if (is.finite(s) && s > 0) Z[, j] <- (Z[, j] - mean(Z[, j])) / s
  }
  qrZ <- qr(Z, tol = 1e-9, LAPACK = FALSE)
  rank <- qrZ$rank
  keep <- qrZ$pivot[seq_len(rank)]
  aliases <- if (rank < ncol(Z)) colnames(Z)[qrZ$pivot[(rank + 1L):ncol(Z)]] else character()
  Z <- Z[, keep, drop = FALSE]
  Q <- qr.Q(qr(Z, tol = 1e-9), complete = FALSE)
  x <- as.numeric(same); y <- as.numeric(D[a, b, drop = FALSE])
  xres <- as.numeric(x - Q %*% crossprod(Q, x))
  e <- as.numeric(y - Q %*% crossprod(Q, y))
  ss <- sum(xres * xres); df <- length(y) - rank - 1L
  if (!is.finite(ss) || ss < 1e-9 || df < 2L) stop("Insufficient residual family-exposure variation or degrees of freedom")
  xunit <- xres / sqrt(ss); num <- sum(xunit * e)
  rss <- sum(e * e) - num * num
  if (!is.finite(rss) || rss <= 0) stop("Non-positive observed residual variance")
  list(a = a, b = b, same = same, y = y, x = x, Z = Z, Q = Q,
    residual = e, fitted = y - e, xunit = xunit, exposure_residual_ss = ss,
    nuisance_rank = rank, nuisance_aliases = aliases,
    constants = c(paste0("older:", old$constant_covariates), paste0("younger:", young$constant_covariates)),
    df = df, estimate = num / sqrt(ss), t = num / sqrt(rss / df))
}

fl_msa_blocks <- function(m, width) {
  key <- interaction(as.character(m$sex), floor(m$age / width), drop = TRUE, lex.order = TRUE)
  groups <- unname(split(seq_len(nrow(m)), key))
  movable <- cross <- rep(FALSE, nrow(m))
  for (z in groups) {
    movable[z] <- length(z) > 1L
    cross[z] <- length(unique(m$family_id[z])) > 1L
  }
  list(groups = groups, movable = movable, cross = cross,
    blocks = length(groups), singletons = sum(lengths(groups) == 1L))
}

fl_msa_node_permutations <- function(blocks, n, B) {
  out <- matrix(rep(seq_len(n), each = B), nrow = B)
  for (i in seq_len(B)) for (g in blocks$groups) if (length(g) > 1L)
    out[i, g] <- g[sample.int(length(g))]
  out
}

fl_msa_lm_check <- function(z, rows, cols, values) {
  E <- matrix(z$residual, length(z$a), length(z$b))
  errors <- numeric(min(5L, nrow(rows)))
  for (i in seq_along(errors)) {
    Ystar <- z$fitted + as.numeric(E[rows[i, ], cols[i, ], drop = FALSE])
    # Deliberately use an ordinary lm fit solely to check the statistic; its
    # independent-row P/CI are never used as inferential output.
    fit <- stats::lm(Ystar ~ z$Z + z$x - 1)
    ordinary_t <- unname(summary(fit)$coefficients[ncol(z$Z) + 1L, "t value"])
    errors[i] <- abs(ordinary_t - values[i])
  }
  if (any(!is.finite(errors)) || max(errors) > 1e-7)
    stop("MRQAP studentized t differs from independent lm fit: ", max(errors))
  max(errors)
}

fl_msa_kernel <- function(z, rows, cols, cfg, prefix) {
  files <- paste0(prefix, c("_E.bin", "_Q.bin", "_X.bin", "_ROW.bin", "_COL.bin"))
  binary <- list(as.numeric(z$residual), as.numeric(z$Q), as.numeric(z$xunit),
    as.integer(rows), as.integer(cols))
  for (i in seq_along(files)) {
    con <- file(files[i], "wb"); tryCatch(writeBin(binary[[i]], con, endian = "little"), finally = close(con))
  }
  args <- c(shQuote(file.path(cfg$project_root, "scripts/mrqap_permutation_kernel.py")),
    "--prefix", shQuote(prefix), "--na", length(z$a), "--nb", length(z$b),
    "--k", ncol(z$Q), "--permutations", nrow(rows), "--df", z$df)
  status <- system2(cfg$python, args = args, stdout = paste0(prefix, ".log"), stderr = paste0(prefix, ".log"))
  if (status != 0L) stop("MRQAP Python kernel failed; see ", prefix, ".log")
  values <- utils::read.csv(paste0(prefix, "_t.csv"))$t_null
  if (length(values) != nrow(rows) || any(!is.finite(values))) stop("Incomplete MRQAP kernel output")
  error <- fl_msa_lm_check(z, rows, cols, values)
  # Keep first five node-index draws for independent replication, without IDs.
  utils::write.csv(data.frame(permutation = seq_len(min(5L, nrow(rows))),
    rows[seq_len(min(5L, nrow(rows))), , drop = FALSE]), paste0(prefix, "_first5_rows.csv"), row.names = FALSE)
  utils::write.csv(data.frame(permutation = seq_len(min(5L, nrow(cols))),
    cols[seq_len(min(5L, nrow(cols))), , drop = FALSE]), paste0(prefix, "_first5_cols.csv"), row.names = FALSE)
  utils::write.csv(data.frame(permutation = seq_len(min(5L, length(values))),
    t_null = head(values, 5L)), paste0(prefix, "_first5_t.csv"), row.names = FALSE)
  # Exact files created above only; never recursively remove a computed path.
  unlink(c(files, paste0(prefix, "_t.csv")))
  list(values = values, error = error)
}

fl_ms_adjusted_compute <- function(m, relative, cfg, private_dir) {
  if (!requireNamespace("vegan", quietly = TRUE)) stop("vegan is required")
  if (anyDuplicated(m$sample_id) || !identical(colnames(relative), as.character(m$sample_id)) ||
      any(!is.finite(relative)) || any(relative < 0) || max(abs(colSums(relative) - 1)) > 1e-10)
    stop("Invalid or misaligned complete-denominator relative abundance input")
  B <- as.integer(cfg$ms_permutations)
  if (length(B) != 1L || is.na(B) || B < 19L || B != cfg$ms_permutations)
    stop("Invalid number of permutations")
  if (!is.finite(cfg$ms_age_band_years) || cfg$ms_age_band_years <= 0) stop("Invalid age band width")
  dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)
  D <- as.matrix(vegan::vegdist(t(relative), method = "bray"))
  scenarios <- fl_cv_scenarios(m, "similarity")
  selected <- cfg$ms_review_scenarios
  if (is.null(selected) || !length(selected) || !all(selected %in% names(scenarios)))
    stop("Select known ms_review_scenarios before running adjusted 07")
  results <- summaries <- exchanges <- nulls <- list()
  for (s in seq_along(selected)) {
    scenario <- selected[s]; sc <- scenarios[[scenario]]
    if (length(sc$mask) != nrow(m)) stop("Scenario mask length mismatch: ", scenario)
    required <- unique(c("age", "sex", "prep_batch", "log_total", sc$covariates))
    if (!all(required %in% names(m))) stop("Missing scenario variables: ", paste(setdiff(required, names(m)), collapse = ","))
    valid <- !is.na(sc$mask) & sc$mask & stats::complete.cases(m[, required, drop = FALSE])
    for (v in required) if (is.numeric(m[[v]])) valid <- valid & is.finite(m[[v]])
    idx <- which(valid); mm <- m[idx, , drop = FALSE]; dd <- D[idx, idx, drop = FALSE]
    for (p in seq_along(cfg$ms_pairs)) {
      g <- cfg$ms_pairs[[p]]; comparison <- paste(g, collapse = "-")
      id <- paste(scenario, comparison, sep = ":")
      cat("07 adjusted MRQAP: ", id, "; selected people=", nrow(mm), "\n", sep = ""); flush.console()
      # Catch model-identifiability issues only. Kernel or numerical validation
      # failures propagate, so computation failures cannot masquerade as results.
      z <- tryCatch(fl_msa_design(mm, dd, g, sc$covariates), error = function(e) e)
      common <- data.frame(scenario = scenario, scenario_label = sc$label, comparison = comparison,
        n_people = sum(as.character(mm$generation) %in% g),
        n_people_before = sum(as.character(m$generation) %in% g),
        n_people_excluded = sum(as.character(m$generation) %in% g) - sum(as.character(mm$generation) %in% g),
        families = length(unique(mm$family_id[as.character(mm$generation) %in% g])),
        covariates = paste(required, collapse = ";"),
        adjustment_design = "both endpoint covariates; age ns3; same prep batch; same sequencing pool; same center when selected",
        p_adjustment = "none; raw absolute studentized t two-sided MRQAP",
        interval_type = "No independent-dyad confidence interval; null reference is not effect CI")
      if (inherits(z, "error")) {
        common$status <- "NOT_ESTIMABLE"; common$reason <- conditionMessage(z)
        results[[id]] <- common; next
      }
      ob <- fl_msa_blocks(mm[z$a, , drop = FALSE], cfg$ms_age_band_years)
      yb <- fl_msa_blocks(mm[z$b, , drop = FALSE], cfg$ms_age_band_years)
      for (kind in c("IF", "BF")) {
        values <- z$y[if (kind == "IF") as.numeric(z$same) == 1 else as.numeric(z$same) == 0]
        summaries[[paste(id, kind)]] <- data.frame(scenario = scenario, comparison = comparison,
          pair_type = kind, pairs = length(values), mean_distance = mean(values),
          median_distance = stats::median(values), Q1 = unname(stats::quantile(values, .25)),
          Q3 = unname(stats::quantile(values, .75)))
      }
      common$within_pairs <- sum(z$same); common$between_pairs <- sum(!z$same)
      common$observed_within_mean <- mean(z$y[as.numeric(z$same) == 1])
      common$observed_between_mean <- mean(z$y[as.numeric(z$same) == 0])
      common$observed_within_minus_between <- common$observed_within_mean - common$observed_between_mean
      common$adjusted_distance_difference <- z$estimate; common$observed_t <- z$t
      common$nuisance_rank <- z$nuisance_rank; common$residual_df <- z$df
      common$nuisance_aliases <- paste(z$nuisance_aliases, collapse = ";")
      common$constant_covariates <- paste(z$constants, collapse = ";")
      common$exposure_residual_ss <- z$exposure_residual_ss
      exchanges[[id]] <- data.frame(scenario = scenario, comparison = comparison,
        older_n = length(z$a), younger_n = length(z$b), older_blocks = ob$blocks, younger_blocks = yb$blocks,
        older_singletons = ob$singletons, younger_singletons = yb$singletons,
        older_movable = sum(ob$movable), younger_movable = sum(yb$movable),
        older_cross_family_movable = sum(ob$cross), younger_cross_family_movable = sum(yb$cross),
        within_edges_both_endpoints_not_cross_family_movable = sum(z$same & outer(!ob$cross, !yb$cross, "&")),
        exchangeability = paste0("separate generation roles; sex x ", cfg$ms_age_band_years, "-year age band"))
      if (!any(ob$cross) && !any(yb$cross)) {
        common$status <- "NOT_ESTIMABLE"; common$reason <- "No nodes can exchange between families within blocks"
        results[[id]] <- common; next
      }
      set.seed(cfg$ms_seed + 1000L * s + p)
      rows <- fl_msa_node_permutations(ob, length(z$a), B)
      cols <- fl_msa_node_permutations(yb, length(z$b), B)
      prefix <- file.path(private_dir, paste0("mrqap_", gsub("[^A-Za-z0-9_]", "_", id)))
      kernel <- fl_msa_kernel(z, rows, cols, cfg, prefix)
      tt <- kernel$values
      effective <- stats::sd(tt) > 1e-10
      common$p_raw_two_sided <- if (effective) (1 + sum(abs(tt) >= abs(z$t) - 1e-12)) / (B + 1) else NA_real_
      common$permutations <- B; common$null_t_mean <- mean(tt); common$null_t_sd <- stats::sd(tt)
      common$null_t_reference_low <- unname(stats::quantile(tt, .025))
      common$null_t_reference_high <- unname(stats::quantile(tt, .975))
      common$first5_lm_statistic_max_error <- kernel$error
      common$status <- if (effective) "ESTIMATED" else "NOT_ESTIMABLE"
      common$reason <- if (effective) "" else "No effective permutation statistic variation"
      results[[id]] <- common
      nulls[[id]] <- data.frame(scenario = scenario, comparison = comparison,
        permutation = seq_len(B), t_null = tt)
      rm(z, rows, cols); invisible(gc(FALSE))
    }
  }
  bind <- function(items) {
    if (!length(items)) return(data.frame())
    columns <- unique(unlist(lapply(items, names)))
    items <- lapply(items, function(x) { for (nm in setdiff(columns, names(x))) x[[nm]] <- NA; x[, columns, drop = FALSE] })
    value <- do.call(rbind, items); rownames(value) <- NULL; value
  }
  list(adjusted_results = bind(results), adjusted_pair_summary = bind(summaries),
    adjusted_exchangeability = bind(exchanges), adjusted_null_distribution = bind(nulls))
}
