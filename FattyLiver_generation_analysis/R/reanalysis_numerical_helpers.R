# Frozen audited numerical helpers. Active designs are defined in reanalysis_20261024.R.

ov_bind <- function(items) {
  items<-Filter(function(z)!is.null(z)&&nrow(z)>0,items)
  if(!length(items))return(data.frame())
  fields<-unique(unlist(lapply(items,names)))
  items<-lapply(items,function(z){for(n in setdiff(fields,names(z)))z[[n]]<-NA;z[,fields,drop=FALSE]})
  z<-do.call(rbind,items);rownames(z)<-NULL;z
}

ov_design <- function(x,mask) {
  ii<-which(mask);m<-droplevels(x$metadata[ii,,drop=FALSE]);nuis<-x$nuisance
  if(!length(ii)||length(unique(m$fatty_liver))!=2||anyNA(m[,nuis]))stop("Invalid sample/design")
  constants<-nuis[vapply(m[,nuis,drop=FALSE],function(z)length(unique(z))<2,logical(1))]
  Z<-stats::model.matrix(stats::reformulate(setdiff(nuis,constants)),m)
  q<-qr(Z);keep<-sort(q$pivot[seq_len(q$rank)]);aliases<-setdiff(colnames(Z),colnames(Z)[keep])
  Z<-Z[,keep,drop=FALSE];Z<-Z[,colnames(Z)!="(Intercept)",drop=FALSE]
  # Scale nuisance design columns only for conditioning; NEVER scale AST outcomes.
  centers<-colMeans(Z);scales<-apply(Z,2,stats::sd);stopifnot(all(is.finite(scales)&scales>0))
  map<-data.frame(model_column=sprintf("nuis_%03d",seq_len(ncol(Z))),original_column=colnames(Z),center=centers,scale=scales)
  Z<-sweep(sweep(Z,2,centers,"-"),2,scales,"/");colnames(Z)<-map$model_column
  md<-data.frame(family_id=m$family_id,fatty_liver=factor(m$fatty_liver,levels=c(0,1),labels=c("No","Yes")),Z,row.names=m$sample_id,check.names=FALSE)
  bases<-"fatty_liver"
  for(v in bases)contrasts(md[[v]])<-stats::contr.treatment(c("No","Yes"))
  M<-stats::model.matrix(stats::reformulate(c(bases,map$model_column)),md)
  if(qr(M)$rank!=ncol(M))stop("Disease contrast is rank deficient")
  a<-x$relative[ii,,drop=FALSE];pv<-colMeans(a>0);included<-pv>=.10
  filter<-data.frame(feature_id=colnames(a),prevalence=pv,retained=included,
    reason=ifelse(included,"RETAINED_PREVALENCE_GE_10_PERCENT","FILTERED_PREVALENCE_LT_10_PERCENT"))
  list(index=ii,metadata=m,design=md,bases=bases,terms=paste0(bases,"Yes"),nuisance=map$model_column,mapping=map,
    formula=paste("~",paste(c(bases,map$model_column,"(1 | family_id)"),collapse=" + ")),
    dropped=c(constants,aliases),filter=filter,n=nrow(m),cases=sum(m$fatty_liver),families=nlevels(m$family_id),generation=FALSE)
}

ov_extract <- function(ans,features,d,scenario,expected) {
  effects<-diagnostics<-tests<-list()
  for(id in features) {
    fit<-ans$fits[[id]];messages<-character();error<-"";cs<-NULL;ft<-NULL
    optcode<-0L;hmin<-NA_real_;singular<-NA;skew<-icc<-NA_real_
    if(inherits(fit,"merMod")) {
      mf<-stats::model.frame(fit);at<-match(rownames(mf),rownames(expected))
      if(anyNA(at)||max(abs(as.numeric(stats::model.response(mf))-expected[at,id]))>1e-10)stop("Official AST response/order does not match asin(sqrt(relative)): ",id)
      messages<-as.character(fit@optinfo$conv$lme4$messages);optcode<-unlist(fit@optinfo$conv$opt)
      H<-fit@optinfo$derivs$Hessian
      if(is.matrix(H)&&all(is.finite(H)))hmin<-min(eigen(H,symmetric=TRUE,only.values=TRUE)$values)
      singular<-lme4::isSingular(fit,tol=1e-4)
      # Reconstitute the three local names in the official MaAsLin2 stored call
      # so public lmerTest conversion can reconstruct its deviance function for
      # derivative diagnostics. This does not optimize/refit the model or alter
      # its coefficients, response or participants.
      data<-mf;formula<-stats::formula(fit);na.action<-stats::na.fail
      ft<-tryCatch(withCallingHandlers(lmerTest::as_lmerModLmerTest(fit),warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e){error<<-conditionMessage(e);NULL})
      if(!is.null(ft))cs<-tryCatch(coef(summary(ft)),error=function(e){error<<-conditionMessage(e);NULL})
      rv<-stats::residuals(fit);sr<-stats::sd(rv)
      if(is.finite(sr)&&sr>0)skew<-mean((rv-mean(rv))^3)/sr^3
      vv<-as.data.frame(lme4::VarCorr(fit));vf<-sum(vv$vcov[vv$grp!="Residual"]);icc<-vf/(vf+stats::sigma(fit)^2)
    }else error<-"Official MaAsLin2 did not return a merMod fit"
    messages<-unique(messages[nzchar(messages)]);conv_messages<-messages[!grepl("boundary.*singular|singular fit",messages,ignore.case=TRUE)]
    bad<-nzchar(error)||length(conv_messages)>0L||any(optcode!=0)||(!is.na(hmin)&&hmin< -1e-6)
    diagnostic_status<-if(bad)"UNRELIABLE_CONVERGENCE"else if(isTRUE(singular))"ESTIMATED_SINGULAR"else"ESTIMATED"
    for(k in seq_along(d$terms)) {
      term<-d$terms[k];b<-se<-df<-p<-lo<-hi<-NA_real_;finite<-FALSE;package_check<-FALSE
      if(!is.null(cs)&&term%in%rownames(cs)&&all(c("Estimate","Std. Error","df","Pr(>|t|)")%in%colnames(cs))) {
        b<-cs[term,"Estimate"];se<-cs[term,"Std. Error"];df<-cs[term,"df"];p<-cs[term,"Pr(>|t|)"]
        finite<-all(is.finite(c(b,se,df,p)))&&se>0&&df>0&&p>=0&&p<=1
        if(finite){critical<-stats::qt(.975,df);lo<-b-critical*se;hi<-b+critical*se}
        rr<-ans$results;rr<-rr[rr$feature==id&rr$metadata==d$bases[k]&as.character(rr$value)=="Yes",,drop=FALSE]
        if(nrow(rr)==1L&&finite)package_check<-isTRUE(all.equal(unname(b),rr$coef[1],tolerance=1e-7))&&isTRUE(all.equal(unname(p),rr$pval[1],tolerance=1e-7))
      }
      state<-if(!finite)"NOT_ESTIMABLE"else if(!package_check)"OFFICIAL_RESULT_MISMATCH"else diagnostic_status
      valid<-state%in%c("ESTIMATED","ESTIMATED_SINGULAR")
      effects[[length(effects)+1L]]<-data.frame(feature_id=id,scenario=scenario,term=term,
        generation="Overall",estimate=b,SE=se,df=df,
        CI_low=lo,CI_high=hi,p_raw=p,status=state,analysis_valid=valid,singular=singular,
        unit="Difference in asin(sqrt(relative abundance)); radians; no outcome SD scaling")
    }
    diagnostics[[length(diagnostics)+1L]]<-data.frame(feature_id=id,scenario=scenario,status=diagnostic_status,
      optimizer_code=paste(optcode,collapse=";"),minimum_optimizer_Hessian_eigenvalue=hmin,singular=singular,
      residual_skewness=skew,family_ICC=icc,warnings=paste(messages,collapse="; "),error=error)

  }
  list(effects=ov_bind(effects),diagnostics=ov_bind(diagnostics),heterogeneity=ov_bind(tests))
}

ov_expand <- function(z,d,features,scenario) {
  parts<-lapply(seq_along(d$terms),function(k) {
    v<-z[z$term==d$terms[k],,drop=FALSE];at<-match(features,v$feature_id);r<-v[at,,drop=FALSE]
    r$feature_id<-features;r$scenario<-scenario;r$term<-d$terms[k]
    r$generation<-"Overall"
    r$status[is.na(at)]<-"FILTERED_PREVALENCE_LT_10_PERCENT";r$analysis_valid[is.na(at)]<-FALSE;r
  });ov_bind(parts)
}

ov_species_label <- function(taxonomy) {
  vapply(taxonomy,function(x) {
    if(is.na(x)||!nzchar(x))return("")
    paths<-strsplit(x,",",fixed=TRUE)[[1]];s<-paths[grepl("s__",paths,fixed=TRUE)]
    s<-unique(sub(".*s__([^|,]+).*","\\1",s));s<-gsub("_"," ",s,fixed=TRUE)
    if(length(s)==1L)s else if(length(s)>1L)"Multiple taxonomic labels; see taxonomy field"else""
  },character(1))
}

ov_export <- function(plot,stem,out,cfg,width=183,height=235) {
  plot<-fl_model_plot_labels(plot)
  # Native facet geometry: unlike a single-panel patchwork wrapper, this audits
  # every comparable facet rectangle at its final physical dimensions.
  audit<-file.path(out,"audit");probe<-file.path(audit,paste0(stem,".layout-probe.pdf"))
  grDevices::cairo_pdf(probe,width=width/25.4,height=height/25.4,family="Arial")
  opened<-TRUE;on.exit(if(opened)grDevices::dev.off(),add=TRUE)
  grob<-ggplot2::ggplotGrob(plot);grob<-fl_model_label_grob(grob);grid::grid.newpage();grid::grid.draw(grob);grid::grid.force()
  rows<-grob$layout[grepl("^panel-[0-9]+-[0-9]+$",grob$layout$name),,drop=FALSE]
  if(!nrow(rows))rows<-grob$layout[grob$layout$name=="panel",,drop=FALSE]
  if(!nrow(rows))stop("No rendered plot area found")
  vp<-grid::grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)$name
  measure<-function(name) {
    hit<-vp[startsWith(vp,paste0(name,"."))];if(length(hit)!=1L)stop("Cannot locate final facet viewport: ",name)
    grid::seekViewport(hit);on.exit(grid::upViewport(0))
    lo<-grid::deviceLoc(grid::unit(0,"npc"),grid::unit(0,"npc"),valueOnly=TRUE)
    hi<-grid::deviceLoc(grid::unit(1,"npc"),grid::unit(1,"npc"),valueOnly=TRUE)
    72*c(lo$x,lo$y,hi$x,hi$y)
  }
  panels<-lapply(seq_len(nrow(rows)),function(i){r<-rows[i,];list(id=r$name,bbox_pt=measure(r$name),grid_id="native-facets",
    row_start=r$t-1L,row_stop=r$b,col_start=r$l-1L,col_stop=r$r)})
  grDevices::dev.off();opened<-FALSE
  layout<-file.path(audit,paste0(stem,".alignment-layout.json"))
  jsonlite::write_json(list(schema_version=1L,backend="r-ggplot2-cairo",measurement="Final native facet viewports",
    figure=list(width_pt=width/25.4*72,height_pt=height/25.4*72),panels=panels,exemptions=list()),layout,auto_unbox=TRUE,pretty=TRUE,digits=NA)
  status<-0L
  if(nrow(rows)>1L)status<-system2(cfg$python,c(shQuote(file.path(nf_skill_root(),"scripts/audit_panel_alignment.py")),shQuote(layout),
    "--json-out",shQuote(file.path(audit,paste0(stem,".alignment.json"))),"--strict"))else
    jsonlite::write_json(list(verdict="NOT_APPLICABLE_SINGLE_PANEL",panels=1L),file.path(audit,paste0(stem,".alignment.json")),auto_unbox=TRUE,pretty=TRUE)
  if(status!=0L)stop("Facet alignment gate blocked export: ",stem)
  pdf<-file.path(out,"figures",paste0(stem,".pdf"))
  draw<-function(){grid::grid.newpage();grid::grid.draw(grob)}
  grDevices::cairo_pdf(pdf,width=width/25.4,height=height/25.4,family="Arial",bg="white");tryCatch(draw(),finally=grDevices::dev.off())
  adapter<-file.path(cfg$project_root,"scripts/audit_nature_cairo_pdf.py")
  status<-system2(cfg$python,c(shQuote(adapter),shQuote(pdf),"--skill",shQuote(nf_skill_root()),"--prefix",shQuote(file.path(audit,stem))))
  if(status!=0L)stop("Rendered text/collision audit blocked export: ",stem)
  svglite::svglite(file.path(out,"figures",paste0(stem,".svg")),width=width/25.4,height=height/25.4,bg="white");tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::png(file.path(out,"figures",paste0(stem,".png")),width=width,height=height,units="mm",res=600,type="cairo",bg="white");tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::tiff(file.path(out,"figures",paste0(stem,".tiff")),width=width,height=height,units="mm",res=600,type="cairo",compression="lzw",bg="white");tryCatch(draw(),finally=grDevices::dev.off())
  data.frame(figure=stem,width_mm=width,height_mm=height,facets=nrow(rows),formats="PDF;SVG;PNG;TIFF",dpi=600)
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
