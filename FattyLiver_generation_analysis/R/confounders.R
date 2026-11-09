# Audited, immutable covariate package shared by standalone entries 06-08.
fl_cv_coldata <- function(path) {
  x <- suppressWarnings(readRDS(path))
  cd <- attributes(x)[["colData"]]; ld <- attributes(cd)[["listData"]]
  if (is.null(cd) && identical(attributes(x)[["class"]], "data.frame")) return(x)
  if (is.null(ld)) stop("Source does not contain sample metadata: ", path)
  rn <- attributes(cd)[["rownames"]]
  # Older S4 serialization stores a NULL symbol instead of a row-name vector;
  # its explicit SampleID column is the validated key used below.
  if (!is.character(rn) && !is.integer(rn)) rn <- seq_len(length(ld[[1L]]))
  structure(ld, class = "data.frame", row.names = rn)
}

fl_cv_rules <- function() {
  data.frame(variable = c("survey_center", "urbanicity", "abx_recent", "constipation", "loose_stools", "probiotic",
    "diarrhea", "gastro_recent", "season_sin", "season_cos", "freeze_delay_log", "stool_gap_log", "log_DNA",
    "smoking", "alcohol", "diet_score", "physical_activity", "energy_log", "fiber", "sugar",
    "diabetes", "diabetes_med", "lipid_med", "bp_med", "rheum_med", "BMI18"),
    raw_field = c("alue18", "kaumaa18", "ANTBV18", "UMMETUS18", "LOYSAVATSA18", "MAITOHAPPO18",
      "RIPULI18", "VATSATAUTI18", "ULOSTE_PVM18", "ULOSTE_PVM18", "UPAKASTUS_PVM18-ULOSTE_PVM18",
      "stool_ultrasound_gap_days", "dna_pitoisuus", "smoke18", "alkodpd18", "dietscore_g18", "PAI18",
      "ENERCK18", "FIBC18", "SUGAR18", "diab218", "DIALAA18", "FHLAAK18", "VPLAAK18", "REULAAK18", "BMI18"),
    role = c(rep("population_or_sampling_core", 6), rep("perturbation_sensitivity", 2),
      rep("sampling_core", 4), "DNA_conditioned_sensitivity", rep("adult_lifestyle", 5), rep("adult_nutrient_alternative", 2),
      rep("adult_metabolic_medication_condition", 5), "adult_adiposity_condition"),
    definition = c("2018 survey municipality, not verified residential city", "Postal-code urban/rural category; original codes retained as factor",
      "Antibiotics in three months before stool", "Constipation in week before stool", "Loose stools in week before stool",
      "Regular lactic-acid bacteria preparation in prior month; not total fermented-food intake",
      "Diarrhea in week before stool", "Gastroenteritis in three months before stool", "Annual cycle sine of stool month",
      "Annual cycle cosine of stool month", "log1p(days from stool collection to recorded aliquoting/freezing); not verified room-temperature transit",
      "log1p(absolute stool-ultrasound interval days)", "log10 measured extracted DNA concentration; may also reflect biological biomass",
      "Current tobacco-product use binary in harmonized source", "Alcohol servings/day; not grams/day", "Diet score 0-27",
      "Adult PAI18 index; children's different scale not pooled", "log10 dietary kcal excluding fiber energy",
      "Dietary fiber grams", "Total mono/disaccharides grams", "2018 type 2 diabetes status, not a core confounder by default",
      "Any diabetes medication; not specific metformin", "Cholesterol medication; not verified statin class",
      "Blood-pressure medication; not a PPI", "Rheumatism medication; not a specific drug", "Adult raw BMI; not applied to minors"),
    stringsAsFactors = FALSE)
}

fl_cv_core <- function() c("survey_center", "urbanicity", "abx_recent", "constipation", "loose_stools", "probiotic",
  "season_sin", "season_cos", "freeze_delay_log", "stool_gap_log")
fl_cv_life <- function() c("smoking", "alcohol", "diet_score", "physical_activity", "energy_log")
fl_cv_drugs <- function() c("diabetes", "diabetes_med", "lipid_med", "bp_med", "rheum_med")

fl_cv_scenarios <- function(m) {
  core <- fl_cv_core(); life <- fl_cv_life(); drugs <- fl_cv_drugs()
  baseline <- c("age", "sex", "generation", "prep_batch", "log_total")
  cc <- function(fields) stats::complete.cases(m[, unique(c(baseline, fields)), drop = FALSE])
  core_mask <- cc(core)
  adult_mask <- m$age >= 18 & cc(c(core, life, "BMI18", drugs))
  nutrient_terms <- c("smoking", "alcohol", "physical_activity", "energy_log", "fiber", "sugar")
  nutrient_mask <- m$age >= 18 & cc(c(core, nutrient_terms))
  unperturbed <- core_mask & !is.na(m$gastro_recent) & m$gastro_recent == 0 & m$abx_recent == 0
  make <- function(mask, covariates, label) list(mask = mask, covariates = covariates, label = label)
  list(core_base = make(core_mask, character(), "Core-complete same-sample technical reference"),
    core_adjusted = make(core_mask, core, "All-age population, bowel, antibiotic, probiotic and sampling adjusted"),
    adult_core = make(adult_mask, core, "Adult complete-case core reference"),
    adult_lifestyle = make(adult_mask, c(core, life), "Adult lifestyle incl physical activity and energy"),
    adult_BMI = make(adult_mask, c(core, life, "BMI18"), "Adult adiposity-conditioned association"),
    adult_medications = make(adult_mask, c(core, life, "BMI18", drugs), "Adult metabolic/medication-conditioned sensitivity"),
    unperturbed_core = make(unperturbed, core, "No recent antibiotics or gastroenteritis; same core adjustment"),
    nutrient_base = make(nutrient_mask, core, "Adult nutrient-complete core reference"),
    nutrient_adjusted = make(nutrient_mask, c(core, nutrient_terms), "Adult energy/fiber/sugar alternative to composite diet score"),
    core_DNA = make(core_mask & !is.na(m$log_DNA), c(core, "log_DNA"), "DNA-conditioned technical sensitivity, not primary"))
}

fl_cv_build <- function(cfg) {
  context <- fl_output_start_run(cfg); cfg <- context$cfg
  pub <- file.path(context$run_dir, "covariates"); audit <- file.path(context$run_dir, "audit"); private <- file.path(context$run_dir, "private")
  dir.create(pub); dir.create(audit); dir.create(private)
  writeLines("RUNNING", file.path(context$run_dir, "run_status.txt"))
  tryCatch({
    d <- fl_load_processed_microbiome(cfg$project_root); m <- d$metadata; microbiome_version <- d$version; rm(d); gc(FALSE)
    primary_path <- cfg$microbiome_file
    nutrient_path <- "E:/YFS/workspace/microbiome/tse_insulin.rds"
    src <- fl_cv_coldata(primary_path)
    if (anyDuplicated(rownames(src)) || anyNA(match(m$sample_id, rownames(src)))) stop("Primary metadata key mismatch")
    src <- src[match(m$sample_id, rownames(src)), , drop = FALSE]
    num <- function(name) {
      i <- match(tolower(name), tolower(names(src)))
      if (is.na(i)) stop("Audited field absent: ", name)
      x <- src[[i]]; if (!is.numeric(x)) stop("Non-numeric field: ", name); as.numeric(x)
    }
    binary <- function(name) {
      v <- num(name); if (any(!is.na(v) & !v %in% c(0, 1))) stop("Unverified binary coding: ", name); v
    }
    dates <- function(name) {
      v <- src[[match(name, names(src))]]
      if (!inherits(v, "Date")) stop("Unverified source date class: ", name)
      v
    }
    m$survey_center <- factor(num("alue18")); m$urbanicity <- factor(num("kaumaa18"))
    for (pair in list(c("abx_recent", "ANTBV18"), c("constipation", "UMMETUS18"),
      c("loose_stools", "LOYSAVATSA18"), c("probiotic", "MAITOHAPPO18"),
      c("diarrhea", "RIPULI18"), c("gastro_recent", "VATSATAUTI18"),
      c("smoking", "smoke18"), c("diabetes", "diab218"), c("diabetes_med", "DIALAA18"),
      c("lipid_med", "FHLAAK18"), c("bp_med", "VPLAAK18"), c("rheum_med", "REULAAK18"))) m[[pair[1]]] <- binary(pair[2])
    # Existing outcome, age, sex and family membership are never replaced from older sources.
    m$alcohol <- m$alkodpd18; m$diet_score <- m$dietscore_g18
    m$physical_activity <- num("PAI18"); m$physical_activity[m$age < 18] <- NA_real_
    if (any(!is.na(m$physical_activity) & !m$physical_activity %in% 5:15)) stop("Adult activity scale invalid")
    stool <- dates("ULOSTE_PVM18"); frozen <- dates("UPAKASTUS_PVM18")
    delay <- as.numeric(frozen - stool)
    if (any(delay < 0, na.rm = TRUE)) stop("Negative stool-freezing interval; source review required")
    month <- as.integer(format(stool, "%m"))
    m$season_sin <- sin(2 * pi * (month - 1) / 12); m$season_cos <- cos(2 * pi * (month - 1) / 12)
    m$freeze_delay_log <- log1p(delay)
    m$stool_gap_log <- log1p(m$stool_ultrasound_gap_days)
    concentration <- num("dna_pitoisuus")
    if (any(concentration <= 0, na.rm = TRUE)) stop("Non-positive DNA concentration")
    m$log_DNA <- log10(concentration); m$log_total <- log10(m$matrix_total)
    nsrc <- fl_cv_coldata(nutrient_path)
    if (!"SampleID" %in% names(nsrc) || anyDuplicated(nsrc$SampleID)) stop("Nutrient ID key invalid")
    j <- match(m$sample_id, as.character(nsrc$SampleID)); if (anyNA(j)) stop("Missing nutrient source identities")
    m$fiber <- as.numeric(nsrc$FIBC18[j]); m$sugar <- as.numeric(nsrc$SUGAR18[j]); kcal <- as.numeric(nsrc$ENERCK18[j])
    if (any(kcal <= 0, na.rm = TRUE) || any(m$fiber < 0, na.rm = TRUE) || any(m$sugar < 0, na.rm = TRUE)) stop("Invalid nutrient units/range")
    m$energy_log <- log10(kcal)
    # Medication harmonization is checked against the raw SAS, whose dictionary uses 1=no/2=yes.
    sas_csv <- file.path(private, "raw_medication_check.csv")
    fields <- "tutknom18,DIALAA18,FHLAAK18,VPLAAK18,REULAAK18"
    code <- file.path(cfg$project_root, "scripts/export_sas_for_r.py")
    status <- system2(cfg$python, c(shQuote(code), "--source", shQuote(cfg$clinical_file), "--out", shQuote(sas_csv), "--fields", fields))
    if (status != 0L) stop("Raw medication coding check failed")
    sas <- utils::read.csv(sas_csv, fileEncoding = "UTF-8-BOM", check.names = FALSE)
    names(sas) <- tolower(names(sas)); k <- match(as.character(m$participant_id), as.character(sas$tutknom18))
    if (anyNA(k) || anyDuplicated(sas$tutknom18)) stop("SAS medication key mismatch")
    med_audit <- list()
    for (p in list(c("diabetes_med", "dialaa18"), c("lipid_med", "fhlaak18"), c("bp_med", "vplaak18"), c("rheum_med", "reulaak18"))) {
      if (!p[2] %in% names(sas)) stop("Raw SAS medication field absent: ", p[2])
      raw <- sas[[p[2]]][k]; finite <- raw[!is.na(raw)]
      if (all(finite %in% c(1, 2))) { derived <- raw - 1; coding <- "raw SAS 1=no/2=yes -> 0/1" }
      else if (all(finite %in% c(0, 1))) { derived <- raw; coding <- "raw SAS already binary 0/1" }
      else stop("Unexpected raw SAS medication coding")
      both <- !is.na(derived) & !is.na(m[[p[1]]]); mismatch <- sum(derived[both] != m[[p[1]]][both])
      if (!sum(both)) stop("No paired source values to validate medication: ", p[1])
      med_audit[[p[1]]] <- data.frame(variable = p[1], coding = coding, compared_n = sum(both), mismatches = mismatch)
      if (mismatch) stop("Medication coding conflict: ", p[1])
    }
    fl_write_csv(do.call(rbind, med_audit), file.path(audit, "medication_coding_validation.csv"))
    rules <- fl_cv_rules(); fl_write_csv(rules, file.path(pub, "covariate_definitions.csv"))
    missing <- do.call(rbind, lapply(rules$variable, function(v) do.call(rbind,
      lapply(c("G0", "G1", "G2", "Adult", "Minor"), function(g) {
        j <- if (g == "Adult") m$age >= 18 else if (g == "Minor") m$age < 18 else m$generation == g
        data.frame(variable = v, group = g, n = sum(j), missing_n = sum(is.na(m[[v]][j])),
          available_cases = sum(m$fatty_liver[j & !is.na(m[[v]])]),
          distinct_nonmissing = length(unique(m[[v]][j & !is.na(m[[v]])])))
      }))))
    fl_write_csv(missing, file.path(pub, "covariate_missingness.csv"))
    sets <- fl_cv_scenarios(m)
    flow <- do.call(rbind, lapply(names(sets), function(s) {
      z <- sets[[s]]; j <- z$mask
      data.frame(scenario = s, n = sum(j), cases = sum(m$fatty_liver[j]), minors = sum(m$age[j] < 18),
        families = length(unique(m$family_id[j])), excluded_n = sum(!j), excluded_cases = sum(m$fatty_liver[!j]),
        covariates = paste(z$covariates, collapse = ";"), label = z$label)
    }))
    fl_write_csv(flow, file.path(pub, "scenario_sample_counts.csv"))
    raw_decisions <- data.frame(raw_field = names(src), role = "not_selected", reason = "Not a prespecified cause/measurement covariate for these estimands; availability alone is not justification")
    for (i in seq_len(nrow(rules))) {
      j <- match(tolower(rules$raw_field[i]), tolower(raw_decisions$raw_field))
      if (!is.na(j)) { raw_decisions$role[j] <- rules$role[i]; raw_decisions$reason[j] <- rules$definition[i] }
    }
    ids <- grepl("id$|^tutk|perhe|sample|^sp$|sukupolvi", raw_decisions$raw_field, ignore.case = TRUE)
    raw_decisions$reason[ids] <- "Identity, linkage or generation variable; used for alignment/stratification, never a continuous confounder; perhe2/3 excluded by user"
    fl_write_csv(raw_decisions, file.path(pub, "all_raw_field_decisions.csv"))
    module_review <- data.frame(module = sprintf("%02d", 1:8),
      action = c(rep("No regression adjustment: descriptive, flow, QC or preprocessing; cohort and abundance unchanged", 5),
        "Re-run clinical models with observed-older age/sex composition, geography and nested adult lifestyle/adiposity/medication/nutrition models",
        "Retain total IF-BF similarity; add double-endpoint residual MRQAP adjusted and perturbation analyses",
        "Re-run Shannon GEE and partial dbRDA under matched core/adult/nutrient/perturbation/DNA scenarios"))
    fl_write_csv(module_review, file.path(pub, "analysis_module_review.csv"))
    source_files <- c(primary_path, nutrient_path, cfg$clinical_file,
      file.path(cfg$project_root, "data/processed", microbiome_version, "microbiome_processed.rds"))
    fingerprints <- data.frame(path = source_files, md5 = unname(tools::md5sum(source_files)))
    fl_write_csv(fingerprints, file.path(audit, "source_fingerprints.csv"))
    version <- cfg$run_id; folder <- file.path(cfg$project_root, "data/covariates", version)
    if (file.exists(folder)) stop("Covariate dataset version already exists")
    dir.create(folder, recursive = TRUE)
    pkg <- list(schema_version = 1L, metadata = m, version = version, microbiome_version = microbiome_version,
      rules = rules, source_fingerprints = fingerprints, scope = "fixed three-generation primary families; no imputation")
    path <- file.path(folder, "covariates.rds"); saveRDS(pkg, path, version = 3)
    fl_write_csv(data.frame(file = "covariates.rds", md5 = unname(tools::md5sum(path))), file.path(folder, "dataset_manifest.csv"))
    writeLines("COMPLETED", file.path(folder, "STATUS"))
    writeLines(version, file.path(cfg$project_root, "data/covariates/CURRENT"))
    capture.output(sessionInfo(), file = file.path(audit, "sessionInfo.txt"))
    writeLines(c("COMPLETED", "Covariate inventory and immutable metadata package; no association tests."), file.path(context$run_dir, "run_status.txt"))
    final <- fl_output_complete_run(context)
    cat("Covariate dataset published: ", folder, "\nAudit: ", final, "\n", sep = "")
    invisible(pkg)
  }, error = function(e) { fl_output_fail_run(context, conditionMessage(e)); stop(e) })
}

fl_cv_load <- function(cfg, base_metadata) {
  root <- file.path(cfg$project_root, "data/covariates")
  pointer <- file.path(root, "CURRENT"); if (!file.exists(pointer)) stop("Run scripts/00_prepare_confounders.R once to build audited covariates")
  version <- readLines(pointer, n = 1)
  if (!grepl("^[A-Za-z0-9_.-]+$", version)) stop("Invalid covariate version")
  folder <- file.path(root, version); path <- file.path(folder, "covariates.rds")
  if (readLines(file.path(folder, "STATUS"), n = 1) != "COMPLETED") stop("Covariates not complete")
  manifest <- read.csv(file.path(folder, "dataset_manifest.csv"), fileEncoding = "UTF-8-BOM")
  if (unname(tools::md5sum(path)) != manifest$md5[manifest$file == "covariates.rds"]) stop("Covariate checksum mismatch")
  d <- readRDS(path); m <- d$metadata
  j <- match(base_metadata$sample_id, m$sample_id)
  if (anyNA(j) || anyDuplicated(m$sample_id) || nrow(m) != nrow(base_metadata)) stop("Covariate cohort mismatch")
  m <- m[j, , drop = FALSE]
  m$seq_pool <- factor(m$seq_pool)
  for (v in c("sample_id", "participant_id", "family_id", "generation", "fatty_liver", "age", "sex"))
    if (!identical(as.character(base_metadata[[v]]), as.character(m[[v]]))) stop("Covariates altered protected field: ", v)
  if (any(rowSums(table(m$family_id, m$generation) > 0) != 3L)) stop("Covariate baseline scope changed")
  list(metadata = m, version = version, source_files = c(path, file.path(folder, "dataset_manifest.csv")))
}
