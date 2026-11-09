fl_text <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "N/A", "NULL", ".")] <- NA_character_
  x
}
# One outcome-blind family rule shared by all clinical and microbiome entries.
# Apply after baseline QC/generation/ultrasound eligibility, before model-specific
# complete-case restrictions; later restrictions must not redefine family membership.
fl_primary_family_key <- function(v) {
  v <- fl_text(v)
  numeric_label <- !is.na(v) & grepl("^[0-9]+(\\.0+)?$", v)
  v[numeric_label] <- sub("\\.0+$", "", v[numeric_label])
  v[numeric_label] <- sub("^0+(?=[0-9])", "", v[numeric_label], perl = TRUE)
  v
}

fl_three_generation_scope <- function(x) {
  if (!all(c("family_id", "generation") %in% names(x))) stop("Family/generation fields missing")
  f <- fl_primary_family_key(x$family_id)
  if ("tutkno80perhe" %in% names(x)) {
    raw <- fl_primary_family_key(x$tutkno80perhe)
    both <- !is.na(raw) & !is.na(f)
    if (any(raw[both] != f[both])) stop("Primary family labels conflict")
    f <- raw # A mapped label never fills a missing primary clinical label.
  }
  g <- as.character(x$generation)
  if (anyNA(g) || any(!g %in% c("G0", "G1", "G2"))) stop("Invalid generation for family eligibility")
  labels <- sort(unique(f[!is.na(f)]))
  coverage <- table(factor(f, levels = labels), factor(g, levels = c("G0", "G1", "G2")))
  eligible <- labels[rowSums(coverage > 0L) == 3L]
  keep <- !is.na(f) & f %in% eligible
  if (!any(keep)) stop("No eligible families with G0, G1 and G2 after baseline QC and ultrasound eligibility")
  x$family_id <- f
  excluded_missing <- is.na(f)
  excluded_incomplete <- !is.na(f) & !keep
  summary <- data.frame(stage = c("Baseline ultrasound-eligible participants",
    "Excluded: missing primary family label", "Excluded: family lacks one or more generations",
    "Three-generation family cohort"),
    n = c(nrow(x), sum(excluded_missing), sum(excluded_incomplete), sum(keep)),
    families = c(length(labels), NA_integer_, length(labels) - length(eligible), length(eligible)))
  list(cohort = x[keep, , drop = FALSE], keep = keep,
    excluded_missing = excluded_missing, excluded_incomplete = excluded_incomplete,
    summary = summary,
    contract = list(rule = "G0_G1_G2_each_at_least_one", family_field = "tutkno80perhe",
      basis = "microbiome-linked baseline QC and valid ultrasound; before covariate restrictions",
      required_generations = c("G0", "G1", "G2"), families = length(eligible), participants = sum(keep)))
}

fl_column <- function(x, name, required = FALSE) {
  hit <- which(tolower(names(x)) == tolower(name))
  if (length(hit) > 1L) stop("Ambiguous column: ", name)
  if (!length(hit)) {
    if (required) stop("Required column missing: ", name)
    return(rep(NA_character_, nrow(x)))
  }
  x[[hit]]
}
fl_has_column <- function(x, name) any(tolower(names(x)) == tolower(name))
fl_number <- function(x, field = "value", strict = FALSE) {
  raw <- fl_text(x); out <- suppressWarnings(as.numeric(raw))
  bad <- !is.na(raw) & !is.finite(out)
  if (any(bad)) {
    msg <- paste(field, "has", sum(bad), "unparseable/non-finite values")
    if (strict) stop(msg) else warning(msg, call. = FALSE)
    out[bad] <- NA_real_
  }
  out
}
fl_unique_id <- function(x, label) {
  x <- fl_text(x)
  if (anyNA(x) || anyDuplicated(x)) stop(label, " must be nonmissing and unique")
  x
}
fl_read_tsv <- function(path, wanted, encoding = "UTF-8") {
  if (!file.exists(path)) stop("Input not found: ", path)
  header <- strsplit(readLines(path, n = 1L, encoding = encoding), "\t", fixed = TRUE)[[1L]]
  header <- sub("^\ufeff", "", header); header <- sub('^"(.*)"$', "\\1", header)
  keep <- tolower(header) %in% tolower(wanted)
  if (!any(keep)) stop("No requested columns in: ", path)
  # Mark string encodings without transcoding the entire wide file to the
  # Windows native locale, which can truncate otherwise valid Latin-1 input.
  x <- utils::read.delim(path, colClasses = ifelse(keep, "character", "NULL"),
    encoding = encoding, quote = "\"", comment.char = "", check.names = FALSE,
    na.strings = c("", "NA"), stringsAsFactors = FALSE)
  duplicates <- unique(tolower(names(x))[duplicated(tolower(names(x)))])
  audit <- data.frame(source = character(), field = character(), original_names = character(),
    copies = integer(), checked_rows = integer(), action = character())
  for (field in duplicates) {
    positions <- which(tolower(names(x)) == field)
    reference <- fl_text(x[[positions[1L]]])
    equal <- vapply(positions[-1L], function(i) identical(reference, fl_text(x[[i]])), logical(1))
    if (!all(equal)) stop("Conflicting duplicate TSV columns for ", field, " in ", path)
    audit <- rbind(audit, data.frame(source = path, field = field,
      original_names = paste(names(x)[positions], collapse = " | "), copies = length(positions),
      checked_rows = nrow(x), action = "Identical values and missingness; retained first column"))
  }
  x <- x[, !duplicated(tolower(names(x))), drop = FALSE]
  attr(x, "duplicate_column_audit") <- audit
  x
}
fl_date <- function(x, field = "date") {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x, tz = "UTC"))
  raw <- fl_text(x); out <- as.Date(rep(NA_character_, length(raw)))
  iso <- !is.na(raw) & grepl("^\\d{4}-\\d{2}-\\d{2}", raw)
  out[iso] <- as.Date(substr(raw[iso], 1L, 10L), format = "%Y-%m-%d")
  numeric <- !is.na(raw) & grepl("^-?[0-9]+(\\.[0-9]+)?$", raw)
  if (any(numeric)) {
    n <- as.numeric(raw[numeric])
    out[numeric] <- as.Date(ifelse(abs(n) > 1e6, n / 86400, n), origin = "1960-01-01")
  }
  if (any(!is.na(raw) & is.na(out))) warning("Unparsed dates in ", field, call. = FALSE)
  out
}
fl_birthday <- function(birth, year) {
  out <- as.Date(sprintf("%04d-%s", as.integer(year), format(birth, "%m-%d")), format = "%Y-%m-%d")
  leap <- !is.na(birth) & format(birth, "%m-%d") == "02-29" & is.na(out)
  out[leap] <- as.Date(sprintf("%04d-02-28", as.integer(year[leap])))
  out
}
fl_age <- function(birth, date) {
  out <- rep(NA_real_, length(birth)); valid <- !is.na(birth) & !is.na(date) & date >= birth
  b <- birth[valid]; d <- date[valid]; year <- as.integer(format(d, "%Y"))
  year <- year - as.integer(d < fl_birthday(b, year))
  previous <- fl_birthday(b, year); following <- fl_birthday(b, year + 1L)
  out[valid] <- year - as.integer(format(b, "%Y")) + as.numeric(d - previous) / as.numeric(following - previous)
  out[!is.finite(out) | out < 0 | out > 120] <- NA_real_; out
}
fl_write_csv <- function(x, path) {
  # Stream large permutation audits instead of accumulating hundreds of
  # thousands of lines in capture.output's growing text connection.
  con <- file(path, open = "wt", encoding = "UTF-8"); on.exit(close(con), add = TRUE)
  writeLines("\ufeff", con, sep = "", useBytes = TRUE)
  utils::write.table(x, file = con, sep = ",", row.names = FALSE,
    col.names = TRUE, na = "", qmethod = "double")
  invisible(path)
}
fl_json <- function(x) {
  escape <- function(s) {
    s <- gsub("\\", "\\\\", s, fixed = TRUE); s <- gsub('"', '\\"', s, fixed = TRUE)
    s <- gsub("\n", "\\n", s, fixed = TRUE); s <- gsub("\r", "\\r", s, fixed = TRUE)
    s <- gsub("\t", "\\t", s, fixed = TRUE); paste0('"', s, '"')
  }
  if (is.null(x)) return("null")
  if (is.list(x)) {
    parts <- vapply(x, fl_json, character(1))
    if (!is.null(names(x)) && all(nzchar(names(x)))) return(paste0("{", paste(paste0(escape(names(x)), ":", parts), collapse = ","), "}"))
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  if (length(x) != 1L) return(paste0("[", paste(vapply(as.list(x), fl_json, character(1)), collapse = ","), "]"))
  if (is.na(x)) return("null")
  if (is.logical(x)) return(if (x) "true" else "false")
  if (is.numeric(x)) return(if (is.finite(x)) format(x, scientific = FALSE, trim = TRUE, digits = 15) else "null")
  escape(as.character(x))
}
fl_read_inputs <- function(cfg, run_dir) {
  fields <- unique(c(fl_requested_fields(), cfg$extra_requested_fields))
  if (identical(cfg$clinical_format, "sas")) {
    if (!file.exists(cfg$python)) stop("Python executable not found: ", cfg$python)
    cache <- file.path(run_dir, "private"); dir.create(cache, recursive = TRUE)
    clinical_csv <- file.path(cache, "clinical_selected.csv")
    helper <- file.path(cfg$project_root, "scripts", "export_sas_for_r.py")
    status <- system2(cfg$python, c(shQuote(helper), "--source", shQuote(cfg$clinical_file),
      "--out", shQuote(clinical_csv), "--fields", shQuote(paste(fields, collapse = ","))))
    if (status != 0L) stop("SAS export failed; no clinical analysis was performed")
  } else if (identical(cfg$clinical_format, "csv")) clinical_csv <- cfg$clinical_file
  else stop("clinical_format must be sas or csv")
  clinical <- utils::read.csv(clinical_csv, colClasses = "character", check.names = FALSE,
    fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
  manifest <- fl_read_tsv(cfg$manifest_file, unique(c("sample_id", "SampleID", fields,
    "prep_batch", "seq_pool")), cfg$manifest_encoding)
  qc <- if (cfg$require_taxonomy_qc) fl_read_tsv(cfg$qc_file,
    c("sample_id", "total", "detected_features", "missing_fraction")) else NULL
  family <- if (file.exists(cfg$family_file)) fl_read_tsv(cfg$family_file, c("sample_id", "family_id")) else NULL
  duplicate_audit <- do.call(rbind, lapply(list(manifest, qc, family), attr, which = "duplicate_column_audit"))
  list(clinical = clinical, manifest = manifest, qc = qc, family = family, duplicate_audit = duplicate_audit)
}
fl_build_cohort <- function(clinical, manifest, qc, family, cfg) {
  cid <- fl_unique_id(fl_column(clinical, "tutknom18", TRUE), "Clinical participant ID")
  mid <- fl_unique_id(fl_column(manifest, "tutknom18", TRUE), "Microbiome participant ID")
  sid <- fl_unique_id(fl_column(manifest, "sample_id", TRUE), "Microbiome sample ID")
  index <- match(mid, cid); matched <- !is.na(index)
  x <- data.frame(sample_id = sid, participant_id = mid, stringsAsFactors = FALSE); audit <- list()
  for (field in setdiff(fl_requested_fields(), "tutknom18")) {
    has_c <- fl_has_column(clinical, field); has_m <- fl_has_column(manifest, field)
    cval <- fl_text(fl_column(clinical, field))[index]; mval <- fl_text(fl_column(manifest, field))
    # Metadata fills only entirely absent columns, never clinical missing entries.
    x[[field]] <- if (has_c) cval else mval
    comparable <- !is.na(cval) & !is.na(mval)
    if (grepl("pvm", field, ignore.case = TRUE) || field == "syntymapvm") {
      conflict <- comparable & as.character(fl_date(cval)) != as.character(fl_date(mval))
    } else {
      cn <- suppressWarnings(as.numeric(cval)); mn <- suppressWarnings(as.numeric(mval))
      both <- comparable & is.finite(cn) & is.finite(mn); conflict <- comparable & cval != mval
      conflict[both] <- abs(cn[both] - mn[both]) > 1e-6
    }
    nc <- sum(conflict, na.rm = TRUE)
    audit[[field]] <- data.frame(field = field,
      source = if (has_c) "clinical" else if (has_m) "manifest" else "not available",
      comparable_n = sum(comparable), conflicting_n = nc, selected_nonmissing_n = sum(!is.na(x[[field]])))
    if (field %in% c("sukupolvi", "sp", "RASVAMAKSA18") && nc > 0L) stop("Critical source conflict in ", field, ": ", nc)
  }
  x$generation <- toupper(fl_text(x$sukupolvi)); x$family_id <- fl_primary_family_key(x$tutkno80perhe)
  if (!is.null(family)) {
    fid <- fl_unique_id(fl_column(family, "sample_id", TRUE), "Family-map sample ID")
    f <- fl_primary_family_key(fl_column(family, "family_id", TRUE))[match(sid, fid)]
    both <- !is.na(f) & !is.na(x$family_id)
    if (any(f[both] != x$family_id[both])) stop("Primary clinical and mapped family labels disagree")
  }
  x$ultrasound_date <- fl_date(x$MAKSAPVM18, "MAKSAPVM18"); x$birth_date <- fl_date(x$syntymapvm, "syntymapvm")
  x$stool_date <- fl_date(x$ULOSTE_PVM18, "ULOSTE_PVM18"); x$age <- fl_age(x$birth_date, x$ultrasound_date)
  x$age_source <- ifelse(is.na(x$age), NA_character_, "birth date + liver ultrasound date")
  for (fallback in c("ika_status18", "TARKKAIKA18")) {
    a <- fl_number(x[[fallback]], fallback); use <- is.na(x$age) & is.finite(a) & a >= 0 & a <= 120
    x$age[use] <- a[use]; x$age_source[use] <- fallback
  }
  # ika18 is age reached in calendar year 2018 and is deliberately not used.
  x$age_group <- ifelse(is.na(x$age), "Age unavailable", ifelse(x$age < 18, "Minor", "Adult"))
  x$stool_ultrasound_gap_days <- abs(as.numeric(x$stool_date - x$ultrasound_date))
  x$fat_code <- fl_number(x$RASVAMAKSA18, "RASVAMAKSA18", strict = TRUE)
  if (any(!is.na(x$fat_code) & !x$fat_code %in% 1:5)) stop("Unrecognised ultrasound code")
  x$fatty_liver <- ifelse(x$fat_code %in% 1:3, as.integer(x$fat_code %in% c(2, 3)), NA_integer_)
  numeric_fields <- setdiff(fl_requested_fields(), c("tutknom18", "sukupolvi", "tutkno80perhe",
    "syntymapvm", "MAKSAPVM18", "TUTKPVM18", "ULOSTE_PVM18"))
  for (field in numeric_fields) x[[field]] <- fl_number(x[[field]], field)
  if (any(!is.na(x$sp) & !x$sp %in% c(1, 2))) stop("Unrecognised sex code")
  known_codes <- list(NOWOB_CAT18 = 0:2, smoke18 = 0:1, ANTBV18 = 0:1,
    VPLAAK18 = 0:1, FHLAAK18 = 0:1, DIALAA18 = 0:1)
  for (field in names(known_codes)) {
    values <- x[[field]]
    if (any(!is.na(values) & !values %in% known_codes[[field]]))
      stop("Unexpected coding in ", field, "; inspect the dictionary and transformation provenance")
  }
  x$sex <- ifelse(is.na(x$sp), NA_character_, ifelse(x$sp == 1, "Female", "Male"))
  qc_ok <- rep(TRUE, nrow(x))
  if (cfg$require_taxonomy_qc) {
    qid <- fl_unique_id(fl_column(qc, "sample_id", TRUE), "QC sample ID"); qi <- match(sid, qid)
    total <- fl_number(fl_column(qc, "total", TRUE))[qi]; features <- fl_number(fl_column(qc, "detected_features", TRUE))[qi]
    missing <- fl_number(fl_column(qc, "missing_fraction", TRUE))[qi]
    qc_ok <- !is.na(qi) & is.finite(total) & total > 0 & is.finite(features) & features > 0 &
      is.finite(missing) & missing >= 0 & missing <= cfg$max_qc_missing_fraction
  }
  reason <- rep("Included", nrow(x))
  exclude <- function(condition, label) {
    use <- reason == "Included" & condition; use[is.na(use)] <- FALSE; reason[use] <<- label
  }
  exclude(!matched, "No matching clinical record"); exclude(!qc_ok, "Microbiome QC failed or unavailable")
  exclude(is.na(x$generation) | !x$generation %in% cfg$generation_levels, "Generation unavailable or invalid")
  exclude(is.na(x$fat_code), "Ultrasound outcome missing")
  exclude(!is.na(x$fat_code) & x$fat_code == 4, "Ultrasound missing (code 4)")
  exclude(!is.na(x$fat_code) & x$fat_code == 5, "Ultrasound unassessable (code 5)")
  baseline <- which(reason == "Included")
  scope <- fl_three_generation_scope(x[baseline, , drop = FALSE])
  reason[baseline[scope$excluded_missing]] <- "Primary family label missing"
  reason[baseline[scope$excluded_incomplete]] <- "Family does not cover G0/G1/G2"
  x$inclusion_reason <- reason
  flow <- data.frame(stage = c("Clinical source records", "Microbiome participants", "Matched participants",
    "After microbiome QC", "After generation check", "Final ultrasound cohort", "Three-generation family cohort"),
    n = c(nrow(clinical), nrow(manifest), sum(matched), sum(matched & qc_ok),
      sum(matched & qc_ok & x$generation %in% cfg$generation_levels), length(baseline), sum(reason == "Included")))
  exclusions <- as.data.frame(table(factor(reason, levels = c("No matching clinical record", "Microbiome QC failed or unavailable",
    "Generation unavailable or invalid", "Ultrasound outcome missing", "Ultrasound missing (code 4)",
    "Ultrasound unassessable (code 5)", "Primary family label missing",
    "Family does not cover G0/G1/G2", "Included"))), stringsAsFactors = FALSE)
  names(exclusions) <- c("reason", "n")
  exclusions <- rbind(data.frame(reason = "Clinical records without microbiome (outside target cohort)", n = sum(!cid %in% mid)), exclusions)
  cohort <- x[reason == "Included", , drop = FALSE]
  if (!nrow(cohort)) stop("No eligible participants")
  stopifnot(sum(exclusions$n[-1L]) == nrow(manifest), all(cohort$fat_code %in% 1:3), !anyNA(cohort$fatty_liver))
  list(cohort = cohort, flow = flow, exclusions = exclusions, field_audit = do.call(rbind, audit),
    family_eligibility = scope$summary, cohort_contract = scope$contract,
    full_inclusion = x[c("sample_id", "participant_id", "inclusion_reason")])
}
fl_groups <- function(x, supplement = FALSE) {
  if (supplement) return(list(G2 = x[x$generation == "G2", , drop = FALSE],
    `G2 <18` = x[x$generation == "G2" & x$age_group == "Minor", , drop = FALSE],
    `G2 >=18` = x[x$generation == "G2" & x$age_group == "Adult", , drop = FALSE],
    `G2 age unavailable` = x[x$generation == "G2" & x$age_group == "Age unavailable", , drop = FALSE]))
  c(list(Overall = x), setNames(lapply(c("G0", "G1", "G2"), function(g) x[x$generation == g, , drop = FALSE]), c("G0", "G1", "G2")))
}
fl_summary_counts <- function(x) {
  groups <- c(fl_groups(x), fl_groups(x, TRUE)[-1L])
  do.call(rbind, lapply(names(groups), function(g) {
    d <- groups[[g]]
    data.frame(group = g, n = nrow(d), cases = sum(d$fatty_liver == 1), controls = sum(d$fatty_liver == 0),
      mild = sum(d$fat_code == 2), clear = sum(d$fat_code == 3), age_nonmissing = sum(!is.na(d$age)),
      proportion = if (nrow(d)) mean(d$fatty_liver) else NA_real_)
  }))
}
