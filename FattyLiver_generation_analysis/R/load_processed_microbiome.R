# Base-R reader: source once in an analysis script, then load the current data.
# All exported matrices are features x samples; metadata rows follow sample order.
fl_validate_processed_microbiome <- function(d, deep = FALSE, require_three_generation = TRUE) {
  if (!identical(d$schema_version, 1L) || !isTRUE(d$analysis_ready) ||
      !identical(d$readiness_scope, "counts_based_analysis")) stop("Unsupported or incomplete processed dataset")
  if (require_three_generation) {
    if (!identical(d$cohort_contract$rule, "G0_G1_G2_each_at_least_one") ||
        !identical(d$cohort_contract$family_field, "tutkno80perhe"))
      stop("Legacy dataset lacks the three-generation family restriction. Run scripts/05_prepare_microbiome.R to create the new version.")
    m <- d$metadata
    if (!all(c("family_id", "generation") %in% names(m)) || anyNA(m$family_id) ||
        any(!nzchar(trimws(m$family_id))) || anyNA(m$generation) ||
        any(!m$generation %in% c("G0", "G1", "G2"))) stop("Invalid three-generation family metadata")
    coverage <- table(m$family_id, factor(m$generation, levels = c("G0", "G1", "G2")))
    if (any(rowSums(coverage > 0L) != 3L) ||
        d$cohort_contract$participants != nrow(m) || d$cohort_contract$families != nrow(coverage))
      stop("Three-generation family contract mismatch")
  }
  validate_set <- function(s) {
    ids <- s$metadata$sample_id; features <- s$features$feature_id
    if (!length(ids) || anyNA(ids) || anyDuplicated(ids) || !length(features) ||
        anyNA(features) || anyDuplicated(features)) stop("Invalid sample/feature identifiers")
    for (nm in c("counts", "relative", "relative_closed", "log10", "clr")) {
      m <- s$matrices[[nm]]
      if (!is.matrix(m) || !is.numeric(m) || !identical(colnames(m), ids) ||
          !identical(rownames(m), features)) stop("Processed matrix alignment failure: ", nm)
      if (deep && any(!is.finite(m))) stop("Nonfinite processed matrix: ", nm)
    }
    if (deep) {
      if (any(s$matrices$counts < 0) || any(s$matrices$counts != floor(s$matrices$counts))) stop("Invalid counts")
      expected <- sweep(s$matrices$counts, 2L, s$metadata$matrix_total, "/")
      if (max(abs(expected - s$matrices$relative)) > 1e-12) stop("Relative abundance denominator mismatch")
      if (max(abs(colSums(s$matrices$relative_closed) - 1)) > 1e-10 ||
          max(abs(colMeans(s$matrices$clr))) > 1e-10) stop("Closure/CLR validation failed")
    }
    invisible(TRUE)
  }
  validate_set(d)
  if (any(!vapply(d$analysis_sets, function(i) is.logical(i) && length(i) == nrow(d$metadata) &&
      !anyNA(i), logical(1)))) stop("Invalid analysis subset masks")
  for (g in names(d$by_generation)) {
    validate_set(d$by_generation[[g]])
    if (!all(d$by_generation[[g]]$metadata$generation == g)) stop("Generation subset mismatch")
  }
  if (deep) {
    full <- d$full
    if (!identical(colnames(full$counts), d$metadata$sample_id) ||
        !identical(dimnames(full$counts), dimnames(full$relative)) ||
        any(!is.finite(full$relative)) || max(abs(colSums(full$relative) - 1)) > 1e-10)
      stop("Full nonzero-feature matrix validation failed")
    expected <- sweep(full$counts, 2L, colSums(full$counts), "/")
    if (max(abs(expected - full$relative)) > 1e-12) stop("Full relative matrix mismatch")
  }
  invisible(TRUE)
}

fl_load_processed_microbiome <- function(project_root = Sys.getenv("FATTYLIVER_ROOT",
    unset = "E:/YFS/workspace/FattyLiver_generation_analysis"), generation = NULL,
    version = NULL, verify_checksum = TRUE, require_three_generation = TRUE) {
  root <- file.path(project_root, "data", "processed")
  if (is.null(version)) {
    pointer <- file.path(root, "CURRENT")
    if (!file.exists(pointer)) stop("Processed data missing: run scripts/05_prepare_microbiome.R once")
    version <- readLines(pointer, n = 1L, warn = FALSE)
  }
  if (length(version) != 1L || !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", version)) stop("Invalid processed version")
  folder <- file.path(root, version)
  status <- file.path(folder, "STATUS")
  if (!file.exists(status) || !identical(readLines(status, n = 1L, warn = FALSE), "COMPLETED"))
    stop("Processed version is not complete: ", version)
  path <- file.path(folder, "microbiome_processed.rds")
  if (verify_checksum) {
    audit <- utils::read.csv(file.path(folder, "dataset_manifest.csv"),
      fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
    expected <- audit$md5[audit$file == "microbiome_processed.rds"]
    if (length(expected) != 1L || !identical(unname(tools::md5sum(path)), expected)) stop("Processed dataset checksum mismatch")
  }
  d <- readRDS(path)
  # FALSE is reserved for legacy metadata extraction in 06; that entry applies
  # the shared restriction before any model. Abundance analyses use TRUE.
  fl_validate_processed_microbiome(d, require_three_generation = require_three_generation)
  if (!identical(d$version, version)) stop("Stored dataset version mismatch")
  if (!is.null(generation)) {
    if (length(generation) != 1L || !generation %in% names(d$by_generation)) stop("generation must be G0, G1 or G2")
    selected <- d$by_generation[[generation]]
    selected$schema_version <- d$schema_version; selected$analysis_ready <- d$analysis_ready
    selected$readiness_scope <- d$readiness_scope; selected$input_contract <- d$input_contract
    selected$version <- d$version; selected$source_paths <- d$source_paths
    selected$cohort_contract <- d$cohort_contract
    return(selected)
  }
  d
}
