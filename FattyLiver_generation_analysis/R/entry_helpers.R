# Shared input/audit utilities. Each analysis script has its own execution entry.
fl_entry_options <- function(args, cfg, analysis) {
  known <- args %in% c("--check", "--run", "--self-test", "--help") | grepl("^--run-id=", args)
  if (any(!known)) stop("Unknown arguments: ", paste(args[!known], collapse = ", "))
  modes <- intersect(args, c("--check", "--run", "--self-test", "--help"))
  if (length(modes) > 1L) stop("Choose only one of --check, --run, --self-test or --help")
  overrides <- args[grepl("^--run-id=", args)]
  if (length(overrides) > 1L) stop("Supply --run-id only once")
  cfg$run_id <- if (length(overrides)) sub("^--run-id=", "", overrides) else unname(cfg$run_ids[analysis])
  if (length(cfg$run_id) != 1L || is.na(cfg$run_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", cfg$run_id)) stop("Unsafe or empty run_id")
  cfg$run_id_explicit <- length(overrides) > 0L
  cfg$output_module <- unname(cfg$output_modules[analysis])
  if (length(cfg$output_module) != 1L || is.na(cfg$output_module) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", cfg$output_module))
    stop("Missing or invalid output module for ", analysis)
  list(mode = if (length(modes)) modes else "--run", cfg = cfg)
}

# Each analysis owns one numbered module. Publication only follows a successful
# run; existing adopted results survive a failed rerun. Old runs are never deleted.
fl_output_start_run <- function(cfg) {
  module <- file.path(cfg$project_root, "outputs", cfg$output_module)
  used <- function(id) any(file.exists(file.path(module, c("current", "history", "_running"), id)))
  if (used(cfg$run_id)) {
    if (isTRUE(cfg$run_id_explicit)) stop("Run ID already exists in this module: ", cfg$run_id)
    stem <- paste0(cfg$run_id, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    cfg$run_id <- stem; i <- 0L
    while (used(cfg$run_id)) { i <- i + 1L; cfg$run_id <- paste0(stem, "_", i) }
  }
  for (folder in c("current", "history", "_running"))
    dir.create(file.path(module, folder), recursive = TRUE, showWarnings = FALSE)
  run_dir <- file.path(module, "_running", cfg$run_id)
  if (!dir.create(run_dir)) stop("Cannot create output run: ", run_dir)
  list(cfg = cfg, module = module, run_dir = run_dir)
}

fl_output_complete_run <- function(context) {
  marker <- file.path(context$run_dir, "run_status.txt")
  if (!file.exists(marker) || !identical(readLines(marker, n = 1L, warn = FALSE), "COMPLETED"))
    stop("Cannot adopt an incomplete run")
  current <- file.path(context$module, "current")
  previous <- list.dirs(current, recursive = FALSE, full.names = TRUE)
  target <- file.path(current, context$cfg$run_id)
  if (file.exists(target)) stop("Current output target already exists")
  archived <- character()
  tryCatch({
    for (old in previous) {
      destination <- file.path(context$module, "history", basename(old))
      if (file.exists(destination)) stop("History collision; prior result retained: ", destination)
      if (!file.rename(old, destination)) stop("Cannot archive prior current result")
      archived <- c(archived, destination)
    }
    if (!file.rename(context$run_dir, target)) stop("Cannot publish completed run")
  }, error = function(e) {
    for (old in archived)
      if (!file.rename(old, file.path(current, basename(old))))
        warning("Previous result remains recoverable in history: ", old)
    stop(e)
  })
  # The directory move is the publication boundary. A pointer-writing problem
  # must not mark already-published results as failed.
  tryCatch(writeLines(paste0("current/", context$cfg$run_id),
    file.path(context$module, "CURRENT.txt")), error = function(e) warning(conditionMessage(e)))
  # All present/future modules using this publication helper refresh the viewer.
  # A viewer failure never changes the successful status of the analysis itself.
  gallery_code <- file.path(context$cfg$project_root, "R", "result_gallery.R")
  if (!isTRUE(context$cfg$synthetic) && file.exists(gallery_code)) {
    tryCatch({
      gallery_ctype <- Sys.getlocale("LC_CTYPE")
      on.exit(Sys.setlocale("LC_CTYPE", gallery_ctype), add = TRUE)
      if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
        if (!nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))) stop("UTF-8 locale required for result viewer")
      }
      gallery_env <- new.env(parent = baseenv())
      sys.source(gallery_code, envir = gallery_env)
      gallery_env$fl_refresh_result_gallery(context$cfg$project_root)
    }, error = function(e) warning("Analysis completed; result viewer refresh failed: ", conditionMessage(e)))
  }
  target
}

fl_output_fail_run <- function(context, message) {
  if (!dir.exists(context$run_dir)) return(invisible(NULL))
  writeLines(c("FAILED", message), file.path(context$run_dir, "run_status.txt"))
  target <- file.path(context$module, "history", context$cfg$run_id)
  if (file.exists(target) || !file.rename(context$run_dir, target))
    warning("Failed run remains in _running for inspection: ", context$run_dir)
  invisible(NULL)
}

fl_entry_help <- function(script, entry) {
  cat("Rscript scripts/", script, " [--check | --run | --self-test] [--run-id=NAME]\n",
    "Default/--run: read inputs and generate this analysis.\n",
    "--check: syntax, input paths and runtime only; no clinical data read.\n",
    "--self-test: synthetic data and this analysis only.\n",
    "RStudio: Source the script to generate results immediately.\n",
    "To load functions only: options(fatty_liver.autorun = FALSE) before source().\n",
    "Then call ", entry, "() to generate results.\n", sep = "")
  invisible(NULL)
}

fl_check_entry <- function(cfg, script) {
  files <- c("config.R", "R/common.R", "R/entry_helpers.R", paste0("scripts/", script), "tests/self_test.R")
  for (f in files) invisible(parse(file.path(cfg$project_root, f), encoding = "UTF-8"))
  fl_required_inputs(cfg)
  if (cfg$clinical_format == "sas") {
    if (!file.exists(cfg$python)) stop("Configured Python executable is absent")
    helper <- file.path(cfg$project_root, "scripts", "export_sas_for_r.py")
    if (!file.exists(helper)) stop("SAS import helper is absent")
    status <- system2(cfg$python, c("-c", shQuote("import pandas; print('pandas '+pandas.__version__)")))
    if (status != 0L) stop("Python/pandas unavailable")
  }
  cat("CHECK PASSED: ", script, "; R syntax, configured paths and runtime.\n",
    "Clinical data were not read; no research figures/tables were generated.\n", sep = "")
  invisible(cfg)
}

fl_required_inputs <- function(cfg) {
  required <- c(clinical = cfg$clinical_file, manifest = cfg$manifest_file,
    if (cfg$require_taxonomy_qc) c(qc = cfg$qc_file))
  if (any(!file.exists(required))) stop("Missing inputs: ", paste(names(required)[!file.exists(required)], collapse = ", "))
  required
}

fl_write_run_audit <- function(prepared, cfg, run_dir, required) {
  audit <- file.path(run_dir, "audit")
  fl_write_csv(prepared$flow, file.path(audit, "sample_flow.csv"))
  fl_write_csv(prepared$exclusions, file.path(audit, "sample_exclusions.csv"))
  fl_write_csv(prepared$family_eligibility, file.path(audit, "three_generation_family_eligibility.csv"))
  fl_write_csv(data.frame(rule = prepared$cohort_contract$rule,
    family_field = prepared$cohort_contract$family_field, basis = prepared$cohort_contract$basis,
    families = prepared$cohort_contract$families, participants = prepared$cohort_contract$participants),
    file.path(audit, "cohort_contract.csv"))
  fl_write_csv(prepared$field_audit, file.path(audit, "field_source_audit.csv"))
  fl_write_csv(fl_summary_counts(prepared$cohort), file.path(audit, "cohort_counts.csv"))
  age_source <- as.data.frame(table(prepared$cohort$generation, prepared$cohort$age_source, useNA = "ifany"))
  names(age_source) <- c("generation", "age_source", "n")
  fl_write_csv(age_source, file.path(audit, "age_source_audit.csv"))
  families <- do.call(rbind, lapply(fl_groups(prepared$cohort), function(d) data.frame(
    participants = nrow(d), families_with_known_id = length(unique(d$family_id[!is.na(d$family_id)])),
    participants_without_family_id = sum(is.na(d$family_id)))))
  families$group <- rownames(families); rownames(families) <- NULL
  fl_write_csv(families, file.path(audit, "family_counts.csv"))
  source_paths <- unique(c(required, if (file.exists(cfg$family_file)) cfg$family_file))
  fingerprints <- data.frame(source = source_paths, md5 = unname(tools::md5sum(source_paths)))
  fl_write_csv(fingerprints, file.path(audit, "input_fingerprints.csv"))
  if (cfg$export_participant_level_cohort) {
    dir.create(file.path(run_dir, "private"), showWarnings = FALSE)
    saveRDS(prepared$cohort, file.path(run_dir, "private", "analysis_cohort.rds"))
    fl_write_csv(prepared$full_inclusion, file.path(run_dir, "private", "participant_inclusion.csv"))
  }
  invisible(NULL)
}

fl_with_prepared_run <- function(cfg, output_folder, build, prepare = NULL) {
  required <- fl_required_inputs(cfg)
  context <- fl_output_start_run(cfg); cfg <- context$cfg; run_dir <- context$run_dir
  for (name in c(output_folder, "audit")) dir.create(file.path(run_dir, name))
  writeLines(c("RUNNING", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), file.path(run_dir, "run_status.txt"))
  label <- if (!is.null(cfg$analysis_label)) cfg$analysis_label else if (output_folder == "tables") "Table 1" else "Figure 1"
  cat("Starting ", label,
    "\nOutput directory: ", normalizePath(run_dir, winslash = "/"), "\n", sep = "")
  tryCatch({
    input <- fl_read_inputs(cfg, run_dir)
    fl_write_csv(input$duplicate_audit, file.path(run_dir, "audit", "input_duplicate_columns.csv"))
    prepared <- fl_build_cohort(input$clinical, input$manifest, input$qc, input$family, cfg)
    if (!is.null(prepare)) prepared <- prepare(prepared, input, cfg)
    cat("Eligible participants: ", nrow(prepared$cohort), "\n", sep = "")
    fl_write_run_audit(prepared, cfg, run_dir, required)
    build(prepared, file.path(run_dir, output_folder), cfg)
    capture.output(utils::sessionInfo(), file = file.path(run_dir, "audit", "sessionInfo.txt"))
    note <- if (!is.null(cfg$completion_note)) cfg$completion_note else
      "Descriptive outputs only; no association models were fitted."
    writeLines(c("COMPLETED", note,
      if (output_folder == "figures") "Inspect the final-size Figure 1 before manuscript use."),
      file.path(run_dir, "run_status.txt"))
    final_dir <- fl_output_complete_run(context)
    cat("Completed: ", normalizePath(final_dir, winslash = "/"), "\n", sep = "")
    invisible(prepared)
  }, error = function(e) {
    fl_output_fail_run(context, conditionMessage(e))
    stop(e)
  })
}
