# Sourcing configuration performs no analysis.
fl_ms_config <- function(project_root) {
  cfg <- fl_config(project_root)
  cfg$run_ids["cross_generation_microbiome"] <- "similarity_approved_design_v1"
  cfg$ms_data_version <- NULL # CURRENT, recorded and checksummed at runtime.
  cfg$ms_pairs <- list(c("G0", "G1"), c("G1", "G2"), c("G0", "G2"))
  cfg$ms_permutations <- 9999L
  cfg$ms_seed <- 20261019L
  cfg$ms_age_band_years <- 5
  cfg$ms_schemes <- c("generation_sex_age", "generation_sex_age_prep_batch")
  # Fixed before examining review estimates. Matched-sample nested scenarios.
  cfg$ms_review_scenarios <- c("core_base", "core_adjusted", "adult_core", "adult_lifestyle", "adult_BMI",
    "bowel_base", "bowel_adjusted", "unperturbed_core", "DNA_base", "core_DNA")
  cfg
}
