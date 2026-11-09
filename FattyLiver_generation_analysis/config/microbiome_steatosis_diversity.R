fl_dv_config <- function(project_root) {
  cfg <- fl_config(project_root)
  cfg$run_ids["microbiome_steatosis_diversity"] <- "diversity_approved_design_v1"
  cfg$dv_data_version <- NULL
  cfg$dv_seed <- 20261020L
  cfg$dv_permutations <- 9999L
  cfg$dv_permutation_chunk <- 250L
  cfg$dv_age_df <- 3L
  cfg
}
