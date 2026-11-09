# User-approved clinical module; no analysis is executed by sourcing this file.
fl_cg_config <- function(project_root) {
  cfg <- fl_config(project_root)
  cfg$run_ids["cross_generation_clinical"] <- "clinical_approved_design_v1"
  cfg$cg_review_nutrition <- FALSE
  # Follow CURRENT (or set an explicit saved version). Only metadata are used;
  # the three-generation restriction is always applied before clinical models.
  cfg$cg_data_version <- NULL
  cfg$cg_pairs <- list(c("G0", "G1"), c("G1", "G2"), c("G0", "G2"))
  cfg$cg_working_correlation <- "independence"
  cfg$cg_alpha <- 0.05
  cfg$cg_max_iterations <- 100L
  cfg$cg_export_private <- FALSE
  cfg
}
