# Sourcing this file reads no data.
fl_config <- function(project_root) {
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  cfg <- list(
    project_root = root,
    clinical_format = "sas",
    clinical_file = "E:/YFS/staualab18_20260902.sas7bdat",
    manifest_file = paste0("E:/YFS/workspace/microbiome_analysis/outputs/",
                           "microbiome_taxonomy_v1/01_manifest/samples.tsv"),
    qc_file = paste0("E:/YFS/workspace/microbiome_analysis/outputs/",
                     "microbiome_taxonomy_v1/02_technical_qc/sample_qc.tsv"),
    family_file = paste0("E:/YFS/workspace/microbiome_analysis/outputs/",
                         "family_design_v1/analysis_samples.tsv"),
    python = paste0("C:/Users/jzqjzq/.cache/codex-runtimes/",
                    "codex-primary-runtime/dependencies/python/python.exe"),
    manifest_encoding = "latin1",
    require_taxonomy_qc = TRUE,
    max_qc_missing_fraction = 0,
    generation_levels = c("G0", "G1", "G2"),
    output_modules = c(sgb_association = "17_SGB_association", harmonized_visualization = "16_Harmonized_visualization", english_journal_tables = "15_English_journal_tables", figure_tables = "14_Figure_tables", adult_covariate_mi = "12_Adult_covariate_MI", clinical_mi = "13_Clinical_MI", nature_composites = "11_Nature_composites", results_visualization = "10_Results_visualization", technical_covariates = "09_Technical_covariates", covariate_review = "00_Covariate_review", table1 = "01_Table1", figure1 = "02_Figure1",
      basic_exploration = "03_Basic_exploration", microbiome_qc = "04_Microbiome_QC",
      processed_microbiome = "05_Microbiome_preprocessing",
      cross_generation_clinical = "06_Cross_generation_clinical",
      cross_generation_microbiome = "07_Cross_generation_microbiome_similarity",
      microbiome_steatosis_diversity = "08_Microbiome_steatosis_diversity"),
    run_ids = c(sgb_association = "sgb_associations_v1_20261022", harmonized_visualization = "harmonized_layout_greys_v4", english_journal_tables = "english_journal_tables_v1", figure_tables = "nature_figures_as_tables_v1", adult_covariate_mi = "adult_multilevel_mice_v1", clinical_mi = "adult_clinical_mi_v1", nature_composites = "nature_composites_v1", results_visualization = "unified_results_visualization_v1", covariate_review = "covariates_review_v1", table1 = "table1_three_generation_v1", figure1 = "figure1_three_generation_v1",
                basic_exploration = "basic_exploration_three_generation_v1",
                microbiome_qc = "microbiome_qc_three_generation_v1",
                processed_microbiome = "processed_microbiome_three_generation_v1",
                cross_generation_microbiome = "microbiome_family_similarity_v1",
                microbiome_steatosis_diversity = "microbiome_steatosis_diversity_v1"),
    processed_microbiome_root = file.path(root, "data", "processed"),
    microbiome_file = "E:/YFS/documents-export-2026-08-31/tse_metaphlan_2025.rds",
    taxonomy_map_file = paste0("E:/YFS/documents-export-2026-08-31/",
      "mpa_vJan25_CHOCOPhlAnSGB_202503_species.txt.bz2"),
    source_scan_root = "E:/YFS",
    qc_prevalence = 0.10,
    qc_sensitivity_prevalence = c(0.01, 0.05, 0.10, 0.20),
    qc_low_mad_multiplier = 3,
    qc_seed = 20261018L,
    export_participant_level_cohort = FALSE,
    figure_width_mm = 183,
    figure_height_mm = 160,
    figure_dpi = 600,
    figure_font = "Arial",
    generation_colors = c(G0 = "#F1B656", G1 = "#397FC7", G2 = "#040676"),
    binary_colors = c(Reference = "#79CAFB", Comparison = "#F7A24F"),
    medication_semantics_verified = FALSE,
    include_followup_diabetes = FALSE,
    synthetic = FALSE
  )
  cfg$output_modules <- c(cfg$output_modules, sgb_abundance = "18_SGB_abundance",
    sgb_abundance_visualization = "19_SGB_abundance_visualization")
  cfg$run_ids <- c(cfg$run_ids, sgb_abundance = "sgb_abundance_v3_20261022",
    sgb_abundance_visualization = "sgb_abundance_figures_v3_20261022")
  cfg$output_modules <- c(cfg$output_modules, ast_sgb = "19_AST_SGB_analysis")
  cfg$run_ids <- c(cfg$run_ids, ast_sgb = "ast_sgb_v1_20261023")
  cfg$output_modules <- c(cfg$output_modules, palette_refresh = "20_Palette_refresh")
  cfg$run_ids <- c(cfg$run_ids, palette_refresh = "binary_blue_orange_20261023")
  cfg$output_modules <- c(cfg$output_modules, ast_overall = "21_AST_overall")
  cfg$run_ids <- c(cfg$run_ids, ast_overall = "ast_overall_v1_20261023")
  cfg$output_modules <- c(cfg$output_modules, ast_volcano = "22_AST_overall_volcano")
  cfg$run_ids <- c(cfg$run_ids, ast_volcano = "ast_overall_volcano_v1_20261024")
  cfg$output_modules <- c(cfg$output_modules, log_overall = "23_MaAsLin2_LOG_overall")
  cfg$run_ids <- c(cfg$run_ids, log_overall = "maaslin2_log_overall_v1_20261024")
  cfg$output_modules <- c(cfg$output_modules, two_generation_reanalysis = "24_Two_generation_reanalysis")
  cfg$run_ids <- c(cfg$run_ids, two_generation_reanalysis = "two_generation_20261024")
  cfg$analysis_protocol <- "two_generation_C1_C2_B1_B2_P1_P2_20261024"
  cfg$current_analysis_entry <- "scripts/24_two_generation_pipeline.R"
  cfg$output_modules <- c(cfg$output_modules,
    taxa_robustness = "25_Taxa_robustness",
    taxon_family_sharing = "26_Taxon_family_sharing")
  cfg$run_ids <- c(cfg$run_ids,
    taxa_robustness = "taxa_effect_comparisons_v2_20261026",
    taxon_family_sharing = "taxon_family_regression_v2_20261026")
  cfg$output_modules <- c(cfg$output_modules,
    family_steatosis_microbiome = "27_Family_steatosis_microbiome")
  cfg$run_ids <- c(cfg$run_ids,
    family_steatosis_microbiome = "unaffected_younger_family_exposure_v1_20261027")
  cfg$output_modules <- c(cfg$output_modules, abc_genus = "28_Family_genus_community",
    abc_within = "29_Within_family_steatosis", abc_external = "30_External_health_sets")
  cfg$run_ids <- c(cfg$run_ids, abc_genus = "family_genus_community_v1_20261027",
    abc_within = "within_family_LDM_v1_20261027", abc_external = "external_health_ranking_v2_20261029")
  cfg$output_modules <- c(cfg$output_modules,
    literature_within = "32_Within_family_literature",
    literature_severity = "33_Steatosis_severity",
    literature_metabolic = "34_Metabolic_phenotype_map")
  cfg$run_ids <- c(cfg$run_ids,
    literature_within = "within_family_Chen_style_v1_20261028",
    literature_severity = "grade_fresh_20261101",
    literature_metabolic = "metabolic_map_Asnicar_style_v1_20261028")
  cfg$output_modules <- c(cfg$output_modules,
    generation_phylogeny = "35_Generation_reference_phylogeny")
  cfg$run_ids <- c(cfg$run_ids,
    generation_phylogeny = "reference_phylogeny_v1_20261028")
  cfg$output_modules <- c(cfg$output_modules,
    generation_genus = "36_Genus_generation_taxonomy")
  cfg$run_ids <- c(cfg$run_ids,
    generation_genus = "genus_generation_taxonomy_v1_20261028")
  cfg$output_modules <- c(cfg$output_modules,
    family_specific_genera = "37_Family_specific_genera")
  cfg$run_ids <- c(cfg$run_ids,
    family_specific_genera = "family_genus_Fig1g_source_style_v1_20261029")
  cfg$output_modules <- c(cfg$output_modules,
    steatosis_genus_generation = "38_Steatosis_selected_genus_generation")
  cfg$run_ids <- c(cfg$run_ids,
    steatosis_genus_generation = "steatosis_selected_genera_v1_20261031")
  cfg$output_modules <- c(cfg$output_modules,
    genus_main_text = "39_Genus_main_text_figures")
  cfg$run_ids <- c(cfg$run_ids,
    genus_main_text = "named_genus_reference_style_v4_20261031")
  cfg$output_modules <- c(cfg$output_modules,
    steatosis_sgb_tree = "40_Steatosis_SGB_circular_tree")
  cfg$run_ids <- c(cfg$run_ids,
    steatosis_sgb_tree = "large_legend_SGB_reference_style_v3_20261101")
  cfg$output_modules <- c(cfg$output_modules,
    generation_heterogeneity = "41_Generation_heterogeneity",
    functional_pathways = "42_Functional_pathways",
    strain_sharing = "43_Strain_sharing")
  cfg$run_ids <- c(cfg$run_ids,
    generation_heterogeneity = "generation_interaction_v1_20261101",
    functional_pathways = "humann_pathways_v1_20261101",
    strain_sharing = "strain_sharing_v1_20261101")
  cfg$output_modules <- c(cfg$output_modules,pathway_redesign="44_Pathway_function_sources")
  cfg$run_ids <- c(cfg$run_ids,pathway_redesign="pathway_function_sources_v1_20261102")
  cfg$output_modules <- c(cfg$output_modules, family_discovery = "45_Family_de_novo_discovery",
    family_strain = "46_Family_strain_typing", family_strain_disease = "47_Strain_sharing_disease",
    family_function = "48_Family_functional_explanation")
  cfg$run_ids <- c(cfg$run_ids, family_discovery = "family_discovery_20261104_v1",
    family_strain = "strain_input_gate_20261104_v1", family_strain_disease = "sharing_disease_20261104_v1",
    family_function = "family_function_20261104_v1")
  cfg$output_modules <- c(cfg$output_modules, family_candidates = "50_Family_candidate_characterization")
  cfg$run_ids <- c(cfg$run_ids, family_candidates = "family_candidates_20261104_v1")
  cfg$output_modules <- c(cfg$output_modules, fig6_family_LMM = "51_Fig6_family_mixed_models")
  cfg$run_ids <- c(cfg$run_ids, fig6_family_LMM = "fig6_family_LMM_20261105_v1")
  cfg
}

fl_requested_fields <- function() c(
  "tutknom18", "sukupolvi", "tutkno80perhe", "sp", "syntymapvm",
  "MAKSAPVM18", "TUTKPVM18", "TARKKAIKA18", "ika_status18", "ika18",
  "RASVAMAKSA18", "PITUUS18", "PAINO18", "BMI18", "WAIST18",
  "NOWOB_CAT18", "syst18", "DKVAI18", "GLUK18", "INSU18", "HOMA18",
  "HbA1c18", "HbA1c_alle36kg_18", "TOTKOL18", "TRIGLY18", "HDLKOL18",
  "LDLKOL18", "CRP18", "ALAT18", "ASAT18", "S_GT18", "GT18",
  "smoke18", "alkodpd18", "dietscore_g18", "DIALAA18", "FHLAAK18",
  "VPLAAK18", "DIAB218", "DIAB218_VUOSI", "ANTBV18", "ULOSTE_PVM18"
)
