# ============================================================================
# 01_prepare_data.R
# Purpose:
#   Prepare 2018 cross-sectional microbiome dataset
# ============================================================================

rm(list = ls())

source("00_config.R")

library(mia)
library(haven)
library(dplyr)
library(tibble)

cat("============================================================\n")
cat("Preparing 2018 cross-sectional dataset ...\n")
cat("============================================================\n")

#------------------------------------------------------------------------------
# 1. Load data
#------------------------------------------------------------------------------

tse <- readRDS(file.path(data_dir, "tse.rds"))

df_host <- read_dta(
  file.path(data_dir, "df_final_combined_2011_2018_MERGED.dta")
)

df_host <- as.data.frame(df_host)
colnames(df_host) <- toupper(colnames(df_host))

#------------------------------------------------------------------------------
# 2. Extract 2018 phenotype
#------------------------------------------------------------------------------

df_clinical_18 <-
  df_host %>%
  filter(
    !is.na(TUTKNOM18),
    TUTKNOM18 != ""
  ) %>%
  select(
    TUTKNOM18,
    IKA18,
    SP,
    
    BMI18,
    WAIST18,
    
    GLUCOSE18,
    TRIGLY18,
    LDLC18,
    HDLC18,
    
    CVHCHOL18,
    DIALAA18
  ) %>%
  mutate(
    
    DIABETES_18 =
      ifelse(
        GLUCOSE18 >= 7 |
          DIALAA18 == 1,
        1,
        0
      ),
    
    HYPERLIPIDEMIA_18 =
      case_when(
        
        CVHCHOL18 == 0 ~ 1,
        
        TRIGLY18 >= 2.3 ~ 1,
        
        LDLC18 >= 4.1 ~ 1,
        
        HDLC18 < 1.0 ~ 1,
        
        TRUE ~ 0
      )
    
  ) %>%
  distinct(TUTKNOM18, .keep_all = TRUE)

df_clinical_18$TUTKNOM18 <-
  toupper(
    trimws(
      as.character(df_clinical_18$TUTKNOM18)
    )
  )

#------------------------------------------------------------------------------
# 3. Extract abundance matrix
#------------------------------------------------------------------------------

abundance <-
  as.data.frame(
    t(
      assay(tse, "counts")
    )
  ) %>%
  rownames_to_column("TUTKNOM18")

abundance$TUTKNOM18 <-
  toupper(
    trimws(
      as.character(abundance$TUTKNOM18)
    )
  )

#------------------------------------------------------------------------------
# 4. Match microbiome with phenotype
#------------------------------------------------------------------------------

df_cross_sectional_18 <-
  inner_join(
    df_clinical_18,
    abundance,
    by = "TUTKNOM18"
  )

#------------------------------------------------------------------------------
# 5. Metadata
#------------------------------------------------------------------------------

metadata <-
  df_cross_sectional_18 %>%
  select(
    TUTKNOM18,
    
    IKA18,
    SP,
    
    BMI18,
    WAIST18,
    
    GLUCOSE18,
    TRIGLY18,
    LDLC18,
    HDLC18,
    
    CVHCHOL18,
    DIALAA18,
    
    DIABETES_18,
    HYPERLIPIDEMIA_18
  )

#------------------------------------------------------------------------------
# 6. Taxa abundance matrix
#------------------------------------------------------------------------------

taxa_names <-
  setdiff(
    colnames(df_cross_sectional_18),
    colnames(metadata)
  )

abundance_matrix <-
  df_cross_sectional_18 %>%
  select(all_of(taxa_names))

abundance_matrix <-
  as.data.frame(abundance_matrix)

rownames(abundance_matrix) <-
  df_cross_sectional_18$TUTKNOM18

#------------------------------------------------------------------------------
# 7. Synchronize TreeSummarizedExperiment
#------------------------------------------------------------------------------

sample_keep <-
  match(
    metadata$TUTKNOM18,
    colnames(assay(tse, "counts"))
  )

tse_cross <- tse[, sample_keep]

colData(tse_cross) <-
  S4Vectors::DataFrame(metadata)

#------------------------------------------------------------------------------
# 8. Summary
#------------------------------------------------------------------------------

cat("\n")
cat("Matched samples :", nrow(metadata), "\n")
cat("Number of taxa  :", ncol(abundance_matrix), "\n\n")

cat("Hyperlipidemia\n")
print(table(metadata$HYPERLIPIDEMIA_18))

cat("\nDiabetes\n")
print(table(metadata$DIABETES_18))

#------------------------------------------------------------------------------
# 9. Save objects
#------------------------------------------------------------------------------

saveRDS(
  df_cross_sectional_18,
  file.path(data_dir, "df_cross_sectional_18.rds")
)

saveRDS(
  metadata,
  file.path(data_dir, "metadata.rds")
)

saveRDS(
  abundance_matrix,
  file.path(data_dir, "abundance_matrix.rds")
)

saveRDS(
  taxa_names,
  file.path(data_dir, "taxa_names.rds")
)

saveRDS(
  tse_cross,
  file.path(data_dir, "tse_cross_sectional.rds")
)

cat("\nObjects saved:\n")
cat("  df_cross_sectional_18.rds\n")
cat("  metadata.rds\n")
cat("  abundance_matrix.rds\n")
cat("  taxa_names.rds\n")
cat("  tse_cross_sectional.rds\n")
cat("============================================================\n")
cat("Finished!\n")


#------------------------------------------------------------------------------
# Save taxonomy information
#------------------------------------------------------------------------------

taxonomy <- as.data.frame(
  rowData(tse_cross)
)

taxonomy$Feature_ID <- rownames(taxonomy)


write.csv(
  taxonomy,
  file.path(
    result_dir,
    "taxonomy_annotation.csv"
  ),
  row.names = FALSE
)