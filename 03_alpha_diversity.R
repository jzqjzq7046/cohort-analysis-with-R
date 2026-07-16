#==============================================================================
# 03_alpha_diversity.R
#
# Purpose:
#   Alpha diversity analysis for 2018 cross-sectional microbiome data
#
# Output:
#   result/Alpha_diversity.csv
#   result/Alpha_statistics.csv
#   figure/Alpha_diversity.pdf
#==============================================================================


rm(list = ls())

source("00_config.R")


library(vegan)
library(ggplot2)
library(dplyr)


#------------------------------------------------------------------------------
# 1. Load data
#------------------------------------------------------------------------------

metadata <- readRDS(
  file.path(data_dir, "metadata.rds")
)


abundance <- readRDS(
  file.path(data_dir, "abundance_matrix.rds")
)



#------------------------------------------------------------------------------
# 2. Data quality check
#------------------------------------------------------------------------------

# NA check

if(sum(is.na(abundance)) > 0){
  
  abundance[is.na(abundance)] <- 0
  
  cat("NA values replaced by 0\n")
  
}


# Remove samples with zero sequencing depth

depth <- rowSums(abundance)


if(any(depth == 0)){
  
  abundance <- abundance[depth > 0, ]
  
  metadata <- metadata %>%
    filter(
      TUTKNOM18 %in% rownames(abundance)
    )
  
}



cat("Samples used:", nrow(abundance), "\n")
cat("Taxa number :", ncol(abundance), "\n")



#------------------------------------------------------------------------------
# 3. Calculate alpha diversity
#------------------------------------------------------------------------------


Observed <- specnumber(abundance)



Chao1 <- estimateR(abundance)["S.chao1", ]



Shannon <- diversity(
  abundance,
  index = "shannon"
)



Simpson <- diversity(
  abundance,
  index = "simpson"
)



InverseSimpson <- diversity(
  abundance,
  index = "invsimpson"
)



alpha_df <- data.frame(
  
  TUTKNOM18 = rownames(abundance),
  
  Observed = as.numeric(Observed),
  
  Chao1 = as.numeric(Chao1),
  
  Shannon = Shannon,
  
  Simpson = Simpson,
  
  InverseSimpson = InverseSimpson
  
)



# Merge phenotype

alpha_df <-
  left_join(
    alpha_df,
    metadata,
    by = "TUTKNOM18"
  )



#------------------------------------------------------------------------------
# 4. Statistical test
#------------------------------------------------------------------------------


alpha_indexes <- c(
  "Observed",
  "Chao1",
  "Shannon",
  "Simpson",
  "InverseSimpson"
)



alpha_statistics <- data.frame()



for(i in alpha_indexes){
  
  test <- wilcox.test(
    alpha_df[[i]] ~ alpha_df[[outcome]]
  )
  
  
  tmp <- data.frame(
    
    Index = i,
    
    Pvalue = test$p.value
    
  )
  
  
  alpha_statistics <-
    rbind(
      alpha_statistics,
      tmp
    )
  
}



alpha_statistics$FDR <-
  p.adjust(
    alpha_statistics$Pvalue,
    method = "fdr"
  )



#------------------------------------------------------------------------------
# 5. Save results
#------------------------------------------------------------------------------


write.csv(
  
  alpha_df,
  
  file.path(
    result_dir,
    "Alpha_diversity.csv"
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  alpha_statistics,
  
  file.path(
    result_dir,
    "Alpha_statistics.csv"
  ),
  
  row.names = FALSE
  
)



#------------------------------------------------------------------------------
# 6. Visualization
#------------------------------------------------------------------------------


plot_alpha <- function(index){
  
  
  ggplot(
    
    alpha_df,
    
    aes(
      
      x = factor(.data[[outcome]]),
      
      y = .data[[index]],
      
      fill = factor(.data[[outcome]])
      
    )
    
  ) +
    
    
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    
    
    geom_jitter(
      width = 0.15,
      size = 1,
      alpha = 0.4
    ) +
    
    
    theme_bw() +
    
    
    labs(
      
      x = "Hyperlipidemia",
      
      y = index
      
    )
  
}



pdf(
  
  file.path(
    figure_dir,
    "Alpha_diversity.pdf"
  ),
  
  width = 8,
  
  height = 10
  
)



print(plot_alpha("Observed"))

print(plot_alpha("Chao1"))

print(plot_alpha("Shannon"))

print(plot_alpha("Simpson"))

print(plot_alpha("InverseSimpson"))



dev.off()



cat("\n=====================================\n")
cat("Alpha diversity analysis completed!\n")
cat("=====================================\n")