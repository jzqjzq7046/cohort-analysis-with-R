#==============================================================================
# 04_beta_diversity.R
#
# Beta diversity analysis
#
# Methods:
#   Bray-Curtis
#   Jaccard
#   Aitchison (CLR)
#   PCoA
#   PERMANOVA
#   betadisper
#
#==============================================================================


rm(list = ls())

source("00_config.R")


library(vegan)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(compositions)



#------------------------------------------------------------------------------
# 1. Load data
#------------------------------------------------------------------------------

metadata <- readRDS(
  file.path(data_dir, "metadata.rds")
)


abundance <- readRDS(
  file.path(data_dir, "abundance_matrix.rds")
)

metadata$HYPERLIPIDEMIA_18 <-
  factor(
    metadata$HYPERLIPIDEMIA_18,
    levels = c(0,1),
    labels = c("Healthy","Hyperlipidemia")
  )

# 保证顺序一致

metadata <- metadata %>%
  filter(
    TUTKNOM18 %in% rownames(abundance)
  )


abundance <- abundance[
  metadata$TUTKNOM18,
]


rownames(metadata) <- metadata$TUTKNOM18



cat(
  "Samples:",
  nrow(abundance),
  "\n"
)



#------------------------------------------------------------------------------
# 2. Remove empty taxa
#------------------------------------------------------------------------------

abundance <- abundance[
  ,
  colSums(abundance) > 0
]


cat(
  "Taxa after filtering:",
  ncol(abundance),
  "\n"
)



#------------------------------------------------------------------------------
# 3. Bray-Curtis distance
#------------------------------------------------------------------------------

bray_dist <- vegdist(
  
  abundance,
  
  method = "bray"
  
)



#------------------------------------------------------------------------------
# 4. Jaccard distance
#------------------------------------------------------------------------------

jaccard_dist <- vegdist(
  
  abundance,
  
  method = "jaccard",
  
  binary = TRUE
  
)



#------------------------------------------------------------------------------
# 5. Aitchison distance
#------------------------------------------------------------------------------

# CLR transformation

abundance_pseudo <-
  abundance + 1


clr_abundance <-
  compositions::clr(
    abundance_pseudo
  )


aitchison_dist <-
  dist(
    clr_abundance,
    method = "euclidean"
  )



#------------------------------------------------------------------------------
# 6. PCoA function
#------------------------------------------------------------------------------


run_pcoa <- function(distance_matrix){
  
  
  pcoa <- cmdscale(
    
    distance_matrix,
    
    k = 2,
    
    eig = TRUE
    
  )
  
  
  pcoa_df <- data.frame(
    
    Sample =
      rownames(abundance),
    
    PC1 =
      pcoa$points[,1],
    
    PC2 =
      pcoa$points[,2]
    
  )
  
  
  pcoa_df <-
    left_join(
      
      pcoa_df,
      
      metadata,
      
      by=c(
        "Sample"="TUTKNOM18"
      )
      
    )
  
  
  return(pcoa_df)
  
}



bray_pcoa <-
  run_pcoa(bray_dist)


jaccard_pcoa <-
  run_pcoa(jaccard_dist)


aitchison_pcoa <-
  run_pcoa(aitchison_dist)



#------------------------------------------------------------------------------
# 7. PERMANOVA
#------------------------------------------------------------------------------


run_permanova <- function(distance){
  
  
  model1 <-
    adonis2(
      
      distance ~
        
        HYPERLIPIDEMIA_18,
      
      data = metadata,
      
      permutations = 999
      
    )
  
  
  
  model2 <-
    adonis2(
      
      distance ~
        
        HYPERLIPIDEMIA_18 +
        
        IKA18 +
        
        SP +
        
        BMI18,
      
      data = metadata,
      
      permutations = 999
      
    )
  
  
  list(
    
    crude=model1,
    
    adjusted=model2
    
  )
  
}



bray_permanova <-
  run_permanova(bray_dist)


jaccard_permanova <-
  run_permanova(jaccard_dist)


aitchison_permanova <-
  run_permanova(aitchison_dist)



#------------------------------------------------------------------------------
# 8. Dispersion test
#------------------------------------------------------------------------------


bray_disp <-
  betadisper(
    
    bray_dist,
    
    metadata$HYPERLIPIDEMIA_18
    
  )


disp_test <-
  anova(bray_disp)



#------------------------------------------------------------------------------
# 9. Save statistics
#------------------------------------------------------------------------------


capture.output(
  
  bray_permanova,
  
  file=file.path(
    result_dir,
    "PERMANOVA_Bray.txt"
  )
  
)


capture.output(
  
  jaccard_permanova,
  
  file=file.path(
    result_dir,
    "PERMANOVA_Jaccard.txt"
  )
  
)


capture.output(
  
  aitchison_permanova,
  
  file=file.path(
    result_dir,
    "PERMANOVA_Aitchison.txt"
  )
  
)


capture.output(
  
  disp_test,
  
  file=file.path(
    result_dir,
    "betadisper_Bray.txt"
  )
  
)



#------------------------------------------------------------------------------
# 10. PCoA plot
#------------------------------------------------------------------------------


plot_pcoa <- function(df,title){
  
  
  ggplot(
    
    df,
    
    aes(
      
      PC1,
      
      PC2,
      
      color=factor(HYPERLIPIDEMIA_18)
      
    )
    
  )+
    
    
    geom_point(
      
      size=2,
      
      alpha=0.7
      
    )+
    
    
    stat_ellipse(
      
      level=0.95,
      
      linetype=2
      
    )+
    
    
    theme_bw()+
    
    
    labs(
      
      color="Hyperlipidemia",
      
      title=title
      
    )
  
}



pdf(
  
  file.path(
    
    figure_dir,
    
    "Beta_diversity_PCoA.pdf"
    
  ),
  
  width=12,
  
  height=4
  
)



print(
  plot_pcoa(
    bray_pcoa,
    "Bray-Curtis PCoA"
  )
)


print(
  plot_pcoa(
    jaccard_pcoa,
    "Jaccard PCoA"
  )
)


print(
  plot_pcoa(
    aitchison_pcoa,
    "Aitchison PCoA"
  )
)


dev.off()



#------------------------------------------------------------------------------
# 11. Save PCoA coordinates
#------------------------------------------------------------------------------


write.csv(
  
  bray_pcoa,
  
  file.path(
    result_dir,
    "PCoA_Bray_coordinates.csv"
  ),
  
  row.names=FALSE
  
)


write.csv(
  
  jaccard_pcoa,
  
  file.path(
    result_dir,
    "PCoA_Jaccard_coordinates.csv"
  ),
  
  row.names=FALSE
  
)


write.csv(
  
  aitchison_pcoa,
  
  file.path(
    result_dir,
    "PCoA_Aitchison_coordinates.csv"
  ),
  
  row.names=FALSE
  
)



cat("\n====================================\n")
cat("Beta diversity analysis completed!\n")
cat("====================================\n")