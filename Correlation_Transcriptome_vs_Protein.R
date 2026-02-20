# Title: Correlation Analysis between Transcriptomic Signatures and Protein Abundance
# Author: Deepika Trehan, Ph.D.
# Description: This script performs cross-platform validation between mRNA (Agilent Microarray) 
# and Protein (IHC Scoring/Western Blot Densitometry) for hub genes identified in Bladder Cancer.

# 1. Load Libraries
library(ggplot2)
library(ggpubr)
library(dplyr)

# 2. Mock Data Structure (Based on AJCR 2023 Publication)
# mRNA_val: Normalized intensity from Agilent/TCGA
# Protein_val: Densitometry score or H-score from IHC
data_correlation <- data.frame(
  Sample_ID = paste0("BC_", 1:10),
  mRNA_VEGFA = c(4.2, 5.1, 3.8, 6.2, 5.5, 4.0, 7.1, 6.8, 5.9, 4.5),
  Protein_VEGFA = c(1.2, 1.5, 1.1, 1.9, 1.6, 1.2, 2.1, 2.0, 1.8, 1.3)
)

# 3. Statistical Correlation (Spearman/Pearson)
# Used to determine the strength of the mRNA-Protein link for hub genes
corr_test <- cor.test(data_correlation$mRNA_VEGFA, 
                      data_correlation$Protein_VEGFA, 
                      method = "spearman")

print(corr_test)

# 4. Visualization: Scatter Plot with Regression Line
# This plot visually demonstrates the 'Wet-to-Dry' validation loop
ggplot(data_correlation, aes(x = mRNA_VEGFA, y = Protein_VEGFA)) +
  geom_point(color = "#2c3e50", size = 3) +
  geom_smooth(method = "lm", color = "#e74c3c", fill = "#ecf0f1") +
  stat_cor(method = "spearman", label.x = 4, label.y = 2) +
  labs(title = "mRNA vs Protein Correlation: VEGFA Hub Gene",
       subtitle = "Validation of Agilent Microarray via Western Blot/IHC",
       x = "mRNA Expression (log2 Intensity)",
       y = "Protein Abundance (Relative Densitometry)") +
  theme_minimal()

# 5. Result Export
# ggsave("VEGFA_Correlation_Plot.png") 
