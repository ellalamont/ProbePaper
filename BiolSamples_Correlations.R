# Correlations between biological samples (ggcorrplot)
# E. Lamont
# 8/4/25


source("Import_data.R") # for GoodBiolSamples_tpm and GoodBiolSamples_w_Rv_tpm


# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=16, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=16), 
        plot.subtitle = element_text(size=12), 
        # plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )


# Log10 transform the data
my_tpm_Log10 <- GoodBiolSamples_tpm %>% 
  # mutate(Gene = rownames(my_tpm)) %>%
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) %>% # Log transform the values
  column_to_rownames("X")
colnames(my_tpm_Log10) <- gsub(pattern = "_S.*", replacement = "", x = colnames(my_tpm_Log10))

# arrange columns alphabatically
my_tpm_Log10 <- my_tpm_Log10[ , sort(colnames(my_tpm_Log10))]



###########################################################
################# PEARSON LOG10 GGCORRPLOT ################

# Make the correlation
corr <- cor(my_tpm_Log10, method = "pearson")

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(my_tpm_Log10)
# head(p.mat[, 1:4])

min(corr) # 0.1863461
my_min <- round(min(corr), 1)

# Plot pearson
ggcorrplot_PearsonLog10 <- corr %>% 
  ggcorrplot(hc.order = F, 
             method = "square", 
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  # scale_fill_gradient2(limit = c(my_min,1), low = "blue", high =  "red", mid = "white", midpoint = (((1-my_min)/2)+my_min)) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Pearson Correlation Log10(TPM+1)", 
       subtitle = NULL, 
       fill = "Correlation")
ggcorrplot_PearsonLog10

# ggsave(ggcorrplot_PearsonLog10,
#        file = "ggcorrplot_PearsonLog10_v3.pdf",
#        path = "ggcorrplot_Figures",
#        width = 7, height = 6, units = "in")


###########################################################
############### BIOL SAMPLES AVERAGES WITH Rv #############
# 9/11/25
# GoodBiolSamples_w_Rv_tpm

# Log10 transform the data
my_tpm_Log10 <- GoodBiolSamples_w_Rv_tpm %>% 
  # mutate(Gene = rownames(my_tpm)) %>%
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) %>% # Log transform the values
  column_to_rownames("X")

# Get the averages of the samples I want (Sputum, Marm, Rabbit, caseum, Lance's pH7 Rv)
my_data <- my_tpm_Log10 %>%
  mutate(
    AVERAGE_Sputum = rowMeans(dplyr::select(., c(W0_11011_S16, W0_12007_S2, W0_12008_S3, W0_12010_S5, W0_12032_S28, W0_12083_S39, W0_13045_S34, W0_15081_S49, W0_15089_S51)), na.rm = TRUE),
    AVERAGE_CaseumMimic = rowMeans(dplyr::select(., c(HN878_mimic_D14_R1_S63, HN878_mimic_D14_R2_S64, HN878_mimic_D28_R1_S65, HN878_mimic_D7_R2_S62)), na.rm = TRUE),
    AVERAGE_Rabbit = rowMeans(dplyr::select(., c(LLL_Cav_L2_18_S58, LLL_Cav_L2_21_S59, LU_Cav_L1_12_S57)), na.rm = TRUE),
    AVERAGE_Marmoset = rowMeans(dplyr::select(., c(BQ12_10_Probe_3A_S29, BQ12_8_Probe_4A_50_S28)), na.rm = TRUE),
    AVERAGE_pH7Rv = rowMeans(dplyr::select(., c(Rv_pH_7_R1, Rv_pH_7_R2)), na.rm = TRUE),
    AVERAGE_Ra = rowMeans(dplyr::select(., c(H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9)), na.rm = TRUE),
         ) %>%
  dplyr::select(AVERAGE_Sputum, AVERAGE_CaseumMimic, AVERAGE_Rabbit, AVERAGE_Marmoset, AVERAGE_pH7Rv, AVERAGE_Ra)

# Make the correlation
corr <- cor(my_data, method = "pearson")

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(my_tpm_Log10)
# head(p.mat[, 1:4])

# Plot pearson
ggcorrplot_PearsonLog10 <- corr %>% 
  ggcorrplot(hc.order = F, 
             method = "square", 
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Pearson Correlation Log10(TPM+1)", 
       subtitle = NULL, 
       fill = "Correlation")
ggcorrplot_PearsonLog10
ggsave(ggcorrplot_PearsonLog10,
       file = "ggcorrplot_PearsonLog10_BiolSamples_Averages.pdf",
       path = "Figures_preNonCodingRemoval/ggcorrplot",
       width = 7, height = 6, units = "in")





