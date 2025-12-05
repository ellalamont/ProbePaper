# Correlations between the individual sputum samples (and Ra)
# 12/4/25
# In response to reviewer comments

source("Import_data.R") # for All_tpm_f

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 0, size=16, vjust=0, hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=16), 
        plot.subtitle = element_text(size=12))


###########################################################
##################### ORGANIZE THE DATA ###################

# Just want my 9 sputum samples and the Rv samples (which will have to be averaged)
my_data <- All_tpm_f %>% select(all_of(Sputum_GoodSampleList), contains("H37Ra"))

# Log10 transform the data
my_tpm_Log10 <- my_data %>% 
  # mutate(Gene = rownames(my_tpm)) %>%
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # %>% # Log transform the values
  # column_to_rownames("X")
colnames(my_tpm_Log10) <- gsub(pattern = "_S.*", replacement = "", x = colnames(my_tpm_Log10))

# arrange columns alphabatically
my_tpm_Log10 <- my_tpm_Log10[ , sort(colnames(my_tpm_Log10))]


# Get the averages of the Ra samples
my_data2 <- my_tpm_Log10 %>%
  mutate(AVERAGE_Ra = rowMeans(dplyr::select(., c(H37Ra_Broth_4, H37Ra_Broth_5, H37Ra_Broth_6)), na.rm = TRUE)) %>%
  dplyr::select(-c(H37Ra_Broth_4, H37Ra_Broth_5, H37Ra_Broth_6))


###########################################################
############### RARITY ALL GRAPHS TOGETHER ################

# https://borisleroy.com/en/2013/06/09/correlation-plots-in-r/
# install.packages("Rarity")
library(Rarity)


# Pearson
pdf("Figures/IndividualSputum_Correlations/rarity_PearsonLog10_v1.pdf", width = 10, height = 10)
corPlot(my_data2, method = "pearson",
        title = "Pearson Correlation Log10(TPM+1)")
dev.off()

# These really don't look that good....


###########################################################
############### NORMAL CORRELATION MATRIX #################

# Make the correlation
corr <- cor(my_data2, method = "pearson")

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(my_data2)
# head(p.mat[, 1:4])

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(my_data2, method = "pearson")
p.mat.adj <- p.mat %>%
  as.vector() %>% 
  p.adjust(method = "fdr") %>%
  matrix(nrow = nrow(p.mat), ncol = ncol(p.mat))
rownames(p.mat.adj) <- rownames(p.mat)
colnames(p.mat.adj) <- colnames(p.mat)

# Plot pearson
ggcorrplot_PearsonLog10 <- corr %>% 
  ggcorrplot(hc.order = F, 
             p.mat = p.mat.adj,
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
       path = "Figures/IndividualSputum_Correlations",
       width = 7, height = 6, units = "in")





