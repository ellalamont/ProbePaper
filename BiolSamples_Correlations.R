# Correlations between biological samples (ggcorrplot)
# E. Lamont
# 8/4/25


source("Import_data.R") # for GoodBiolSamples_tpm


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







