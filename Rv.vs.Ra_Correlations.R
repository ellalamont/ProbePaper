# Correlation (TPM) Graphs comparing my Ra broth with the other Rvs
# E Lamont
# 7/31/25

source("Import_data.R") # 
source("TPM_from_RawReads.R") # To get the calculate TPM function
source("BatchCorrection_PCA.R") # To get the batch corrected stuff

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(2, 2, 2, 2)#
        # plot.margin = margin(10, 10, 10, 20),
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

###########################################################
####################### PROCESS DATA ######################


my_tpm <- my_tpm_bc %>% rownames_to_column("Gene")
names(my_tpm) <- gsub(x = names(my_tpm), pattern = "_S.*", replacement = "")
my_tpm <- my_tpm %>% select(Gene, H37Ra_Broth_4, H37Ra_Broth_5, H37Ra_Broth_6, Rv_pH_7_R1, Rv_pH_7_R2)

### Log10 transform ###
my_tpm_Log10 <- my_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

### Add averages columns ###
my_tpm_Log10 <- my_tpm_Log10 %>% 
  mutate(
    AVERAGE_Ra = rowMeans(select(., c(H37Ra_Broth_4, H37Ra_Broth_5, H37Ra_Broth_6)), na.rm = TRUE),
    AVERAGE_Rv = rowMeans(select(., c(Rv_pH_7_R1, Rv_pH_7_R2)), na.rm = TRUE),
  )


###################################################################
############################ Ra VS Rv #############################

Sample1 <- "AVERAGE_Ra" # Broth Not Captured
Sample2 <- "AVERAGE_Rv" # THP1 spiked Captured
ScatterCorr <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 0.5, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
    subtitle = "Pearson correlation",
    x = paste0("Log10(TPM+1) Ra"),
    y = paste0("Log10(TPM+1) Rv pH7"), ) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = "Ra.vs.LanceRvpH7.pdf",
       path = "Figures_preNonCodingRemoval/Ra.vs.Rv",
       width = 7, height = 5, units = "in")




