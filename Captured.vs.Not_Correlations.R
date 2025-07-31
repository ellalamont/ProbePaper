# Correlation (TPM) Graphs comparing the THP1 spiked samples that were captured and not captured 
# E Lamont
# 7/31/25

source("Import_data.R") # ProbeTest5_tpm_CapturedVsNot_wBroth

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
        # plot.margin = margin(10, 10, 10, 20),
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

###########################################################
####################### PROCESS DATA ######################

my_tpm <- ProbeTest5_tpm_CapturedVsNot_wBroth %>% rename("Gene" = "X")
names(my_tpm) <- gsub(x = names(my_tpm), pattern = "_S.*", replacement = "")

### Log10 transform ###
my_tpm_Log10 <- my_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

### Add averages columns ###
my_tpm_Log10 <- my_tpm_Log10 %>% 
  mutate(
    AVERAGE_THP1Spiked = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    AVERAGE_THP1Spiked_NotCaptured = rowMeans(select(., c(THP1_1e6_1b, THP1_1e6_2a, THP1_1e6_3b)), na.rm = TRUE),
    AVERAGE_BrothNotCaptured = rowMeans(select(., c(H37Ra_Broth_4, H37Ra_Broth_5, H37Ra_Broth_6)), na.rm = TRUE),
  )


###################################################################
####################### THP1Spiked VS BROTH #######################

Sample1 <- "AVERAGE_BrothNotCaptured" # Broth Not Captured
Sample2 <- "AVERAGE_THP1Spiked" # THP1 spiked Captured
ScatterCorr <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; 1e6 Ra THP1 spiked captured VS Broth Not captured (Not scaled)",
       x = paste0("Log10(TPM+1) Not captured broth samples"),
       y = paste0("Log10(TPM+1) Captured mixed samples"), ) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = "THP1SpikedCaptured.vs.Broth.pdf",
       path = "Figures/Captured.vs.Not_Correlations",
       width = 7, height = 5, units = "in")


###################################################################
#################### THP1Spiked CAPTURED VS NOT ###################

Sample1 <- "AVERAGE_THP1Spiked_NotCaptured" # Broth Not Captured
Sample2 <- "AVERAGE_THP1Spiked" # THP1 spiked Captured
ScatterCorr <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; 1e6 Ra THP1 spiked captured VS Not captured spiked samples",
       x = paste0("Log10(TPM+1) Not captured mixed samples"),
       y = paste0("Log10(TPM+1) Captured mixed samples"), ) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = "THP1SpikedCaptured.vs.NotCaptured_fig1.pdf",
       path = "Figures/Captured.vs.Not_Correlations",
       width = 7, height = 5, units = "in")











