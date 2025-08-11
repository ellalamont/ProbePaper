# Compare TPMs of samples sequenced on the Sept and Nov runs (Scatterplot)
# E. Lamont
# 7/31/25


source("Import_data.R") # ProbeTest4_tpm, ProbeTest3_tpm


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

my_ProbeTest4_tpm <- ProbeTest4_tpm %>% rename("Gene" = "X")
names(my_ProbeTest4_tpm) <- gsub(x = names(my_ProbeTest4_tpm), pattern = "_S.*", replacement = "")

my_ProbeTest3_tpm <- ProbeTest3_tpm %>% rename("Gene" = "X")
names(my_ProbeTest3_tpm) <- gsub(x = names(my_ProbeTest3_tpm), pattern = "_S.*", replacement = "")

###########################################################
################### SPIKED THP1 COMPARISON ################

THP1_ProbeTest4_tpm <- my_ProbeTest4_tpm %>% select(THP1_1e6_1, THP1_1e6_2, THP1_1e6_3, THP1_1e6_4, THP1_1e6_5, Gene)
THP1_ProbeTest3_tpm <- my_ProbeTest3_tpm %>% select(THP1_1e6_1_Probe_3D_100, THP1_1e6_2_Probe_3D_50, THP1_1e6_3_Probe_3D_25, THP1_1e6_4_Probe_3D_10, THP1_1e6_5_Probe_1, Gene)

THP1_Combined <- inner_join(THP1_ProbeTest4_tpm, THP1_ProbeTest3_tpm, by = "Gene")

THP1_Combined_Log10 <- THP1_Combined %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Add average columns
THP1_Combined_Log10 <- THP1_Combined_Log10 %>% mutate(
  ProbeTest3_Averages = rowMeans(select(., c(THP1_1e6_1_Probe_3D_100, THP1_1e6_2_Probe_3D_50, THP1_1e6_3_Probe_3D_25, THP1_1e6_4_Probe_3D_10, THP1_1e6_5_Probe_1)), na.rm = TRUE),
  ProbeTest4_Averages = rowMeans(select(., c(THP1_1e6_1, THP1_1e6_2, THP1_1e6_3, THP1_1e6_4, THP1_1e6_5)), na.rm = TRUE),
)

# Compare the averages! 
Sample1 <- "ProbeTest4_Averages" # Captured
Sample2 <- "ProbeTest3_Averages" # Not Captured
ScatterCorr <- THP1_Combined_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 0.5, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(# title = paste0("THP1 ProbeTest 3 vs 4: Not scaled Samples AVERAGED: ", Sample1, " vs ", Sample2),
       # subtitle = "Pearson correlation; 5 samples: THP1 1e6 Ra spiked ",
       x = paste0("Log10(TPM+1)\nmixed samples probe A"), 
       y = paste0("Log10(TPM+1)\nmixed samples probe B")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ScatterCorr_2F <- ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = "Averages_THP1.Spiked_CompareAcrossRuns.pdf",
       path = "Figures/Comparing_across_Runs",
       width = 7, height = 5, units = "in")


###########################################################
#################### S_250754 COMPARISON ##################

# Make a new merged DF
W0_ProbeTest4_tpm <- my_ProbeTest4_tpm %>%
  select(S_575533_MtbrRNA, S_687338_MtbrRNA, S_503557, S_250754, THP1_1e6_3, Gene)
W0_ProbeTest3_tpm <- my_ProbeTest3_tpm %>% 
  select(S_575533_Probe_3A, S_687338_Probe_4A_100, S_503557_Probe_3D_10, S_250754_Probe_4A_50, THP1_1e6_3_Probe_3D_25, Gene)

# Merge the dataframes
W0_Combined <- merge(W0_ProbeTest4_tpm, W0_ProbeTest3_tpm, by = "Gene", all = T)

# Log10 transform the data
W0_Combined_Log10 <- W0_Combined %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values


# S_250754 (Nov) and S_250754_Probe_4A_50 (Sept)
Sample1 <- "S_250754"
Sample2 <- "S_250754_Probe_4A_50"
ScatterCorr <- W0_Combined_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") +
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Not Scaled! Log10 transformed, Pearson correlation",
       x = paste0(Sample1, " Log10(TPM+1)"), 
       y = paste0(Sample2, " Log10(TPM+1)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = "Sputum_W0.250754_CompareAcrossRuns.pdf",
       path = "Figures/Comparing_across_Runs",
       width = 7, height = 5, units = "in")










