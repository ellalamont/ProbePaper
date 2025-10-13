# Figure 2

source("FINAL_ImportData.R")


###########################################################
###################### FIGURE 2A ##########################

my_plot_themes <- theme_void() +
  theme(legend.position = "none",legend.text=element_text(size=14),
        strip.text = element_text(size = 12, face = "bold"), # For the facet
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10, hjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        plot.subtitle = element_text(size=10), 
        plot.margin = margin(2, 2, 2, 2))

### AVERAGE THE REPLICATES ###
Averages_CapturedVsNot_pipeSummary <- CapturedVsNot_pipeSummary %>%
  group_by(Probe) %>%
  summarize(count = n(),
            mean_P_Genomic = mean(P_Genomic_Rv),
            mean_P_RiboClear = mean(P_RiboClear),
            mean_P_NoHit = mean(P_NoHit),
            mean_P_ncRNA = mean(P_NonCodingMtbRNA)) %>%
  pivot_longer(cols = c("mean_P_Genomic", "mean_P_RiboClear", "mean_P_NoHit", "mean_P_ncRNA"), names_to = "Percent_Type", values_to = "Percent") %>%
  mutate(Percent = round(Percent, 1)) 

replacement_values <- c(mean_P_Genomic = "mRNA", mean_P_RiboClear = "rRNA", mean_P_NoHit = "Host RNA", mean_P_ncRNA = "ncRNA")
Averages_CapturedVsNot_pipeSummary <- Averages_CapturedVsNot_pipeSummary %>% 
  mutate(Percent_Type = replacement_values[Percent_Type])

### PIE CHART OF AVERAGES ###
Fig2A <- Averages_CapturedVsNot_pipeSummary %>% 
  arrange(Probe, desc(Percent_Type)) %>%  # Need this so the numbers go to the correct slices
  group_by(Probe) %>% # Need this or it won't be a round pie! 
  mutate(cumulative = cumsum(Percent), midpoint = cumulative - Percent / 2) %>%
  ungroup() %>%
  ggplot(aes(x = "", y = Percent, fill = Percent_Type)) +
  geom_bar(width = 1, stat = "identity", color = "black") + 
  coord_polar(theta = "y", start = 0) + 
  facet_wrap(~Probe, labeller = as_labeller(c("Uncaptured" = "Uncaptured average", "Captured" = "Captured average"))) +
  scale_fill_manual(values = c("#708090", "#00CED1", "#B9AED3", "#E0D8B0")) + 
  geom_text_repel(aes(y = midpoint, label = paste(Percent_Type, "\n", scales::percent(Percent / 100))), size = 4, color = "black", box.padding = 0.3, force = 2, force_pull = 2, min.segment.length = 0.2, segment.size = 0.5) + 
  # labs(title = "AVERAGES THP1 cells spiked with 1e6 H37Ra") + 
  my_plot_themes
Fig2A


###########################################################
###################### FIGURE 2B ##########################

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none", legend.text=element_text(size=10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(2, 2, 2, 2))

Fig2B <- CapturedVsNot_pipeSummary %>%
  ggplot(aes(x = Probe, y = N_Genomic_Rv)) + 
  geom_point(size = 2.5, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_line(aes(group = Replicates), color = "black", linewidth = 0.5, linetype = "dashed") + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,17000000), breaks = seq(0, 17000000, 4000000)) + 
  labs(x = NULL, y = "# reads aligning to Mtb") + 
  my_plot_themes
Fig2B

###########################################################
###################### FIGURE 2C ##########################

Fig2C <- CapturedVsNot_pipeSummary %>% 
  ggerrorplot(x = "Probe", y = "Txn_Coverage", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4, add.params = list(size = 0.4)) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  geom_hline(yintercept = 80, linetype = "dashed", alpha = 0.5) + 
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) + 
  labs(x = NULL, y = "% transcriptional coverage") + 
  my_plot_themes
Fig2C

###########################################################
###################### FIGURE 2D ##########################

Fig2D <- LimitofDetect_pipeSummary %>% 
  ggerrorplot(x = "Ra_cells2", y = "N_Genomic_Rv", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4, add.params = list(size = 0.4)) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,17000000), breaks = seq(0, 17000000, 4000000)) + 
  labs(x = "# Mtb cells in spiked sample", y = "# reads aligning to Mtb") + 
  my_plot_themes
Fig2D

###########################################################
###################### FIGURE 2E ##########################

Fig2E <- LimitofDetect_pipeSummary %>% 
  ggerrorplot(x = "Ra_cells2", y = "P_Genomic_Rv", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4, add.params = list(size = 0.4)) + 
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  labs(x = "# Mtb cells in spiked sample", y = "% reads aligning to Mtb") + 
  my_plot_themes
Fig2E


###########################################################
###################### FIGURE 2F ##########################

Fig2F <- LimitofDetect_pipeSummary %>% 
  ggerrorplot(x = "Ra_cells2", y = "Txn_Coverage", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4, add.params = list(size = 0.4)) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  geom_hline(yintercept = 80, linetype = "dashed", alpha = 0.5) + 
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) + 
  labs(x = "# Mtb cells in spiked sample", y = "% transcriptional coverage") + 
  my_plot_themes
Fig2F

###########################################################
###################### FIGURE 2G ##########################

my_tpm <- All_tpm # %>% rename("Gene" = "X")
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

my_tpm_Log10 <- my_tpm_Log10 %>% rownames_to_column("Gene")

Sample1 <- "AVERAGE_THP1Spiked_NotCaptured" 
Sample2 <- "AVERAGE_THP1Spiked" 
Fig2G <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 0.5, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(x = paste0("Log10(TPM+1)\nUncaptured spiked samples"),
    y = paste0("Log10(TPM+1)\nCaptured spiked samples"), ) + 
  stat_cor(method="pearson") + 
  my_plot_themes
Fig2G


###########################################################
###################### FIGURE 2H ##########################

my_tpm <- All_tpm
names(my_tpm) <- gsub(x = names(my_tpm), pattern = "_S.*", replacement = "")
my_tpm <- my_tpm %>% select(ProbeA_THP1_1e6_1, ProbeA_THP1_1e6_2, ProbeA_THP1_1e6_3, ProbeA_THP1_1e6_4, ProbeA_THP1_1e6_5, 
                            ProbeB_THP1_1e6_1, ProbeB_THP1_1e6_2, ProbeB_THP1_1e6_3, ProbeB_THP1_1e6_4, ProbeB_THP1_1e6_5) %>%
  rownames_to_column("Gene")

my_tpm_Log10 <- my_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Add average columns
my_tpm_Log10 <- my_tpm_Log10 %>% mutate(
  ProbeA_Averages = rowMeans(select(., c(ProbeA_THP1_1e6_1, ProbeA_THP1_1e6_2, ProbeA_THP1_1e6_3, ProbeA_THP1_1e6_4, ProbeA_THP1_1e6_5)), na.rm = TRUE),
  ProbeB_Averages = rowMeans(select(., c(ProbeB_THP1_1e6_1, ProbeB_THP1_1e6_2, ProbeB_THP1_1e6_3, ProbeB_THP1_1e6_4, ProbeB_THP1_1e6_5)), na.rm = TRUE),
)

Sample1 <- "ProbeB_Averages" # Captured
Sample2 <- "ProbeA_Averages" # Not Captured
Fig2H <- my_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 0.5, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(x = paste0("Log10(TPM+1)\nSpiked samples probe B"), 
    y = paste0("Log10(TPM+1)\nSpiked samples probe A")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
Fig2H


###########################################################
##################### COMBINE FIGURES #####################

top <- plot_grid(Fig2A, Fig2B, Fig2C,
                 Fig2D, Fig2E, Fig2F,
                 ncol = 3, align = "hv", axis = "tblr")
bottom <- plot_grid(NULL,  # left spacer
                    plot_grid(Fig2G, Fig2H, ncol = 2, rel_widths = c(1,1)),
                    NULL,  # right spacer
                    ncol = 3, rel_widths = c(0.2, 1, 0.2))  # adjust spacers vs. center width
combined <- plot_grid(
  top, bottom,
  ncol = 1, rel_heights = c(2, 1))  # adjust ratio as needed
combined

ggsave(combined,
       file = paste0("FINAL_Figure2.pdf"),
       path = "Figures/CombinedFigures",
       width = 15, height = 12, units = "in")


###########################################################
######## NUMBER OF GENES WITH 0 READS (UNCAPTURED) ########

# THP1_1e6_1b_S29, THP1_1e6_2a_S30, THP1_1e6_3b_S33

Uncaptured_RawReads <- All_RawReads %>% select(X, THP1_1e6_1b_S29, THP1_1e6_2a_S30, THP1_1e6_3b_S33)

NumNoReads <- colSums(Uncaptured_RawReads == 0)
# X THP1_1e6_1b_S29 THP1_1e6_2a_S30 THP1_1e6_3b_S33 
# 0             150             132             687 
(150 + 132 + 687) / 3
mean(c(150, 132, 687))
sd(c(150, 132, 687))
