# Graphs comparing the THP1 spiked samples that were captured and not captured 
# E Lamont
# 7/31/25

source("Import_data.R") # to get CapturedVsNot_pipeSummary

# Plot basics
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
        plot.margin = margin(2, 2, 2, 2) #
        # plot.margin = margin(10, 10, 10, 20),
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default


###########################################################
########## ALL CELL NUMBER VS AtLeast.10.Reads ############

#ggerrroplot
CapturedvsNot_10Reads_Fig1 <- CapturedVsNot_pipeSummary %>% 
  mutate(Txn_Coverage = (AtLeast.10.Reads/4499)*100) %>%
  ggerrorplot(x = "Probe", y = "Txn_Coverage", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4,  # Size of error bars
              add.params = list(size = 0.4)) +  # Size of mean points
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  # scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, 500)) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  geom_hline(yintercept = 80, linetype = "dashed", alpha = 0.5) + 
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) + 
  # geom_hline(yintercept = 4499*0.8, linetype = "dashed", alpha = 0.5) + 
  # annotate("text", x = 0.6, y = 4499*0.8, label = "80%", 
  #          hjust = 1, vjust = -0.5, color = "black") + 
  # geom_hline(yintercept = 4499*0.5, linetype = "dashed", alpha = 0.5) + 
  # annotate("text", x = 0.6, y = 4499*0.5, label = "50%", 
  #          hjust = 1, vjust = -0.5, color = "black") + 
  labs(# title = "ProbeTest5 THP1 cells spiked with H37Ra", 
       # subtitle = "Mean with standard deviation", 
       # y = "# genes with >=10 reads aligning",
      y = "% transcriptional coverage") + 
  scale_x_discrete(labels = c("None" = "Uncaptured",
                              "JA2" = "Captured")) + 
  my_plot_themes + theme(axis.title.x = element_blank())
CapturedvsNot_10Reads_Fig1
# ggsave(CapturedvsNot_10Reads_Fig1,
#        file = "CapturedvsNot_10Reads_Fig1.pdf",
#        path = "Figures/Captured.vs.Not_ReadCounts",
#        width = 6, height = 4, units = "in")



###########################################################
################# HOW MANY HAVE >1 READ ###################

# 7/31/25 - not run yet, haven't pulled in the TPM

# See how many have at least one read

# counts <- CapturedVsNot_tpm %>%
  # summarise(across(everything(), ~ sum(. > 0)))
# result <- 4499-counts # This is how many genes have 0 reads aligning in each!
# THP1_1e6_1a THP1_1e6_1b THP1_1e6_2a THP1_1e6_2b THP1_1e6_3a THP1_1e6_3b
# 1          14         265         247          15          12         923


###########################################################
###################### N_GENOMIC ##########################
CapturedVsNot_N.Genomic_fig1 <- CapturedVsNot_pipeSummary %>%
  ggplot(aes(x = Probe, y = N_Genomic)) + 
  geom_point(size = 2.5, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_line(aes(group = Replicates), color = "black", linewidth = 0.5, linetype = "dashed") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,19000000), breaks = seq(0, 19000000, 2000000)) +
  scale_y_continuous(limits = c(0,20000000), breaks = seq(0, 20000000, 4000000)) + 
  labs(# title = "THP1 spiked with 1e6 H37Ra",
       subtitle = NULL, 
       x = NULL, 
       y = "# reads aligning to Mtb") + 
  scale_x_discrete(labels = c("None" = "Uncaptured",
                              "JA2" = "Captured")) + 
  my_plot_themes
CapturedVsNot_N.Genomic_fig1
# ggsave(CapturedVsNot_N.Genomic_fig1,
#        file = "CapturedVsNot_N.Genomic_fig1.pdf",
#        path = "Figures/Captured.vs.Not_ReadCounts",
#        width = 6, height = 4, units = "in")


###########################################################
###################### P_GENOMIC ##########################
CapturedVsNot_P.Genomic_fig1 <- CapturedVsNot_pipeSummary %>%
  ggplot(aes(x = Probe, y = P_Genomic)) + 
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_line(aes(group = Replicates), color = "black", size = 0.5, linetype = "dashed") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 3, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) +
  labs(title = "THP1 spiked with 1e6 H37Ra",
       subtitle = NULL, 
       x = NULL, 
       y = "% reads aligning to Mtb genome") + 
  scale_x_discrete(labels = c("None" = "Uncaptured",
                              "JA2" = "Captured")) + 
  my_plot_themes
CapturedVsNot_P.Genomic_fig1
ggsave(CapturedVsNot_P.Genomic_fig1,
       file = "CapturedVsNot_P.Genomic_fig1.pdf",
       path = "Figures/Captured.vs.Not_ReadCounts",
       width = 6, height = 4, units = "in")


###########################################################
###################### PIE CHART ##########################

my_plot_themes <- theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        strip.text = element_text(size = 12, face = "bold"), # For the facet
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10, hjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        plot.subtitle = element_text(size=10), 
        plot.margin = margin(2, 2, 2, 2), #
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_blank(), # Remove facet panel borders
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank())

### AVERAGE THE REPLICATES ###
Averages_CapturedVsNot_pipeSummary <- CapturedVsNot_pipeSummary %>%
  group_by(Probe) %>%
  summarize(count = n(),
            mean_P_Genomic = mean(P_Genomic),
            mean_P_RiboClear = mean(P_RiboClear),
            mean_P_NoHit = mean(P_NoHit)) %>%
  pivot_longer(cols = c("mean_P_Genomic", "mean_P_RiboClear", "mean_P_NoHit"), names_to = "Percent_Type", values_to = "Percent") %>%
  mutate(Percent = round(Percent, 1)) 

replacement_values <- c(mean_P_Genomic = "mRNA", mean_P_RiboClear = "rRNA", mean_P_NoHit = "other RNA")
Averages_CapturedVsNot_pipeSummary <- Averages_CapturedVsNot_pipeSummary %>% 
  mutate(Percent_Type = replacement_values[Percent_Type])

### PIE CHART OF AVERAGES ###
PieChart_Averages_fig1 <- Averages_CapturedVsNot_pipeSummary %>% 
  filter(Probe %in% c("JA2", "None")) %>%
  mutate(Probe = factor(Probe, levels = c("None", "JA2"))) %>% 
  arrange(Probe, desc(Percent_Type)) %>%  # Need this so the numbers go to the correct slices
  group_by(Probe) %>% # Need this or it won't be a round pie! 
  mutate(cumulative = cumsum(Percent), midpoint = cumulative - Percent / 2) %>%
  ungroup() %>%
  ggplot(aes(x = "", y = Percent, fill = Percent_Type)) +
  geom_bar(width = 1, stat = "identity", color = "black") + 
  coord_polar(theta = "y", start = 0) + 
  facet_wrap(~Probe, labeller = as_labeller(c("None" = "Uncaptured average", "JA2" = "Captured average"))) +
  scale_fill_manual(values = c("#00CED1", "#708090", "#E0D8B0")) + 
  geom_text_repel(aes(y = midpoint, label = paste(Percent_Type, "\n", scales::percent(Percent / 100))), size = 4, color = "black", box.padding = 0.3, force = 2, force_pull = 2, min.segment.length = 0.2, segment.size = 0.5) + 
  # labs(title = "AVERAGES THP1 cells spiked with 1e6 H37Ra") + 
  my_plot_themes
PieChart_Averages_fig1
# ggsave(PieChart_Averages_fig1,
#        file = "PieChart_Averages_fig1.pdf",
#        path = "Figures/Captured.vs.Not_ReadCounts",
#        width = 6, height = 4, units = "in")





