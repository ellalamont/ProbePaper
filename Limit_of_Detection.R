# Limit of Detection with the THP1 spiked samples graphs
# E. Lamont
# 7/31/25

source("Import_data.R") # to get LimitofDetect_pipeSummary


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
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )


###########################################################
################ CELL NUMBER VS N_GENOMIC #################

# ggerrorplot
LimitofDetect_NumReads_Fig1 <- LimitofDetect_pipeSummary %>% 
  ggerrorplot(x = "Ra_cells2", y = "N_Genomic", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4,  # Size of error bars
              add.params = list(size = 0.4)) +  # Size of mean points
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # annotate("text", x = 5.6, y = 1000000*0.8, label = "1 million", 
  #          hjust = 1.1, vjust = -0.5, color = "black") + 
  # scale_y_continuous(limits = c(0,19000000), breaks = seq(0, 19000000, 2000000)) + 
  scale_y_continuous(limits = c(0,20000000), breaks = seq(0, 20000000, 4000000)) + 
  labs(# title = "ProbeTest5 THP1 cells spiked with H37Ra", 
       # subtitle = "Mean with standard deviation", 
       x = "# spiked in H37Ra cells", 
       y = "# reads aligning to Mtb") + 
  my_plot_themes
LimitofDetect_NumReads_Fig1
# ggsave(LimitofDetect_NumReads_Fig1,
#        file = "N_Genomic_Limit.of.Detection_fig1.pdf",
#        path = "Figures/Limit_of_Detection",
#        width = 7, height = 5, units = "in")


###########################################################
################ CELL NUMBER VS P_GENOMIC #################

# ggerrorplot
LimitofDetect_PercentReads_Fig1 <- LimitofDetect_pipeSummary %>% 
  ggerrorplot(x = "Ra_cells2", y = "P_Genomic", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4,  # Size of error bars
              add.params = list(size = 0.4)) +  # Size of mean points
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  labs(# title = "ProbeTest5 THP1 cells spiked with H37Ra", 
       # subtitle = "Mean with standard deviation", 
       x = "# spiked in H37Ra cells", 
       y = "% reads aligning to Mtb") + 
  my_plot_themes
LimitofDetect_PercentReads_Fig1
# ggsave(LimitofDetect_PercentReads_Fig1,
#        file = "P_Genomic_Limit.of.Detection_fig1.pdf",
#        path = "Figures/Limit_of_Detection",
#        width = 7, height = 5, units = "in")


###########################################################
############# CELL NUMBER VS AtLeast.10.Reads #############

# ggerrorplot
LimitofDetect_10Reads_Fig1 <- LimitofDetect_pipeSummary %>% 
  mutate(Txn_Coverage = (AtLeast.10.Reads/4499)*100) %>%
  ggerrorplot(x = "Ra_cells2", y = "Txn_Coverage", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4,  # Size of error bars
              add.params = list(size = 0.4)) +  # Size of mean points
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  geom_hline(yintercept = 80, linetype = "dashed", alpha = 0.5) + 
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) + 
  labs(x = "# spiked in H37Ra cells", 
       y = "% transcriptional coverage") + 
  my_plot_themes
LimitofDetect_10Reads_Fig1
# ggsave(LimitofDetect_10Reads_Fig1,
#        file = "TenReads_Limit.of.Detection_fig1.pdf",
#        path = "Figures/Limit_of_Detection",
#        width = 7, height = 5, units = "in")








