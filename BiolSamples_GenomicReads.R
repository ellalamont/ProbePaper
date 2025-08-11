# Look at the P_Genomic and N_Genomic for the biological samples
# E. Lamont
# 7/31/25

source("Import_data.R") # To get BiolSamples_pipeSummary
# 8/4/25: Removed all samples with RawReads < 1M (these failed sequencing)


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_blank(), 
        # axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.text.x = element_text(angle = 45, size=14, vjust=1, hjust=1),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(2, 2, 2, 2)#
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

###########################################################
################ N_Genomic vs SAMPLE TYPE #################

### BOXPLOT ###
N_Genomic_box1 <- BiolSamples_pipeSummary %>% 
  filter(Type != "Broth") %>%
  ggplot(aes(x = Type, y = N_Genomic)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type, shape = Type), alpha = 0.8, size = 2, position = position_jitter(0.2)) + 
  scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = format(SampleID, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_fill_manual(values=my_fav_colors) +  
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,12000000), breaks = seq(0, 12000000, 2000000)) +
  labs(# title = "N_Genomic for all biological sample types",
       x = "Sample type", 
       y = "# reads aligning to Mtb") + 
  my_plot_themes
N_Genomic_box1
# ggsave(N_Genomic_box1,
#        file = paste0("N_Genomic_Box1.pdf"),
#        path = "Figures/GenomicRead_Analyses",
#        width = 5, height = 5, units = "in")

## GGERRORPLOT ###
# N_Genomic_errorplot1 <- BiolSamples_pipeSummary %>% 
#   filter(Type != "Broth") %>%
#   ggerrorplot(x = "Type", y = "N_Genomic", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.6,  # Size of error bars
#               add.params = list(size = 0.7)) +  # Size of mean points
#   geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
#   scale_y_continuous(limits = c(0,12000000), breaks = seq(0, 12000000, 2000000)) +
#   geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
#   labs(title = "N_Genomic for all biological sample types", 
#        subtitle = "Mean with standard deviation", 
#        x = "# spiked in H37Ra cells", 
#        y = "% reads aligning to Mtb") + 
#   my_plot_themes
# N_Genomic_errorplot1
# ggsave(N_Genomic_errorplot1,
#        file = paste0("N_Genomic_errorplot1.pdf"),
#        path = "Figures/GenomicRead_Analyses",
#        width = 7, height = 5, units = "in")


###########################################################
################ P_Genomic vs SAMPLE TYPE #################

P_Genomic_box1 <- BiolSamples_pipeSummary %>% 
  filter(Type != "Broth") %>%
  ggplot(aes(x = Type, y = P_Genomic)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type, shape = Type), alpha = 0.8, size = 2, position = position_jitter(0.2)) + 
  scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = format(SampleID, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_fill_manual(values=my_fav_colors) +  
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  labs(# title = "P_Genomic for all biological sample types",
       x = "Sample type", 
       y = "% reads aligning to Mtb") + 
  my_plot_themes
P_Genomic_box1
# ggsave(P_Genomic_box1,
#        file = paste0("P_Genomic_Box1.pdf"),
#        path = "Figures/GenomicRead_Analyses",
#        width = 5, height = 5, units = "in")


###########################################################
############ AtLeast.10.Reads vs SAMPLE TYPE ##############

TenReads_box1 <- BiolSamples_pipeSummary %>% 
  filter(Type != "Broth") %>%
  mutate(Txn_Coverage = (AtLeast.10.Reads/4499)*100) %>%
  ggplot(aes(x = Type, y = Txn_Coverage)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type, shape = Type), alpha = 0.8, size = 2, position = position_jitter(0.2)) + 
  scale_shape_manual(values = my_fav_shapes) + 
  scale_fill_manual(values=my_fav_colors) +  
  # geom_text_repel(aes(label = format(SampleID, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  geom_hline(yintercept = 80, linetype = "dashed", alpha = 0.5) + 
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) + 
  labs(# title = "Genes with >= 10 reads aligning for all biological sample types",
       # subtitle = "", 
       x = "Sample type", 
       y = "% transcriptional coverage") + 
  my_plot_themes
TenReads_box1
# ggsave(TenReads_box1,
#        file = paste0("TenReads_box1.pdf"),
#        path = "Figures/GenomicRead_Analyses",
#        width = 5, height = 5, units = "in")




