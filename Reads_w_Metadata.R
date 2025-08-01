# Compare the % and number reads to the sample metadata
# E. Lamont
# 7/31/25

# source("Import_data.R") To get SputumSubset_pipeSummary and BiolSamples_pipeSummary
source("Import_SampleMetadata.R")


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=14, vjust=1, hjust=1),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

my_regression_line <- stat_poly_line(method = "lm", se = TRUE, level = 0.95, color = "black", alpha = 0.3)
my_regression_equations <- stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                                          after_stat(rr.label),
                                                          after_stat(p.value.label), 
                                                          sep = "*\", \"*")))

###########################################################
####################### MERGE METADATA ####################

# For the Sputum
# Starting with SputumSubset_pipeSummary and my_metadata

# my_metadata is already a subset of the columns that I want to work on now

my_metadata_2 <- my_metadata %>% 
  filter(Patient %in% SputumSubset_pipeSummary$Patient) %>%
  filter(Visit == "Day 0")
# Missing P_14005 because this isn't in the subMIC sample set so I don't have Ct value for it

my_pipeSummary <- SputumSubset_pipeSummary %>%
  filter(Type == "Week 0 sputum") %>%
  select(SampleID, total_RNA_ng, Week, RawReads, N_RiboClear, P_RiboClear, N_Genomic, P_Genomic, N_NoHit, P_NoHit, AtLeast.10.Reads, AtLeast.100.Reads, Run, Type, Type2, SampleID2, Patient, Outcome)


merged_metadata <- merge(my_metadata_2, my_pipeSummary, by = "Patient", all = T)



###########################################################
################# SPUTUM N+P GENOMIC TO Ct ################

# N_GENOMIC
ctVsReads <- merged_metadata %>% 
  filter(XpertCT_wk0 != "NA") %>% 
  mutate(XpertCT_wk0 = as.numeric(XpertCT_wk0)) %>% 
  ggplot(aes(x = XpertCT_wk0, y = N_Genomic)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,12000000), breaks = seq(0,12000000, 2000000), labels = scales::scientific ) +
  scale_x_continuous(limits = c(12,31), breaks = seq(12,31,2), expand = c(0,0)) +
  labs(title = "Sputum subset: Ct value vs number reads aligned to Mtb",
       subtitle = NULL,
       y = "# reads aligning to Mtb genome") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = F, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ctVsReads # + my_regression_line + my_regression_equations
ggsave(ctVsReads,
       file = paste0("Sputum_ctVsReads_1.pdf"),
       path = "Figures/Reads_w_Metadata",
       width = 7, height = 5, units = "in")


# P_GENOMIC
ctVsPercent <- merged_metadata %>% 
  filter(XpertCT_wk0 != "NA") %>% 
  mutate(XpertCT_wk0 = as.numeric(XpertCT_wk0)) %>% 
  ggplot(aes(x = XpertCT_wk0, y = P_Genomic)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  geom_text_repel(aes(label = format(Patient, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) +
  scale_x_continuous(limits = c(12,31), breaks = seq(12,31,2), expand = c(0,0)) +
  labs(title = "Sputum subset: Ct value vs percent reads aligned to Mtb",
       subtitle = NULL,
       y = "% reads aligning to Mtb genome") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = F, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ctVsPercent # + my_regression_line + my_regression_equations
ggsave(ctVsPercent,
       file = paste0("Sputum_ctVsPercent_2.pdf"),
       path = "Figures/Reads_w_Metadata",
       width = 7, height = 5, units = "in")

###########################################################
################ SPUTUM N+P GENOMIC TO TTD ################

# N_GENOMIC
ttdVsReads <- merged_metadata %>% 
  filter(TTD != "NA") %>% 
  mutate(TTD = as.numeric(TTD)) %>% 
  ggplot(aes(x = TTD, y = N_Genomic)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,12000000), breaks = seq(0,12000000, 2000000), labels = scales::scientific ) +
  scale_x_continuous(limits = c(2,16), breaks = seq(2,16,2), expand = c(0,0)) +
  labs(title = "Sputum subset: TTD vs number reads aligned to Mtb",
       subtitle = NULL,
       y = "# reads aligning to Mtb") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = F, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ttdVsReads # + my_regression_line + my_regression_equations
ggsave(ttdVsReads,
       file = paste0("Sputum_ttdVsReads_1.pdf"),
       path = "Figures/Reads_w_Metadata",
       width = 7, height = 5, units = "in")

# P_GENOMIC
ttdVsPercent <- merged_metadata %>% 
  filter(TTD != "NA") %>% 
  mutate(TTD = as.numeric(TTD)) %>% 
  ggplot(aes(x = TTD, y = P_Genomic)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) +
  scale_x_continuous(limits = c(2,16), breaks = seq(0,16,2), expand = c(0,0)) +
  labs(title = "Sputum subset: TTD vs percent reads aligned to Mtb",
       subtitle = NULL,
       y = "percent reads aligning to Mtb") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = F, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T)
ttdVsPercent # + my_regression_line + my_regression_equations
ggsave(ttdVsPercent,
       file = paste0("Sputum_ttdVsPercent_1.pdf"),
       path = "Figures/Reads_w_Metadata",
       width = 7, height = 5, units = "in")


###########################################################
###################### CFU vs GENOMIC #####################

# N_GENOMIC
CFU.Reads <- BiolSamples_pipeSummary %>% 
  filter(Type %in% c("Marmoset", "Caseum mimic", "Rabbit")) %>% 
  ggplot(aes(x = CFU_per_g.or.mL, y = N_Genomic)) +
  geom_point(aes(fill = Type, shape = Type), size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = my_fav_colors) +  
  scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,10000000), breaks = seq(0, 10000000, 1000000)) +
  # scale_x_continuous(limits = c(12,31), breaks = seq(12,31,2), expand = c(0,0)) +
  scale_x_continuous(trans = "log10") + 
  labs(title = "Marmoset, rabbit, caseum mimic: CFU/g (/mL for mimic) vs number reads aligned to Mtb",
       subtitle = "Marmoset data frome ProbeTest3; two rabbits are replicates from same lesion",
       x = "log10(CFU/g or /mL)", 
       y = "# reads aligning to Mtb") + 
  stat_poly_line(method = "lm", se = F, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T) + 
  my_plot_themes
CFU.Reads
ggsave(CFU.Reads,
       file = paste0("CFU.Reads_1.pdf"),
       path = "Figures/Reads_w_Metadata",
       width = 8, height = 5, units = "in")

# P_GENOMIC
CFU.Percent <- BiolSamples_pipeSummary %>% 
  filter(Type %in% c("Marmoset", "Caseum mimic", "Rabbit")) %>% 
  ggplot(aes(x = CFU_per_g.or.mL, y = P_Genomic)) +
  geom_point(aes(fill = Type, shape = Type), size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = my_fav_colors) +  
  scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = format(N_Genomic, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) +
  # scale_x_continuous(limits = c(12,31), breaks = seq(12,31,2), expand = c(0,0)) +
  scale_x_continuous(trans = "log10") + 
  labs(title = "Marmoset, rabbit, caseum mimic: CFU/g (/mL for mimic) vs percent reads aligned to Mtb",
       subtitle = "Marmoset data frome ProbeTest3; two rabbits are replicates from same lesion",
       x = "log10(CFU/g or /mL)", 
       y = "% reads aligning to Mtb") + 
  stat_poly_line(method = "lm", se = F, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "right", label.y = "top", parse = T) + 
  my_plot_themes
CFU.Percent # + my_regression_line + my_regression_equations
ggsave(CFU.Percent,
       file = paste0("CFU.Percent_1.pdf"),
       path = "Figures/Reads_w_Metadata",
       width = 8, height = 5, units = "in")
