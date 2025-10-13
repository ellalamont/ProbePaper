# Figure 3


source("FINAL_ImportData.R")


###########################################################
###################### FIGURE 3A ##########################

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, size=14, vjust=1, hjust=1),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(2, 2, 2, 2))

my_fav_shapes <- c(`Sputum` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)
my_fav_colors <- c(`Sputum` = "#0072B2", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")

Fig3A <- BiolSamples_pipeSummary %>% 
  filter(Type != "Broth") %>%
  ggplot(aes(x = Type, y = N_Genomic_Rv)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type, shape = Type), alpha = 0.8, size = 2, position = position_jitter(0.2)) + 
  scale_shape_manual(values = my_fav_shapes) + 
  scale_fill_manual(values=my_fav_colors) +  
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  scale_y_continuous(limits = c(0,10000000), breaks = seq(0, 10000000, 2000000)) +
  labs(x = "Sample type", y = "# reads aligning to Mtb") + 
  my_plot_themes
Fig3A


###########################################################
###################### FIGURE 3B ##########################

Fig3B <- BiolSamples_pipeSummary %>% 
  filter(Type != "Broth") %>%
  ggplot(aes(x = Type, y = P_Genomic_Rv)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type, shape = Type), alpha = 0.8, size = 2, position = position_jitter(0.2)) + 
  scale_shape_manual(values = my_fav_shapes) + 
  scale_fill_manual(values=my_fav_colors) +  
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  labs(x = "Sample type", y = "% reads aligning to Mtb") + 
  my_plot_themes
Fig3B

###########################################################
###################### FIGURE 3C ##########################

Fig3C <- BiolSamples_pipeSummary %>% 
  filter(Type != "Broth") %>%
  ggplot(aes(x = Type, y = Txn_Coverage)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type, shape = Type), alpha = 0.8, size = 2, position = position_jitter(0.2)) + 
  scale_shape_manual(values = my_fav_shapes) + 
  scale_fill_manual(values=my_fav_colors) +  
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  geom_hline(yintercept = 80, linetype = "dashed", alpha = 0.5) + 
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) + 
  labs(x = "Sample type", y = "% transcriptional coverage") + 
  my_plot_themes
Fig3C

###########################################################
###################### FIGURE 3D ##########################

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(2, 2, 2, 2))

Fig3D <- BiolSamples_pipeSummary %>% 
  filter(Type %in% c("Marmoset", "Caseum mimic", "Rabbit")) %>% 
  ggplot(aes(x = CFU_per_g.or.mL, y = P_Genomic_Rv)) +
  geom_point(aes(fill = Type, shape = Type), size = 4, alpha = 0.8, stroke = 0.8, color = "black") + 
  scale_fill_manual(values = my_fav_colors) +  
  scale_shape_manual(values = my_fav_shapes) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) +
  scale_x_continuous(trans = "log10") + 
  labs(x = "log10(CFU/g)", y = "% reads aligning to Mtb") + 
  stat_poly_line(method = "lm", se = F, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "left", label.y = 0.99, parse = T, size = 4) + 
  my_plot_themes
Fig3D

###########################################################
###################### FIGURE 3E ##########################

Fig3E <- BiolSamples_pipeSummary %>% 
  filter(Type == "Sputum") %>%
  filter(XpertCT_wk0 != "NA") %>% 
  mutate(XpertCT_wk0 = as.numeric(XpertCT_wk0)) %>% 
  ggplot(aes(x = XpertCT_wk0, y = P_Genomic_Rv)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black", fill = "#0072B2", shape = 21) + 
  # geom_text_repel(aes(label = format(Patient, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) +
  scale_x_continuous(limits = c(12,24.5), breaks = seq(12,24,2), expand = c(0,0)) +
  labs(# title = "Sputum subset: Ct value vs percent reads aligned to Mtb",
    subtitle = NULL,
    y = "% reads aligning to Mtb",
    x = "Xpert Ct value") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = F, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "left", label.y = 0.99, parse = T, size = 4)
Fig3E

###########################################################
###################### FIGURE 3F ##########################

Fig3F <- BiolSamples_pipeSummary %>% 
  filter(Type == "Sputum") %>%
  filter(TTD != "NA") %>% 
  mutate(TTD = as.numeric(TTD)) %>% 
  ggplot(aes(x = TTD, y = P_Genomic_Rv)) +
  geom_point(size = 4, alpha = 0.8, stroke = 0.8, color = "black", fill = "#0072B2", shape = 21) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) +
  scale_x_continuous(limits = c(2,6), breaks = seq(2,6,1)) +
  labs(y = "% reads aligning to Mtb", x = "TTD (Day)") + 
  my_plot_themes + 
  stat_poly_line(method = "lm", se = F, level = 0.95, color = "grey23", alpha = 0.25) + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*")),
               label.x = "left", label.y = 0.99, parse = T, size = 4)
Fig3F

###########################################################
##################### COMBINE FIGURES #####################

combined <- plot_grid(Fig3A, Fig3B, Fig3C, Fig3D, Fig3E, Fig3F, align = "hv", axis = "tblr", nrow = 2, ncol = 3, rel_widths = c(1, 1, 1), rel_heights = c(1, 1))
combined
# ggsave(combined,
#        file = paste0("FINAL_Figure3.pdf"),
#        path = "Figures/CombinedFigures",
#        width = 15, height = 10, units = "in")






