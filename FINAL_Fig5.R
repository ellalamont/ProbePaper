# Figure 5

source("FINAL_ImportData.R")

my_plot_themes <- theme_bw() +
  theme(legend.position = "right",legend.text=element_text(size=7),
        legend.title = element_text(size = 7),
        # legend.title = element_blank(),
        plot.title = element_text(size=7), 
        axis.title.x = element_text(size=7), 
        axis.text.x = element_text(angle = 0, size=7, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=7),
        axis.text.y = element_text(size=7), 
        plot.subtitle = element_text(size=7))
facet_themes <- theme(strip.background=element_rect(fill="white", linewidth = 0.9),
                      strip.text = element_text(size = 7))




SputumVsBroth_iModulons <- DEG_dfs$Sputum.vs.Broth_iModulons

SputumVsBroth_iModulons <- SputumVsBroth_iModulons %>%
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant"))

iModulons_of_interest <- c("PrpR", "BkaR", "Rv0681", "KstR2", "Fatty Acid Biosynthesis", "FasR", "Peptidoglycan Biosynthesis", "Positive Regulation of Growth", "Mycofactocin Synthesis Pathway", "DevR-2", "GroEL-GroES Complex", "DevR-1", "PhoP", "MprA", "Mce3R", "Mce1R", "SigC", "SigD", "SigH", "SigK", "RicR", "IdeR", "Zur", "WhiB1") 
iModulons_of_interest_pattern <- str_c(iModulons_of_interest, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

Fav_Pathways <- SputumVsBroth_iModulons %>% 
  filter(str_detect(PathName, iModulons_of_interest_pattern)) %>% 
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  pull(PathName)

# Put them in categories
SputumVsBroth_iModulons2 <- SputumVsBroth_iModulons %>%
  mutate(iModulonCategory2 = case_when(
    str_detect(PathName, "DevR") ~ "DosR", # specifically putt the DevR in a different category (normally in redox)
    str_detect(PathName, paste(Growth_iModulons_pattern, NucleicAcid_iModulons_pattern, Redox_iModulons_pattern, AminoAcid_iModulons_pattern, sep = "|")) ~ "Growth",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence and Persistence",
    str_detect(PathName, paste(CentralCarbon_iModulons_pattern, FattyAcid.Cholesterol_iModulons_pattern, sep = "|")) ~ "Fatty Acid and Cholesterol",
    TRUE ~ "Other"
  )) %>%
  mutate(iModulonCategory2 = factor(iModulonCategory2, levels = c("Fatty Acid and Cholesterol", "Growth", "DosR", "Metal","Virulence and Persistence", "Other")))


my_plot_themes2 <- theme_bw(base_family = "Arial") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=7),
        legend.title = element_text(size = 7),
        axis.ticks.y = element_blank(),
        # legend.title = element_blank(),
        plot.title = element_text(size=7),
        axis.title.x = element_text(size=12),
        axis.text.x = element_text(angle = 0, size=9, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=7),
        axis.text.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.subtitle = element_text(size=7))
facet_themes2 <- theme(strip.background=element_rect(fill="white", linewidth = 0.9),
                      strip.text = element_text(size = 12),
                      panel.border = element_rect(color = "black", fill = NA, size = 1))


Fig5 <- SputumVsBroth_iModulons2 %>%
  filter(PathName %in% Fav_Pathways) %>%
  filter(iModulonCategory2 != "Virulence and Persistence") %>% 
  filter(!str_detect(PathName, "WhiB4/IdeR")) %>% # To remove this iModulon which is tagging along with IdeR iModulon
  filter(N_Genes >=3) %>% 
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  ggplot(aes(x = LOG2FOLD, y = PathName_2, shape = Type)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = ifelse(FDR_Significance == "significant", ifelse(LOG2FOLD>0, "pos", "neg"), "ns")),
             size = 4, shape = 21, alpha = 0.8) + 
  scale_fill_manual(
    values = c("pos" = "#bb0c00", 
               "neg" = "#00AFBB", 
               "ns"  = "grey"),
    name = "Significance / Direction"
  ) +
  facet_grid(rows = vars(iModulonCategory2), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  scale_x_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(y = NULL, x = "Log2Fold change") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
Fig5
# ggsave(Fig5,
#        file = paste0("Fig5_Thumbnail2.png"),
#        path = "Figures",
#        dpi = 600,
#        width = 4, height = 5, units = "in")
