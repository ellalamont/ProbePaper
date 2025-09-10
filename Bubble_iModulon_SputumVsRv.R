# Make a bubble plot sort of graph from Bob's metagenesets
# E. Lamont
# 8/28/25

# Just comparing the Good sputum subset to the Indigo no drug control Rv to see what is going on with the PhoP iModulon. And comparing to the sputum vs Ra

source("Import_DEG_sets.R")
source("Import_GeneSets.R")

# Here the gene set enrichment analysis has already been done in Bob's meta way and I am just visualizing the result


# Plot basics
my_plot_themes <- theme_bw() +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=7),
        legend.title = element_text(size = 7),
        # legend.title = element_blank(),
        plot.title = element_text(size=7), 
        axis.title.x = element_text(size=7), 
        axis.text.x = element_text(angle = 0, size=7, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=7, vjust=1, hjust=1),
        axis.title.y = element_text(size=7),
        axis.text.y = element_text(size=7), 
        plot.subtitle = element_text(size=7)# , 
        # plot.margin = margin(10, 10, 10, 20)# ,
  )

facet_themes <- theme(strip.background=element_rect(fill="white", linewidth = 0.9),
                      strip.text = element_text(size = 7))


###########################################################
######## iMODULONS: COLLECT AND ORGANIZE SAMPLES ##########

# I know the UP and DOWN files are a little different, just working with the UP files for now, should maybe ask Bob if that is okay at some point


SputumVsRa_iModulons <- MetaGeneSets_GoodSputumSubset.vs.Broth_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(Sputum.Ra_FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% # Because Bob's pipeline doesn't actually do this...
  rename("Sputum.Ra_LOG2FOLD" = "LOG2FOLD", "Sputum.Ra_AVG.PVALUE" = "AVG_PVALUE", "Sputum.Ra_AVG.RANK" = "AVG_RANK") %>%
  mutate()

SputumVsINDIGORv_iModulons <- MetaGeneSets_GoodSputumSubset.vs.Indigo_Rv %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(Sputum.INDIGO_Rv_FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% # Because Bob's pipeline doesn't actually do this...
  rename("Sputum.INDIGO_Rv_LOG2FOLD" = "LOG2FOLD", "Sputum.INDIGO_Rv_AVG.PVALUE" = "AVG_PVALUE", "Sputum.INDIGO_Rv_AVG.RANK" = "AVG_RANK") %>%
  mutate()

SputumVsLanceRv_iModulons <- MetaGeneSets_GoodSputumSubset.vs.LancepH7_Rv_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(Sputum.LancepH7Rv_FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% # Because Bob's pipeline doesn't actually do this...
  rename("Sputum.LancepH7Rv_LOG2FOLD" = "LOG2FOLD", "Sputum.LancepH7Rv_AVG.PVALUE" = "AVG_PVALUE", "Sputum.LancepH7Rv_AVG.RANK" = "AVG_RANK") %>%
  mutate()

IndigoRvVsRa_iModulons <- MetaGeneSets_Indigo_Rv.vs.Ra_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(INDIGORv.Ra_FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% # Because Bob's pipeline doesn't actually do this...
  rename("INDIGORv.Ra_LOG2FOLD" = "LOG2FOLD", "INDIGORv.Ra_AVG.PVALUE" = "AVG_PVALUE", "INDIGORv.Ra_AVG.RANK" = "AVG_RANK") %>%
  mutate()

Lanceph7RvVsRa_iModulons <- MetaGeneSetsLancepH7_Rv.vs.Ra_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(LancepH7Rv.Ra_FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% # Because Bob's pipeline doesn't actually do this...
  rename("LancepH7Rv.Ra_LOG2FOLD" = "LOG2FOLD", "LancepH7Rv.Ra_AVG.PVALUE" = "AVG_PVALUE", "LancepH7Rv.Ra_AVG.RANK" = "AVG_RANK") %>%
  mutate()

Lanceph7RvVsIndigoRv_iModulons <- MetaGeneSets_LancepH7_Rv.vs.Indigo_Rv_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(LancepH7Rv.IndigoRv_FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% # Because Bob's pipeline doesn't actually do this...
  rename("LancepH7Rv.IndigoRv_LOG2FOLD" = "LOG2FOLD", "LancepH7Rv.IndigoRv_AVG.PVALUE" = "AVG_PVALUE", "LancepH7Rv.IndigoRv_AVG.RANK" = "AVG_RANK") %>%
  mutate()


###########################################################
############### iMODULONS: MERGE DATAFRAMES ###############

# Merge the sputum and the marmoset dataframes
merged_df <- full_join(SputumVsRa_iModulons, SputumVsRv_iModulons, by = c("PathName", "CellType"))
merged_df <- full_join(merged_df, SputumVsLanceRv_iModulons, by = c("PathName", "CellType"))
merged_df <- full_join(merged_df, IndigoRvVsRa_iModulons, by = c("PathName", "CellType"))
merged_df <- full_join(merged_df, Lanceph7RvVsRa_iModulons, by = c("PathName", "CellType"))
merged_df <- full_join(merged_df, Lanceph7RvVsIndigoRv_iModulons, by = c("PathName", "CellType"))

merged_long <- merged_df %>%
  pivot_longer(cols = -c(1:3, 8, 13, 18, 23, 28), # Have to take out all the extra N_Gene columns too...
               names_to = c("Type", ".value"),
               names_sep = "_") %>% 
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  mutate(Significance = ifelse(AVG.PVALUE < 0.05, "significant", "not significant")) %>%
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant"))

# Make new column with Groups, from Import_DEG_sets.R
merged_long <- merged_long %>%
  mutate(iModulonCategory = case_when(
    str_detect(PathName, CentralCarbon_iModulons_pattern) ~ "Central Carbon",
    str_detect(PathName, AminoAcid_iModulons_pattern) ~ "Amino Acid",
    str_detect(PathName, NucleicAcid_iModulons_pattern) ~ "Nucleic Acid",
    str_detect(PathName, FattyAcid.Cholesterol_iModulons_pattern) ~ "Fatty Acid_Cholesterol",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, SulfurMetabolism_iModulons_pattern) ~ "Sulfur",
    str_detect(PathName, Growth_iModulons_pattern) ~ "Growth",
    str_detect(PathName, Redox_iModulons_pattern) ~ "Redox",
    str_detect(PathName, AcidStress_iModulons_pattern) ~ "Acid Stress",
    str_detect(PathName, Antibiotic_iModulons_pattern) ~ "Antibiotic",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence_Persistence",
    TRUE ~ "Other"
  ))



###########################################################
############### FDR ADJUSTED SPECIFIC GROUPS ##############

# Lipid Metabolism and CCM: PrpR, BkaR, Rv0681, KstR2, Fatty Acid Biosynthesis, FasR, Peptidoglycan Biosynthesis
# Growth: Positive Regulation of Growth, Mycofactocin Synthesis Pathway, Nucleic Acid Hydrolysis, DevR-2, DevR-1
# Virulence/Persistence: All of them?
# Metal: Zur, RicR, IdeR

iModulons_of_interest <- c("PrpR", "BkaR", "Rv0681", "KstR2", "Fatty Acid Biosynthesis", "FasR", "Peptidoglycan Biosynthesis", "Positive Regulation of Growth", "Mycofactocin Synthesis Pathway", "DevR-2", "GroEL-GroES Complex", "DevR-1", "PhoP", "MprA", "Mce3R", "Mce1R", "SigC", "SigD", "SigH", "SigK", "RicR", "IdeR", "Zur", "WhiB1") 
iModulons_of_interest_pattern <- str_c(iModulons_of_interest, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

Fav_Pathways <- merged_long %>% 
  filter(str_detect(PathName, iModulons_of_interest_pattern)) %>% 
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  pull(PathName)

# Make new column with Groups, from Import_DEG_sets.R
merged_long2 <- merged_long %>%
  mutate(iModulonCategory2 = case_when(
    str_detect(PathName, paste(Growth_iModulons_pattern, Redox_iModulons_pattern, NucleicAcid_iModulons_pattern, AminoAcid_iModulons_pattern, sep = "|")) ~ "Growth",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence and Persistence",
    str_detect(PathName, paste(CentralCarbon_iModulons_pattern, FattyAcid.Cholesterol_iModulons_pattern, sep = "|")) ~ "Fatty Acid and Cholesterol",
    TRUE ~ "Other"
  ))


iModulons_newPathways <- merged_long2 %>%
  filter(PathName %in% Fav_Pathways) %>%
  filter(N_Genes.x >=3) %>% 
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes.x, ")")) %>%
  ggplot(aes(x = Type, y = PathName_2, fill = LOG2FOLD, shape = Type)) + 
  geom_point(aes(fill = LOG2FOLD, shape = Type, stroke = ifelse(FDR_Significance == "significant", 0.8, 0)), size = 4, alpha = 1) + # FDR adjusted P-values
  scale_shape_manual(values = c(21, 21, 21, 21, 21, 21)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_grid(rows = vars(iModulonCategory2), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  labs(title = "iModulons (N>=3, Other removed)",
       subtitle = "FDR Adjusted P-values! circles without outlines means not significant",
       y = NULL, x = NULL) + 
  my_plot_themes + facet_themes
iModulons_newPathways
# ggsave(iModulons_newPathways,
#        file = paste0("Sputum_vs_Ra.Rv.pdf"),
#        path = "Figures_preNonCodingRemoval/Bubbles/iModulons/FDR",
#        width = 6.5, height = 8.5, units = "in")


###########################################################
################# JUST LANCE Rv VS SPUTUM #################

iModulons_SputumVsLanceRv <- merged_long2 %>%
  filter(Type == "Sputum.LancepH7Rv") %>%
  filter(PathName %in% Fav_Pathways) %>%
  # filter(iModulonCategory2 != "Virulence and Persistence") %>% # Removing this because using Ra!!
  filter(!str_detect(PathName, "WhiB4/IdeR")) %>% # To remove this iModulon which is tagging along with IdeR iModulon
  filter(N_Genes.y >=3) %>% 
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes.y, ")")) %>%
  ggplot(aes(x = LOG2FOLD, y = PathName_2, shape = Type)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = ifelse(FDR_Significance == "significant", ifelse(LOG2FOLD>0, "pos", "neg"), "ns")),
             size = 4, shape = 21, alpha = 0.8) + # FDR adjusted P-values
  # scale_alpha_manual(values = c("significant" = 1, "not significant" = 0.7)) + 
  # scale_shape_manual(values = c(21, 21, 21, 21)) +
  scale_fill_manual(
    values = c("pos" = "#bb0c00", 
               "neg" = "#00AFBB", 
               "ns"  = "grey"),
    name = "Significance / Direction"
  ) +
  # scale_fill_manual(values = c("TRUE" = "#bb0c00", "FALSE" = "#00AFBB")) + 
  facet_grid(rows = vars(iModulonCategory2), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  # scale_x_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = "iModulons Sputum only vs LanceRv7",
       subtitle = "FDR Adjusted P-values! circles without outlines means not significant",
       y = NULL, 
       x = "Log2Fold change") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
iModulons_SputumVsLanceRv
ggsave(iModulons_SputumVsLanceRv,
       file = paste0("Sputum_vs_LanceRv.pdf"),
       path = "Figures_preNonCodingRemoval/Bubbles/iModulons/FDR",
       width = 6.5, height = 8.5, units = "in")




