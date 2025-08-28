# Make a bubble plot sort of graph from Bob's metagenesets
# E. Lamont
# 8/1/25

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

SputumVsRv_iModulons <- MetaGeneSets_GoodSputumSubset.vs.Indigo_Rv %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(Sputum.Rv_FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% # Because Bob's pipeline doesn't actually do this...
  rename("Sputum.Rv_LOG2FOLD" = "LOG2FOLD", "Sputum.Rv_AVG.PVALUE" = "AVG_PVALUE", "Sputum.Rv_AVG.RANK" = "AVG_RANK") %>%
  mutate()


###########################################################
############### iMODULONS: MERGE DATAFRAMES ###############

# Merge the sputum and the marmoset dataframes
merged_df <- full_join(SputumVsRa_iModulons, SputumVsRv_iModulons, by = c("PathName", "CellType", "N_Genes"))

merged_long <- merged_df %>%
  pivot_longer(cols = -c(1:3),
               names_to = c("Type", ".value"),
               names_sep = "_") %>% 
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  mutate(Significance = ifelse(AVG.PVALUE < 0.05, "significant", "not significant")) %>%
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant"))

ordered_type <- c("Sputum", "Marmoset", "Rabbit", "Caseum mimic")
merged_long$Type <- factor(merged_long$Type, levels = ordered_type)

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
