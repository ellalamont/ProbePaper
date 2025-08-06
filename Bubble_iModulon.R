# Make a bubble plot sort of graph from Bob's metagenesets
# E. Lamont
# 8/1/25


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

SputumVsBroth_iModulons <- MetaGeneSets_GoodSputumSubset.vs.Broth_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(Sputum_FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% # Because Bob's pipeline doesn't actually do this...
  rename("Sputum_LOG2FOLD" = "LOG2FOLD", "Sputum_AVG.PVALUE" = "AVG_PVALUE", "Sputum_AVG.RANK" = "AVG_RANK") %>%
  mutate()


MarmosetVsBroth_iModulons <- MetaGeneSets_Marmoset.vs.Broth_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(Marmoset_FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>%
  rename("Marmoset_LOG2FOLD" = "LOG2FOLD", "Marmoset_AVG.PVALUE" = "AVG_PVALUE", "Marmoset_AVG.RANK" = "AVG_RANK")

RabbitVsBroth_iModulons <- MetaGeneSets_Rabbit.vs.Broth_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(Rabbit_FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>%
  rename("Rabbit_LOG2FOLD" = "LOG2FOLD", "Rabbit_AVG.PVALUE" = "AVG_PVALUE", "Rabbit_AVG.RANK" = "AVG_RANK")

MimicVsBroth_iModulons <- MetaGeneSets_CaseumMimic.vs.Broth_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(`Caseum mimic_FDR.pvalue`  = p.adjust(AVG_PVALUE, method = "fdr")) %>%
  rename("Caseum mimic_LOG2FOLD" = "LOG2FOLD", "Caseum mimic_AVG.PVALUE" = "AVG_PVALUE", "Caseum mimic_AVG.RANK" = "AVG_RANK")

###########################################################
############### iMODULONS: MERGE DATAFRAMES ###############

# Merge the sputum and the marmoset dataframes
merged_df <- full_join(SputumVsBroth_iModulons, MarmosetVsBroth_iModulons, by = c("PathName", "CellType", "N_Genes"))
merged_df <- full_join(merged_df, RabbitVsBroth_iModulons, by = c("PathName", "CellType", "N_Genes"))
merged_df <- full_join(merged_df, MimicVsBroth_iModulons, by = c("PathName", "CellType", "N_Genes"))

merged_long <- merged_df %>%
  pivot_longer(cols = -c(1:3),
               names_to = c("Type", ".value"),
               names_sep = "_") %>% 
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  mutate(Significance = ifelse(AVG.PVALUE < 0.05, "significant", "not significant")) %>%
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant"))

ordered_type <- c("Sputum", "Marmoset", "Rabbit", "Caseum mimic")
merged_long$Type <- factor(merged_long$Type, levels = ordered_type)


###########################################################
########### iMODULONS: CHOOSE GENE SETS TO SHOW ###########

# Get the most significant in Sputum and have that be the list for all because that's what we really care about
Fav_Pathways <- SputumVsBroth_iModulons %>% 
  # arrange(Sputum_AVG.PVALUE) %>% slice_head(n=20) %>% 
  filter(Sputum_AVG.PVALUE < 0.05) %>%
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  pull(PathName)


merged_subset_2 <- merged_long %>% filter(PathName %in% Fav_Pathways) 
# merged_subset_2$Significance <- ifelse(merged_subset_2$AVG.PVALUE < 0.05, "significant", "not significant")


###########################################################
################ iMODULONS: MAKE BUBBLE PLOT ##############

iModulons_bubble_1 <- merged_subset_2 %>%
  ggplot(aes(x = Type, y = PathName, fill = LOG2FOLD, shape = Type)) + 
  geom_point(aes(fill = LOG2FOLD, shape = Type, stroke = ifelse(Significance == "significant", 0.8, 0)), size = 5, alpha = 1) +
  # geom_point(data = filter(merged_subset_2, Significance == "significant"),
  #            size = 5, stroke = 0.8, alpha = 1) +
  # geom_point(data = filter(merged_subset_2, Significance != "significant"),
  #            size = 5, stroke = 0, alpha = 1) +
  scale_shape_manual(values = c(21, 21, 21, 21)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-4.2, 4.2)) +
  geom_text(aes(x = 0.5, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  guides(shape = "none") + 
  labs(title = "All iModulons that are significant in the Sputum",
       subtitle = "circles without outlines means not significant",
       y = NULL, x = NULL) + 
  my_plot_themes
iModulons_bubble_1
# ggsave(iModulons_bubble_1,
#        file = paste0("iModulons_bubble_1.pdf"),
#        path = "Figures/Bubbles/iModulons",
#        width = 8, height = 12, units = "in")


###########################################################
############ iMODULONS: MAKE LISTS OF GROUPS ##############


# Make new column with Groups, from Import_DEG_sets.R
merged_long <- merged_long %>%
  mutate(iModulonCategory = case_when(
    str_detect(PathName, CentralCarbon_iModulons_pattern) ~ "Central Carbon",
    str_detect(PathName, AminoAcid_iModulons_pattern) ~ "Amino Acid",
    str_detect(PathName, NucleicAcid_iModulons_pattern) ~ "Nucleic Acid",
    str_detect(PathName, FattyAcid.Cholesterol_iModulons_pattern) ~ "Fatty Acid/Cholesterol",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, SulfurMetabolism_iModulons_pattern) ~ "Sulfur",
    str_detect(PathName, Growth_iModulons_pattern) ~ "Growth",
    str_detect(PathName, Redox_iModulons_pattern) ~ "Redox",
    str_detect(PathName, AcidStress_iModulons_pattern) ~ "Acid Stress",
    str_detect(PathName, Antibiotic_iModulons_pattern) ~ "Antibiotic",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence/Persistence",
    TRUE ~ "Other"
  ))


###########################################################
############# iModulons: LOOP FOR ALL BUBBLES #############

BubblePlot_Function <- function(df, title_text = NULL) {
  ggplot(df, aes(x = Type, y = PathName, fill = LOG2FOLD, shape = Type)) + 
    geom_point(aes(fill = LOG2FOLD, shape = Type, stroke = ifelse(Significance == "significant", 0.8, 0)), size = 6, alpha = 1) +
    scale_shape_manual(values = c(21, 21, 21, 21)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-4.2, 4.2)) +
    geom_label(aes(x = 0.3, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.7, fill = "white", label.size = NA) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
    guides(shape = "none") + 
    labs(title = title_text,
         subtitle = "circles without outlines means not significant",
         y = NULL, x = NULL) + 
    my_plot_themes
}

# for (cat in unique(na.omit(merged_long$iModulonCategory))) {
#   df_cat <- merged_long %>% filter(iModulonCategory == cat)
#   
#   p <- BubblePlot_Function(df_cat, title_text = paste(cat, " iModulon"))
#   
#   ggsave(p,
#          file = paste0(cat,".pdf"),
#          path = "Figures/Bubbles/iModulons",
#          width = 7.5, height = 6, units = "in")
#   }




###########################################################
################# iModulons: FACETED BUBBLES ##############

# Get the most significant in Sputum and have that be the list for all because that's what we really care about
Fav_Pathways <- SputumVsBroth_iModulons %>% 
  filter(Sputum_AVG.PVALUE < 0.05) %>%
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  pull(PathName)

# First try make a graph with all the significant for sputum like bubble 1
merged_subset_2 <- merged_long %>% filter(PathName %in% Fav_Pathways) 

iModulons_bubble_facet <- merged_subset_2 %>%
  # filter(iModulonCategory %in% c("Acid Stress", "Fatty Acid/Cholesterol", "Metal", "Nucleic Acid", "Redox", "Sulfur", "Virulence/Persistence")) %>% 
  filter(iModulonCategory != "Other") %>%
  filter(N_Genes >=3) %>% 
  # mutate(PathName = sub(":.*", "", PathName)) %>%
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  ggplot(aes(x = Type, y = PathName_2, fill = LOG2FOLD, shape = Type)) + 
  geom_point(aes(fill = LOG2FOLD, shape = Type, stroke = ifelse(Significance == "significant", 0.8, 0)), size = 4, alpha = 1) + # unadjusted P-values
  scale_shape_manual(values = c(21, 21, 21, 21)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-4.2, 4.2)) +
  # facet_col(facets = ~iModulonCategory, scales = "free_y", space = "free") + 
  # facet_wrap(facets = ~iModulonCategory, scales = "free_y", ncol = 2) + 
  facet_grid(rows = vars(iModulonCategory), scales = "free_y", space = "free", labeller = labeller(iModulonCategory = label_wrap_gen(width = 10))) + 
  guides(shape = "none") + 
  labs(title = "iModulons significant in Sputum (N>=3, Other removed)",
       subtitle = "circles without outlines means not significant",
       y = NULL, x = NULL) + 
  my_plot_themes + facet_themes
iModulons_bubble_facet
ggsave(iModulons_bubble_facet,
       file = paste0("iModulons_bubble_facet_testing3.pdf"),
       path = "Figures/Bubbles/iModulons",
       width = 5.8, height = 9, units = "in")




# Combining some different things here, but all are the iModulons
merged_long2 <- merged_long %>% 
  mutate(Categories2 = case_when(
    str_detect(PathName, "MarR") ~ "Acid stress",
    str_detect(PathName, "Growth") ~ "Growth",
    str_detect(PathName, "IdeR|Zur|RicR") ~ "Metal",
    TRUE ~ "other"
  ))


###########################################################
########## iModulons: FACETED BUBBLES FDR ADJUSTED ########
# Look at the FDR adusted p-values!!!

# Get the most significant in Sputum and have that be the list for all because that's what we really care about
FDRAdjusted_Fav_Pathways <- SputumVsBroth_iModulons %>% 
  filter(Sputum_FDR.pvalue < 0.05) %>%
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  pull(PathName)

# First try make a graph with all the significant for sputum like bubble 1
merged_subset_FDRAdjusted <- merged_long %>% filter(PathName %in% FDRAdjusted_Fav_Pathways) 

iModulons_bubble_facet_FDR <- merged_subset_FDRAdjusted %>%
  filter(iModulonCategory != "Other") %>%
  filter(N_Genes >=3) %>% 
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  ggplot(aes(x = Type, y = PathName_2, fill = LOG2FOLD, shape = Type)) + 
  geom_point(aes(fill = LOG2FOLD, shape = Type, stroke = ifelse(FDR_Significance == "significant", 0.8, 0)), size = 4, alpha = 1) + # FDR adjusted P-values
  scale_shape_manual(values = c(21, 21, 21, 21)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-4.2, 4.2)) +
  facet_grid(rows = vars(iModulonCategory), scales = "free_y", space = "free", labeller = labeller(iModulonCategory = label_wrap_gen(width = 10))) + 
  guides(shape = "none") + 
  labs(title = "iModulons significant in Sputum (N>=3, Other removed)",
       subtitle = "FDR Adjusted P-values! circles without outlines means not significant",
       y = NULL, x = NULL) + 
  my_plot_themes + facet_themes
iModulons_bubble_facet_FDR
ggsave(iModulons_bubble_facet_FDR,
       file = paste0("iModulons_bubble_facet_testing3_FDR.pdf"),
       path = "Figures/Bubbles/iModulons",
       width = 5.8, height = 8.5, units = "in")


###########################################################
########### iModulons: MAKE INDIVIDUAL BUBBLES ############

### CENTRAL CARBON METABOLISM ###
# iModulons_bubble_CCM <- merged_long %>% 
#   filter(str_detect(PathName, CentralCarbon_iModulons_pattern)) %>% 
#   ggplot(aes(x = Type, y = PathName, fill = LOG2FOLD, shape = Type)) + 
#   geom_point(aes(fill = LOG2FOLD, shape = Type, stroke = ifelse(Significance == "significant", 0.8, 0)), size = 6, alpha = 1) +
#   scale_shape_manual(values = c(21, 21, 21, 21)) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-4.2, 4.2)) +
#   geom_label(aes(x = 0.3, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.7, fill = "white", label.size = NA) +
#   scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
#   guides(shape = "none") + 
#   labs(title = "CCM iModulons that are significant in the Sputum",
#        subtitle = "circles without outlines means not significant",
#        y = NULL, x = NULL) + 
#   my_plot_themes
# iModulons_bubble_CCM
# ggsave(iModulons_bubble_CCM,
#        file = paste0("iModulons_bubble_CCM.pdf"),
#        path = "Figures/Bubbles/iModulons",
#        width = 7.5, height = 6, units = "in")
# 
# ### AMINO ACID ###
# iModulons_bubble_AA <- merged_long %>% 
#   filter(str_detect(PathName, AminoAcid_iModulons_pattern)) %>% 
#   ggplot(aes(x = Type, y = PathName, fill = LOG2FOLD, shape = Type)) + 
#   geom_point(aes(fill = LOG2FOLD, shape = Type, stroke = ifelse(Significance == "significant", 0.8, 0)), size = 6, alpha = 1) +
#   scale_shape_manual(values = c(21, 21, 21, 21)) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-4.2, 4.2)) +
#   geom_label(aes(x = 0.3, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.7, fill = "white", label.size = NA) +
#   scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
#   guides(shape = "none") + 
#   labs(title = "Amino Acid iModulons that are significant in the Sputum",
#        subtitle = "circles without outlines means not significant",
#        y = NULL, x = NULL) + 
#   my_plot_themes
# iModulons_bubble_AA
# ggsave(iModulons_bubble_AA,
#        file = paste0("iModulons_bubble_AA.pdf"),
#        path = "Figures/Bubbles/iModulons",
#        width = 7.5, height = 6, units = "in")
# 
# 


