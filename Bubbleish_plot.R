# Make a bubble plot sort of graph from Bob's metagenesets
# E. Lamont


source("Import_DEG_sets.R")

# Here the gene set enrichment analysis has already been done in Bob's meta way and I am just visualizing the result

# Plot basics
my_plot_themes <- theme_bw() +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=9),
        legend.title = element_text(size = 10),
        # legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=8), 
        plot.subtitle = element_text(size=10), 
        plot.margin = margin(10, 10, 10, 20)# ,
  )



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
  # mutate(Type = "Sputum")
  rename("Sputum_LOG2FOLD" = "LOG2FOLD", "Sputum_AVG.PVALUE" = "AVG_PVALUE", "Sputum_AVG.RANK" = "AVG_RANK")

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
  # mutate(Type = "Marmoset")
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
  # mutate(Type = "Marmoset")
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
  # mutate(Type = "Marmoset")
  rename("Caseum mimic_LOG2FOLD" = "LOG2FOLD", "Caseum mimic_AVG.PVALUE" = "AVG_PVALUE", "Caseum mimic_AVG.RANK" = "AVG_RANK")

###########################################################
############### iMODULONS: MERGE DATAFRAMES ###############

# Merge the sputum and the marmoset dataframes
merged_df <- full_join(SputumVsBroth_iModulons, MarmosetVsBroth_iModulons, by = c("PathName", "CellType", "N_Genes"))
merged_df <- full_join(merged_df, RabbitVsBroth_iModulons, by = c("PathName", "CellType", "N_Genes"))
merged_df <- full_join(merged_df, MimicVsBroth_iModulons, by = c("PathName", "CellType", "N_Genes"))

# Pivot_longer the data (Chat GPT helped here)
# merged_long <- merged_df %>%
#   pivot_longer(cols = c(Sputum_LOG2FOLD, Marmoset_LOG2FOLD,
#                         Sputum_AVG.PVALUE, Marmoset_AVG.PVALUE,
#                         Sputum_AVG.RANK, Marmoset_AVG.RANK),
#                names_to = c("Type", ".value"),
#                names_sep = "_") %>%
#   mutate(Significance = ifelse(AVG.PVALUE < 0.05, "significant", "not significant"))

merged_long <- merged_df %>%
  pivot_longer(cols = -c(1:3),
               names_to = c("Type", ".value"),
               names_sep = "_") %>%
  mutate(Significance = ifelse(AVG.PVALUE < 0.05, "significant", "not significant"))

ordered_type <- c("Sputum", "Marmoset", "Rabbit", "Caseum mimic")
merged_long$Type <- factor(merged_long$Type, levels = ordered_type)

# Just taking the top 60 here, will need to subset better later
# merged_subset <- merged_long %>% slice_head(n=60) %>% 
#   mutate(PathName = str_wrap(PathName, width = 50))
# 
# merged_subset$Significance <- ifelse(merged_subset$AVG.PVALUE < 0.05, "significant", "not significant")

###########################################################
########### iMODULONS: CHOOSE GENE SETS TO SHOW ###########

# Get the most significant in Sputum and have that be the list for all because that's what we really care about
Fav_Pathways <- SputumVsBroth_iModulons %>% 
  # arrange(Sputum_AVG.PVALUE) %>% slice_head(n=20) %>% 
  filter(Sputum_AVG.PVALUE < 0.05) %>%
  pull(PathName)

merged_subset_2 <- merged_long %>% filter(PathName %in% Fav_Pathways) %>% mutate(PathName = str_wrap(PathName, width = 50))
# merged_subset_2$Significance <- ifelse(merged_subset_2$AVG.PVALUE < 0.05, "significant", "not significant")


###########################################################
################ iMODULONS: MAKE BUBBLE PLOT ##############

iModulons_bubble_1 <- merged_subset_2 %>%
  ggplot(aes(x = Type, y = PathName, fill = LOG2FOLD, shape = Type)) + 
  geom_point(data = filter(merged_subset_2, Significance == "significant"),
             size = 5, stroke = 0.8, alpha = 1) +
  geom_point(data = filter(merged_subset_2, Significance != "significant"),
             size = 4.5, stroke = 0, alpha = 1) +
  scale_shape_manual(values = c(21, 21, 21, 21)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(x = 0.5, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  guides(shape = "none") + 
  labs(title = "All iModulons that are significant in the Sputum",
       subtitle = "circles without outlines means not significant",
       y = NULL, x = NULL) + 
  my_plot_themes
iModulons_bubble_1
ggsave(iModulons_bubble_1,
       file = paste0("iModulons_bubble_1.pdf"),
       path = "Figures/Bubbles",
       width = 8, height = 12, units = "in")








# Make the plot
# test <- merged_subset %>%
#   ggplot(aes(x = Type, y = PathName, fill = LOG2FOLD)) + 
#   geom_point(aes(fill = LOG2FOLD, shape = Type, stroke = ifelse(Significance == "significant", 0.8, 0)),
#                  size = 5, alpha = 0.7) +
#   scale_shape_manual(values = c(21, 21)) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
#   my_plot_themes
# test


