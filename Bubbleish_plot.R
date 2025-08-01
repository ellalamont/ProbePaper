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
######## iMODULONS: SPUTUM and MARMOSET vs BROTH ##########

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

# Merge the sputum and the marmoset dataframes
merged_df <- full_join(SputumVsBroth_iModulons, MarmosetVsBroth_iModulons, by = c("PathName", "CellType", "N_Genes"))

# Pivot_longer the data (Chat GPT helped here)
merged_long <- merged_df %>%
  pivot_longer(cols = c(Sputum_LOG2FOLD, Marmoset_LOG2FOLD,
                        Sputum_AVG.PVALUE, Marmoset_AVG.PVALUE,
                        Sputum_AVG.RANK, Marmoset_AVG.RANK),
               names_to = c("Type", ".value"),
               names_sep = "_")

# Just taking the top 60 here, will need to subset better later
merged_subset <- merged_long %>% slice_head(n=60) %>% 
  mutate(PathName = str_wrap(PathName, width = 50))

merged_subset$Significance <- ifelse(merged_subset$AVG.PVALUE < 0.05, "significant", "not significant")


# Make the plot
# test <- merged_subset %>%
#   ggplot(aes(x = Type, y = PathName, fill = LOG2FOLD)) + 
#   geom_point(aes(fill = LOG2FOLD, shape = Type, stroke = ifelse(Significance == "significant", 0.8, 0)),
#                  size = 5, alpha = 0.7) +
#   scale_shape_manual(values = c(21, 21)) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
#   my_plot_themes
# test


test <- merged_subset %>%
  ggplot(aes(x = Type, y = PathName, fill = LOG2FOLD, shape = Type)) + 
  geom_point(data = filter(merged_subset, Significance == "significant"),
             size = 5, stroke = 0.8, alpha = 0.7) +
  geom_point(data = filter(merged_subset, Significance != "significant"),
             size = 4, stroke = 0, alpha = 0.7) +
  scale_shape_manual(values = c(21, 21)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  guides(shape = "none") + 
  my_plot_themes
test



