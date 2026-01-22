# Graph of the gene numbers for the literature data
# E. Lamont
# 1/21/26


source("Quantile_Normalization.R") # RANK_Average_RawReads_GoodSputumSubset
source("Import_literature_data.R") # combined_lit_df


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=16, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=16), 
        plot.subtitle = element_text(size=12)# , 
        # plot.margin = margin(10, 10, 10, 20),
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

###########################################################
##################### COMBINE THE DATA ####################

All_RankExpression <- full_join(RANK_Average_RawReads_GoodSputumSubset %>% rownames_to_column("Gene"),
                                combined_lit_df,
                                by = "Gene") %>%
  rename(CurrentStudy = RANK_Average) %>%
  rename_with(~ str_replace(., "_.*", "")) %>%
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", Gene)) %>% # Keep only the protein coding Rv genes
  mutate(Lai2021 = na_if(as.numeric(Lai2021), 0))

GeneCounts <- All_RankExpression %>%
  summarize(across(-Gene, ~ sum(!is.na(.x)))) %>%
  pivot_longer(cols = everything(), names_to = "Study", values_to = "GeneCount") %>%
  mutate(Type = c("RNAseq", "RT-qPCR", "RT-qPCR", "Microarray", "RNAseq")) # add these by hand and check they are correct! 

Study_ordered <- c("Walter2015", "Garcia2016", "Sharma2017", "Lai2021", "CurrentStudy")
GeneCounts$Study <- factor(GeneCounts$Study, levels = Study_ordered)

Gene_plot <- GeneCounts %>%
  mutate(Type = factor(Type)) %>% 
  ggplot(aes(x = reorder(Study, GeneCount), y = GeneCount, fill = Type)) + 
  geom_col(position = position_dodge()) + 
  labs(x = "Study", y = "# of genes measured") + 
  my_plot_themes
Gene_plot






