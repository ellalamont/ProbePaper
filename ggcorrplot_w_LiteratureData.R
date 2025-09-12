# Correlations with the literature data from Coppola et al. 2021
# E. Lamont
# 8/14/25

source("Quantile_Normalization.R") # RANK_Average_RawReads_GoodSputumSubset
source("Import_literature_data.R") # combined_lit_df


# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2

# Coppola 2021 just uses Spearman's, so just do this one (is it because of the way they normalized?)
# Is it bad that there are lots of NAs? (Some of the sputum datasets have a lot fewer genes)

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
  rename_with(~ str_replace(., "_.*", ""))


###########################################################
############### GGCORRPLOT NA ROWS REMOVED ################

# From Coppola2021 spearman correlation analysis: "The entire datasets-comparison consisted of the ranked expression of Mtb genes (n=1813) that were commonly investigated in all nine Mtb transcriptomes."
# So I think I have to keep just the genes that are in ALL that datasets (remove any row that contains NA) 
# I think this is the best way to do it, but not sure.....

# Remove any row that contains NA
All_RankExpression_noNA <- All_RankExpression %>% na.omit()
# Now there are 1932 genes... still more than they have, but they had included non-sputum datasets so this may be where the descrepancy is....

# See how many genes (not including NA's are in each column)
All_RankExpression %>%
  summarise(across(everything(), ~sum(!is.na(.))))
# Gene PredictTB Walter2015 Garcia2016 Sharma2017 Lai2021
# 1 4597      4499       2406       1970       3924    4111

# Not sure why I had to process the data so much here
All_RankExpression_noNA <- All_RankExpression_noNA %>%
  as.data.frame()
rownames(All_RankExpression_noNA) <- NULL
All_RankExpression_noNA <- All_RankExpression_noNA %>%
  column_to_rownames(var = "Gene")

# Make the correlation
corr <- cor(All_RankExpression_noNA, method = "spearman")
# Needs the use argument so it will ignore the NAs between each pairwise comparison

# Order the columns
# my_order <- c("Walter2015", "Garcia2016", "Sharma2017", "Lai2021", "CurrentStudy")
# corr <- corr[my_order, my_order]

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(All_RankExpression_noNA, method = "spearman")
# head(p.mat[, 1:4])

p.mat.adj <- p.mat %>%
  as.vector() %>% 
  p.adjust(method = "fdr") %>%
  matrix(nrow = nrow(p.mat), ncol = ncol(p.mat))
rownames(p.mat.adj) <- rownames(p.mat)
colnames(p.mat.adj) <- colnames(p.mat)

# Plot
ggcorrplot_Spearman_noNA <- corr %>% 
  ggcorrplot(hc.order = F, 
             p.mat = p.mat.adj,
             lab = TRUE, lab_size = 6,
             type = c("lower"),
             # colors = c("#1B86B4", "white", "#AC204B"),
             legend.title = "Spearman\ncorrelation",
             outline.col = "white") + 
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Spearman Correlation Literature sputum. All NA rows removed", 
       subtitle = "n = 1932 genes, X = not significant", 
       x = NULL, y = NULL)
ggcorrplot_Spearman_noNA

ggsave(ggcorrplot_Spearman_noNA,
       file = "LiteratureData_Spearman_v1.pdf",
       path = "Figures/ggcorrplot",
       width = 7, height = 6, units = "in")
