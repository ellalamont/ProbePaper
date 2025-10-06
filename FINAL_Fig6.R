# Figure 6

source("FINAL_ImportData.R")

combined_lit_df <- DEG_dfs$Literature.Data


################################################
###### FUNCTION FROM COPPOLA ET AL. 2021 #######

quantile_normalisation <- function(df)
  
{df_rank <- apply (df,2, rank, ties.method="min")

df_sorted <- data.frame (apply (df, 2, sort))

df_mean <- apply (df_sorted, 1, mean)

index_to_mean <- function (my_index, my_mean) {return (my_mean [my_index])}

df_final <- apply (df_rank, 2, index_to_mean, my_mean = df_mean)

rownames (df_final) <- rownames(df)

return(df_final)}


#########################################################
######### DO THIS ALL WITH THE SPUTUM RAW READS #########
GoodSampleList_Sputum <- BiolSamples_pipeSummary %>%
  filter(N_Genomic_Rv >= 1000000 & Txn_Coverage >= 80) %>% 
  filter(Type == "Sputum") %>%
  pull(SampleID)

GoodSputum_RawReads <- All_RawReads %>% 
  select(all_of(c("X", GoodSampleList_Sputum))) %>% 
  column_to_rownames(("X"))

GoodSputumSubset_RawReads_QuantileNormalization <- as.data.frame(quantile_normalisation(GoodSputum_RawReads))

# Change the Quantile Normalization to a 0-100 scale
rank_cols_to_percentile <- function(df) {
  apply(df, 2, function(col) {
    ranks <- rank(col, ties.method = "min")
    percentiles <- (ranks - 1) / (length(ranks) - 1) * 100
    return(percentiles)
  }) |> as.data.frame()
}

GoodSputumSubset_RawReads_QuantileNormalization_0.100 <- rank_cols_to_percentile(GoodSputumSubset_RawReads_QuantileNormalization)

# Get the average for my 9 sputum samples
RANK_Average_RawReads_GoodSputumSubset <- GoodSputumSubset_RawReads_QuantileNormalization_0.100 %>% 
  mutate(RANK_Average = rowMeans(across(where(is.numeric)))) %>%
  select(RANK_Average)

###########################################################
##################### COMBINE THE DATA ####################

All_RankExpression <- full_join(RANK_Average_RawReads_GoodSputumSubset %>% rownames_to_column("Gene"), combined_lit_df,  by = "Gene") %>%
  rename(CurrentStudy = RANK_Average) %>%
  rename_with(~ str_replace(., "_.*", ""))

###########################################################
############### GGCORRPLOT NA ROWS REMOVED ################

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=16, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=16), 
        plot.subtitle = element_text(size=12))


# Remove any row that contains NA
All_RankExpression_noNA <- All_RankExpression %>% na.omit()

# Not sure why I had to process the data so much here
All_RankExpression_noNA <- All_RankExpression_noNA %>%
  as.data.frame()
rownames(All_RankExpression_noNA) <- NULL
All_RankExpression_noNA <- All_RankExpression_noNA %>%
  column_to_rownames(var = "Gene")

# Make the correlation
corr <- cor(All_RankExpression_noNA, method = "spearman")

# Order the columns
# my_order <- c("Walter2015", "Garcia2016", "Sharma2017", "Lai2021", "CurrentStudy")
# corr <- corr[my_order, my_order]

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(All_RankExpression_noNA, method = "spearman")
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
             legend.title = "Spearman\ncorrelation",
             outline.col = "white") + 
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Spearman Correlation Literature sputum. All NA rows removed", 
       subtitle = "n = 1932 genes, X = not significant", 
       x = NULL, y = NULL)
ggcorrplot_Spearman_noNA

