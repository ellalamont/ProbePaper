# Figure 4

source("FINAL_ImportData.R")

###########################################################
###################### FIGURE 4A ##########################

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20))

my_fav_colors2 <- c(`Sputum L4` = "#004C73", `Sputum L2` = "#99CCE8", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00") 
my_fav_shapes2 <- c(`Sputum L4` = 21, `Sputum L2` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)

GoodSampleList <- BiolSamples_pipeSummary %>%
  filter(N_Genomic_Rv >= 1000000 & Txn_Coverage >= 80) %>% 
  pull(SampleID)

# Just keep the samples passing filter
GoodBiolSamples_tpm <- All_tpm %>% select(all_of(GoodSampleList))
my_tpm <- GoodBiolSamples_tpm 

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 34.7% of variance
summary_PCA[2,1] # PC2 explains 12.6% of variance
summary_PCA[3,1] # PC3 explains 9.0% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, BiolSamples_pipeSummary, by = "SampleID", )

Fig4A <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type3, shape = Type3)) + 
  geom_point(aes(fill = Type3, shape = Type3), size = 5, alpha = 0.7, stroke = 0.8) +
  # geom_text_repel(aes(label = Lineage), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors2) +  
  scale_shape_manual(values = my_fav_shapes2) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(# title = "PCA: >1M reads and >80% genes with at least 10 reads",
    # subtitle = "TPM",
    x = paste0("PC1: ", summary_PCA[1,1], "%"),
    y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
Fig4A
# ggsave(Fig4A,
#        file = paste0("Fig4A.pdf"),
#        path = "Figures/PCA/",
#        width = 6, height = 6, units = "in")


###########################################################
###################### FIGURE 4B ##########################

# Make an empty list to hold the count summaries
count_list <- list()

# Loop through list of dataframes
for (list_name in names(DEG_dfs)) {
  
  # Get the current dataframe
  df <- DEG_dfs[[list_name]]
  
  # Count the DE values (Significant, not, etc)
  GeneCounts <- df %>%
    count(DE, name = "Count") %>% # Does the actual counting 
    pivot_wider(names_from = DE, values_from = Count, values_fill = 0)
  
  # Add back the name of the current list
  GeneCounts <- GeneCounts %>%
    mutate(DEG_Group = list_name) %>%
    select(DEG_Group, everything()) # Puts list name first
  
  # Add the results to my list
  count_list[[list_name]] <- GeneCounts
  
}

# Convert the list of dataframes into one dataframe
DEG_count_df <- bind_rows(count_list)

DEG_count_df <- DEG_count_df %>% mutate(All_Signifiant_DEG = `significant down` + `significant up`)



###########################################################
###################### FIGURE 4C ##########################

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=10), 
        plot.margin = margin(10, 10, 10, 20))

make_volcano_function <- function(my_df, graph_title) {
  
  my_volcano <- my_df %>%
    ggplot(aes(x = LOG2FOLD, y = -log10(FDR_PVALUE), col = DE, label = Labels, text = GENE_NAME, label2 = GENE_ID)) + # text is for plotly, could be GENE_ID
    geom_point(alpha = 0.7) + 
    labs(title = graph_title) + 
    geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
    geom_text_repel(max.overlaps = 10, size = 3) +  
    
    labs(title = NULL, y = "-log10(FDR-adjusted p-value)") +
  
    scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00")) # +
  
  
  # Determine the max and min axes values for labeling 
  plot_build <- ggplot_build(my_volcano)
  y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
  x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
  x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
  
  # Add the gene number annotations
  text_up <- my_df %>% filter(DE == "significant up") %>% nrow()
  text_down <- my_df %>% filter(DE == "significant down") %>% nrow()
  my_volcano_annotated <- my_volcano +
    annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
    annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
  
  final_volcano <- my_volcano_annotated + my_plot_themes
  
}

Sputum.vs.Broth_volcano <- make_volcano_function(DEG_dfs[[1]], names(DEG_dfs)[1])
Sputum.vs.Broth_volcano
# ggsave(Sputum.vs.Broth_volcano,
#        file = paste0("FINAL_Figure4C.png"),
#        path = "Figures/CombinedFigures",
#        dpi = 600,
#        width = 6, height = 4, units = "in")

###########################################################
################# SUPPLEMENTAL FIGURE 2 ###################

L2.vs.L4_Volcano <- make_volcano_function(DEG_dfs[[11]], names(DEG_dfs)[11])
L2.vs.L4_Volcano



