# Make a volcano plot 
# 7/31/25

source("Import_DEG_sets.R") 



# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=10), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank())


###########################################################
########### FUNCTION TO MAKE ALL VOLCANO PLOTS ############

make_volcano_function <- function(my_df, graph_title) {
  
  ## Make a volcano plot using output from Bob's pipeline
  
  my_volcano <- my_df %>%
    ggplot(aes(x = LOG2FOLD, y = -log10(AVG_PVALUE), col = DE, label = DE_labels, text = GENE_NAME, label2 = GENE_ID)) + # text is for plotly, could be GENE_ID
    geom_point(alpha = 0.7) + 
    labs(title = graph_title) + 
    geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
    geom_text_repel(max.overlaps = 10, size = 3) +  # Can do geom_text_repel or geom_label_rebel
    
    # Need it this way so the colors aren't messed up by not having significant up or down
    # scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) + 
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

###########################################################
############### MAKE A SINGLE VOLCANO PLOT ################

my_path <- "Figures/Volcano_plot"
single_plot <- make_volcano_function(list_dfs_2[[1]], df_names[1])
single_plot
# ggsave(single_plot,
#        file = paste0(df_names[3], ".pdf"),
#        path = my_path,
#        width = 6, height = 4, units = "in")
# ggplotly(single_plot)



###########################################################
################# LOOP ALL VOLCANO PLOTS ##################

# df_names # Remember I have this

my_path <- "Figures/Volcano_plot/Log2Fold_at_1"

for (i in 1:length(list_dfs_2)) {

  current_df_name <- df_names[i]
  filename <- paste0(current_df_name, ".pdf")

  my_plot <- make_volcano_function(list_dfs_2[[i]], df_names[i])

  ggsave(my_plot,
         file = filename,
         path = my_path,
         width = 7, height = 5, units = "in")
}

# Make sure to check what the warnings are... okay if they are just the geom_text_repel but make sure not data is being lost by defined axes!


###########################################################
######### FUNCTION FOR VOLCANO WITH LOG2FOLD OF 2 #########

make_volcano_function_2 <- function(my_df, graph_title) {
  
  ## Make a volcano plot using output from Bob's pipeline
  
  my_volcano <- my_df %>%
    ggplot(aes(x = LOG2FOLD, y = -log10(AVG_PVALUE), col = DE_2, label = DE_2_labels, text = GENE_NAME, label2 = GENE_ID)) + # text is for plotly, could be GENE_ID
    geom_point(alpha = 0.7) + 
    labs(title = graph_title) + 
    geom_vline(xintercept = c(-2,2), col = "grey", linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
    geom_text_repel(max.overlaps = 10, size = 3) +  # Can do geom_text_repel or geom_label_rebel
    
    # Need it this way so the colors aren't messed up by not having significant up or down
    # scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) + 
    scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00")) # +
  
  
  # Determine the max and min axes values for labeling 
  plot_build <- ggplot_build(my_volcano)
  y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
  x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
  x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
  
  # Add the gene number annotations
  text_up <- my_df %>% filter(DE_2 == "significant up") %>% nrow()
  text_down <- my_df %>% filter(DE_2 == "significant down") %>% nrow()
  my_volcano_annotated <- my_volcano +
    annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
    annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
  
  final_volcano <- my_volcano_annotated + my_plot_themes
}


single_plot <- make_volcano_function_2(list_dfs_2[[1]], df_names[1])
single_plot
ggplotly(single_plot)

my_path <- "Figures/Volcano_plot/Log2Fold_at_2"
for (i in 1:length(list_dfs_2)) {
  
  current_df_name <- df_names[i]
  filename <- paste0(current_df_name, ".pdf")
  
  my_plot <- make_volcano_function_2(list_dfs_2[[i]], df_names[i])
  
  ggsave(my_plot,
         file = filename,
         path = my_path,
         width = 7, height = 5, units = "in")
}

###########################################################
############ FUNCTION FOR VOLCANO FDR CORRECTED ###########

make_volcano_function_FDR <- function(my_df, graph_title) {
  
  ## Make a volcano plot using output from Bob's pipeline
  
  my_volcano <- my_df %>%
    ggplot(aes(x = LOG2FOLD, y = -log10(FDR_PVALUE), col = DE2_FDR, label = DE2_FDR_labels, text = GENE_NAME, label2 = GENE_ID)) + # text is for plotly, could be GENE_ID
    geom_point(alpha = 0.7) + 
    labs(title = graph_title) + 
    geom_vline(xintercept = c(-2,2), col = "grey", linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
    geom_text_repel(max.overlaps = 10, size = 3) +  # Can do geom_text_repel or geom_label_rebel
    
    # Need it this way so the colors aren't messed up by not having significant up or down
    # scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) + 
    scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00")) # +
  
  
  # Determine the max and min axes values for labeling 
  plot_build <- ggplot_build(my_volcano)
  y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
  x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
  x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
  
  # Add the gene number annotations
  text_up <- my_df %>% filter(DE2_FDR == "significant up") %>% nrow()
  text_down <- my_df %>% filter(DE2_FDR == "significant down") %>% nrow()
  my_volcano_annotated <- my_volcano +
    annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
    annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
  
  final_volcano <- my_volcano_annotated + my_plot_themes
}


single_plot <- make_volcano_function_FDR(list_dfs_2[[1]], df_names[1])
single_plot
# ggplotly(single_plot)

# Loop for all the volcano
my_path <- "Figures/Volcano_plot/Log2Fold2_FDR"
for (i in 1:length(list_dfs_2)) {
  
  current_df_name <- df_names[i]
  filename <- paste0(current_df_name, "_FDR.pdf")
  
  my_plot <- make_volcano_function_FDR(list_dfs_2[[i]], df_names[i])
  
  ggsave(my_plot,
         file = filename,
         path = my_path,
         width = 7, height = 5, units = "in")
}




###########################################################
############ MAKE VOLCANO WITH YELLOW POINTS ##############

YellowAddition_volcano_function <- function(my_df, graph_title) {
  
  ## Make a volcano plot using output from Bob's pipeline
  
  single_gene <- my_df %>% 
    # filter(GENE_ID == "Rv0494")
    filter(GENE_ID %in% allGeneSets$`gluconeogenesis I`)
  
  my_volcano <- my_df %>%
    ggplot(aes(x = LOG2FOLD, y = -log10(AVG_PVALUE), col = DE, label = DE_labels, text = GENE_NAME, label2 = GENE_ID)) + # text is for plotly, could be GENE_ID
    geom_point() + 
    
    # Add a differently colored point
    geom_point(data = single_gene, color = "yellow", aes(col = DE, label = DE_labels, text = GENE_ID)) + 
    
    labs(title = graph_title) + 
    geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
    
    # Need it this way so the colors aren't messed up by not having significant up or down
    # scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) + 
    scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00")) # +
  
  # geom_text_repel(max.overlaps = 10, size = 3) # Can do geom_text_repel or geom_label_rebel
  
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

single_plot <- YellowAddition_volcano_function(list_dfs_2[[2]], df_names[2])
single_plot
ggplotly(single_plot)









