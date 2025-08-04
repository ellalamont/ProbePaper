# Try to make my own forest plot in ggplot, using the data that Bob's forest plots generate
# E. Lamont 
# 8/2/25

# ggforestplot package looks like the best

# Lets just focus on iModulons again

# This would be easier if I had the underlying math going on in the forest plots

source("Import_DEG_sets.R")
source("Import_GeneSets.R")

# Plot basics
my_plot_themes <- theme_bw() +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = "white", color = NA)) + 
  theme(legend.position = "right",legend.text=element_text(size=9),
        legend.title = element_text(size = 10),
        # legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=9, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=8), 
        plot.subtitle = element_text(size=10), 
        plot.margin = margin(10, 10, 10, 20)# ,
  )

facet_themes <- theme(strip.background=element_rect(fill="white", linewidth = 0.9),
                      strip.text = element_text(size = 10))


###########################################################
########### SAVE DATA FROM BOB'S FOREST PLOTS #############

Sputum_data <- plotGeneSetForest(file = list_dfs_2$GoodSputumSubset.ComparedTo.Broth,
                             geneSets = allGeneSetList$MTb.iModulons,
                             main = "H37Ra vs W0 sputum subset",
                             min.genes.per.set = 3,
                             max.show = 80,
                             text.cex = 1.1, pt.cex = 1.25, lwd = 3.5)
write.csv(Sputum_data,
          file = "Figures/ForestPlots/EllaNew/Sputum_iModulonsForest.csv")

Marmoset_data <- plotGeneSetForest(file = list_dfs_2$Marmoset.ComparedTo.Broth,
                                 geneSets = allGeneSetList$MTb.iModulons,
                                 main = "H37Ra vs W0 sputum subset",
                                 min.genes.per.set = 3,
                                 max.show = 80,
                                 text.cex = 1.1, pt.cex = 1.25, lwd = 3.5)

Rabbit_data <- plotGeneSetForest(file = list_dfs_2$Rabbit.ComparedTo.Broth,
                                   geneSets = allGeneSetList$MTb.iModulons,
                                   main = "H37Ra vs W0 sputum subset",
                                   min.genes.per.set = 3,
                                   max.show = 80,
                                   text.cex = 1.1, pt.cex = 1.25, lwd = 3.5)
Mimic_data <- plotGeneSetForest(file = list_dfs_2$CaseumMimic.ComparedTo.Broth,
                                 geneSets = allGeneSetList$MTb.iModulons,
                                 main = "H37Ra vs W0 sputum subset",
                                 min.genes.per.set = 3,
                                 max.show = 80,
                                 text.cex = 1.1, pt.cex = 1.25, lwd = 3.5)

# Combine all the dataframes
combined_df <- bind_rows(
  Sputum_data %>% mutate(Type = "Sputum"),
  Marmoset_data %>% mutate(Type = "Marmoset"),
  Rabbit_data %>% mutate(Type = "Rabbit"),
  Mimic_data %>% mutate(Type = "Caseum mimic")
) %>% select(-Unused) %>% 
  mutate(Label = str_wrap(Label, width = 50))

# Split the N_mean_SD column 
combined_df_2 <- combined_df %>%
  separate(N_Mean_SD, into = c("N", "mean_sd"), sep = ", ") %>%
  mutate(Mean = str_extract(mean_sd, "-?\\d+\\.\\d+") %>% as.numeric(),
    SD = str_extract(mean_sd, "\\(.*?\\)") %>% str_remove_all("[()]") %>% as.numeric) %>% 
  select(-mean_sd)

# Get the single CI value for error bars
combined_df_3 <- combined_df_2 %>%
  mutate(
    CI_clean = str_remove_all(CI, "[()]"),
    CI_lower = as.numeric(str_extract(CI_clean, "^-?\\d*\\.?\\d+")),
    CI_upper = as.numeric(str_extract(CI_clean, "(?<=, )-?\\d*\\.?\\d+")),
    CI2 = abs(Mean - CI_lower)
  )


# Add the iModulon category
combined_df_4 <- combined_df_3 %>%
  mutate(iModulonCategory = case_when(
    str_detect(Label, CentralCarbon_iModulons_pattern) ~ "Central Carbon",
    str_detect(Label, AminoAcid_iModulons_pattern) ~ "Amino Acid",
    str_detect(Label, NucleicAcid_iModulons_pattern) ~ "Nucleic Acid",
    str_detect(Label, FattyAcid.Cholesterol_iModulons_pattern) ~ "Fatty Acid_Cholesterol",
    str_detect(Label, Metal_iModulons_pattern) ~ "Metal",
    str_detect(Label, SulfurMetabolism_iModulons_pattern) ~ "Sulfur Metabolism",
    str_detect(Label, Growth_iModulons_pattern) ~ "Growth",
    str_detect(Label, Redox_iModulons_pattern) ~ "Redox",
    str_detect(Label, AcidStress_iModulons_pattern) ~ "Acid Stress",
    str_detect(Label, Antibiotic_iModulons_pattern) ~ "Antibiotic",
    str_detect(Label, Virulence.Persistence_iModulons_pattern) ~ "Virulence_Persistence",
    TRUE ~ "Other"
  ))

# Make a new Label category with size
combined_df_4 <- combined_df_4 %>%
  mutate(Label_2 = paste0(Label, " (n=", N, ")"))

###########################################################
################### GGFORESTPLOT PACKAGE ##################

# https://nightingalehealth.github.io/ggforestplot/articles/ggforestplot.html

test <- combined_df_4 %>% arrange(Label)

ggforestplot::forestplot(
  df = test[1:30,],
  name = Label,
  colour = Type,
  estimate = Mean,
  se = CI2
)

test_graph <- test[1:30,] %>% 
  forestplot(name = Label, colour = Type, estimate = Mean, pvalue = P_Value, 
             se = CI2/1.96) + # Need to divide by 1.96 because I am giving it the true CI, not the SE.
  scale_color_manual(values = my_fav_colors) + 
  labs(x = "LOG2FOLD from broth") + 
  my_plot_themes 
# test_graph$layers[[1]]$aes_params$odd <- "#00000000" # Run this code to remove the stripes (https://stackoverflow.com/questions/71745719/how-to-control-stripe-transparency-using-ggforestplot-geom-stripes)
test_graph



###########################################################
############# LOOP FOR ALL iMODULON CATEGORIES ############

ForestPlot_Function <- function(df, title_text = NULL){
  
  forestplot(df, name = Label_2, colour = Type, estimate = Mean, pvalue = P_Value, 
             se = CI2/1.96) + # Need to divide by 1.96 because I am giving it the true CI, not the SE.
    scale_color_manual(values = my_fav_colors) + 
    labs(title = title_text,
         subtitle = "mean, error bars are 95% CI",
         x = "LOG2FOLD from broth") + 
    my_plot_themes
}

# Run just one with the function
# ForestPlot_Function(combined_df_4 %>% filter(iModulonCategory == "Redox"))

for (category in unique(na.omit(combined_df_4$iModulonCategory))) {
  df_cat <- combined_df_4 %>% filter(iModulonCategory == category)
  
  p <- ForestPlot_Function(df_cat, title_text = paste(category, " iModulon"))
  
  ggsave(p,
         file = paste0(category,"_Forest.pdf"),
         path = "Figures/ForestPlots/EllaNew/iModulons",
         width = 7.5, height = 6, units = "in")
}


###########################################################
###################### FACET BY GROUP #####################
# https://nightingalehealth.github.io/ggforestplot/articles/ggforestplot.html

# Choose a random smaller subset
test_df <- combined_df_4 %>% filter(iModulonCategory %in% c("Acid Stress", "Fatty Acid_Cholesterol", "Metal"))

test_graph <- test_df %>% 
  forestplot(name = Label_2, colour = Type, estimate = Mean, pvalue = P_Value, 
             se = CI2/1.96) + # Need to divide by 1.96 because I am giving it the true CI, not the SE.
  scale_color_manual(values = my_fav_colors) + 
  labs(x = "LOG2FOLD from broth") + 
  # geom_text(aes(x = -3.55, label = paste0("n = ", N)), hjust = 1, size = 2.7) +
  facet_col(facets = ~iModulonCategory, scales = "free_y", space = "free") + 
  my_plot_themes + facet_themes
test_graph
ggsave(test_graph,
       file = "Groups_testing.pdf",
       path = "Figures/ForestPlots/EllaNew/iModulons",
       width = 7.5, height = 6, units = "in")


