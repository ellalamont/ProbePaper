# Try to make my own forest plot in ggplot, using the data that Bob's forest plots generate
# E. Lamont 
# 12/05/25

# Adjusted to just be sputum and updated for all the other stuff I have done

# ggforestplot package looks like the best

# Lets just focus on iModulons again

# This would be easier if I had the underlying math going on in the forest plots

source("Import_DEG_sets.R")
source("Import_GeneSets.R")

install.packages("tidytext")
library(tidytext)

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

Sputum_data <- plotGeneSetForest(file = list_dfs_f2$GoodSputumSubset.ComparedTo.Broth,
                                 geneSets = allGeneSetList$MTb.iModulons,
                                 main = "H37Ra vs W0 sputum subset",
                                 min.genes.per.set = 3,
                                 max.show = 80, # There will only be 71 because of the 3 threshold above
                                 text.cex = 1.1, pt.cex = 1.25, lwd = 3.5)
write.csv(Sputum_data,
          file = "Figures/ForestPlots/EllaNew/Sputum_iModulonsForest.csv")


# Combine all the dataframes
combined_df <- bind_rows(
  Sputum_data %>% mutate(Type = "Sputum"),
  Marmoset_data %>% mutate(Type = "Marmoset"),
  Rabbit_data %>% mutate(Type = "Rabbit"),
  Mimic_data %>% mutate(Type = "Caseum mimic")
) %>% select(-Unused) %>% 
  mutate(Label = str_wrap(Label, width = 50))

# Split the N_mean_SD column 
Sputum_data_2 <- Sputum_data %>%
  separate(N_Mean_SD, into = c("N", "mean_sd"), sep = ", ") %>%
  mutate(Mean = str_extract(mean_sd, "-?\\d+\\.\\d+") %>% as.numeric(),
         SD = str_extract(mean_sd, "\\(.*?\\)") %>% str_remove_all("[()]") %>% as.numeric) %>% 
  select(-mean_sd)

# Get the single CI value for error bars
Sputum_data_3 <- Sputum_data_2 %>%
  mutate(
    CI_clean = str_remove_all(CI, "[()]"),
    CI_lower = as.numeric(str_extract(CI_clean, "^-?\\d*\\.?\\d+")),
    CI_upper = as.numeric(str_extract(CI_clean, "(?<=, )-?\\d*\\.?\\d+")),
    CI2 = abs(Mean - CI_lower)
  )


# Add the iModulon category
Sputum_data_4 <- Sputum_data_3 %>%
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
Sputum_data_4 <- Sputum_data_4 %>%
  mutate(Label_2 = paste0(Label, " (n=", N, ")"))

###########################################################
################### GGFORESTPLOT PACKAGE ##################

# https://nightingalehealth.github.io/ggforestplot/articles/ggforestplot.html

test <- Sputum_data_4 %>% arrange(Label)

ggforestplot::forestplot(
  df = test[1:30,],
  name = Label,
  # colour = Type,
  estimate = Mean,
  se = CI2
)

test_graph <- test[1:30,] %>% 
  ggforestplot::forestplot(name = Label, estimate = Mean, pvalue = P_Value, 
             se = CI2/1.96) + # Need to divide by 1.96 because I am giving it the true CI, not the SE.
  # scale_color_manual(values = my_fav_colors) + 
  labs(x = "LOG2FOLD from broth") + 
  my_plot_themes 
# test_graph$layers[[1]]$aes_params$odd <- "#00000000" # Run this code to remove the stripes (https://stackoverflow.com/questions/71745719/how-to-control-stripe-transparency-using-ggforestplot-geom-stripes)
test_graph



###########################################################
############# LOOP FOR ALL iMODULON CATEGORIES ############

# ForestPlot_Function <- function(df, title_text = NULL){
#   
#   forestplot(df, name = Label_2, colour = Type, estimate = Mean, pvalue = P_Value, 
#              se = CI2/1.96) + # Need to divide by 1.96 because I am giving it the true CI, not the SE.
#     scale_color_manual(values = my_fav_colors) + 
#     labs(title = title_text,
#          subtitle = "mean, error bars are 95% CI",
#          x = "LOG2FOLD from broth") + 
#     my_plot_themes
# }
# 
# # Run just one with the function
# # ForestPlot_Function(combined_df_4 %>% filter(iModulonCategory == "Redox"))
# 
# for (category in unique(na.omit(combined_df_4$iModulonCategory))) {
#   df_cat <- combined_df_4 %>% filter(iModulonCategory == category)
#   
#   p <- ForestPlot_Function(df_cat, title_text = paste(category, " iModulon"))
#   
#   ggsave(p,
#          file = paste0(category,"_Forest.pdf"),
#          path = "Figures/ForestPlots/EllaNew/iModulons",
#          width = 7.5, height = 6, units = "in")
# }


###########################################################
###################### FACET BY GROUP #####################
# https://nightingalehealth.github.io/ggforestplot/articles/ggforestplot.html


Sputum_data_5 <- Sputum_data_4 %>%
  mutate(Label_2 = str_wrap(Label_2, width = 50))

iModulons_of_interest <- c("PrpR", "BkaR", "Rv0681", "KstR2", "Fatty Acid Biosynthesis", "FasR", "Peptidoglycan Biosynthesis", "Positive Regulation of Growth", "Mycofactocin Synthesis Pathway", "DevR-2", "GroEL-GroES Complex", "DevR-1", "PhoP", "MprA", "Mce3R", "Mce1R", "SigC", "SigD", "SigH", "SigK", "RicR", "IdeR", "Zur", "WhiB1") 
iModulons_of_interest_pattern <- str_c(iModulons_of_interest, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

Fav_Pathways <- Sputum_data_5 %>% 
  filter(str_detect(Label_2, iModulons_of_interest_pattern)) %>% 
  mutate(Label_2 = str_wrap(Label_2, width = 50)) %>%
  pull(Label_2)

# Put them in categories
Sputum_data_6 <- Sputum_data_5 %>%
  mutate(iModulonCategory2 = case_when(
    str_detect(Label_2, "DevR") ~ "DosR", # specifically putt the DevR in a different category (normally in redox)
    str_detect(Label_2, paste(Growth_iModulons_pattern, NucleicAcid_iModulons_pattern, Redox_iModulons_pattern, AminoAcid_iModulons_pattern, sep = "|")) ~ "Growth",
    str_detect(Label_2, Metal_iModulons_pattern) ~ "Metal",
    str_detect(Label_2, Virulence.Persistence_iModulons_pattern) ~ "Virulence and Persistence",
    str_detect(Label_2, paste(CentralCarbon_iModulons_pattern, FattyAcid.Cholesterol_iModulons_pattern, sep = "|")) ~ "Fatty Acid and Cholesterol",
    TRUE ~ "Other"
  )) %>%
  mutate(iModulonCategory2 = factor(iModulonCategory2, levels = c("Fatty Acid and Cholesterol", "Growth", "DosR", "Metal","Virulence and Persistence", "Other"))) %>%
  filter(Label_2 %in% Fav_Pathways) %>%
  filter(iModulonCategory2 != "Virulence and Persistence") %>% 
  filter(!str_detect(Label_2, "WhiB4/IdeR")) # To remove this iModulon which is tagging along with IdeR iModulon

# Re-order to match the bubble plot... not sure why it is different
my_labels <- c("Rv0681: Transcription with previously unknown\nfunction, regulates core fatty acid response\n(n=27)", "PrpR: Transcription factor that is involved\nin catabolism of short chain fatty acids via\ngloxylate and methylcitrate cycle (n=5)", "Peptidoglycan Biosynthesis: iModulon enriched with\nthe Peptidoglycan biosyntehsis KEGG pathway (n=8)", "Mycofactocin Synthesis Pathway: Manually annotated\niModulon that captures multiple genes associated\nwith mycofactocin metabolism (n=6)", "KstR2: Transcription factor that regulates a small\nregulon related to cholesterol utilization (n=14)", "Fatty Acid Biosynthesis: iModulon that maps to the\nintersection of the Rv1033c, Rv1776c, and Rv3681c\nregulons (n=21)", "FasR: Transcription factor that controls the\nexpression of fatty acid metabolism genes (n=12)", "BkaR: Transcription factor that regulates genes\nrelated to branched-chain keto-acid metabolism\n(n=8)", "WhiB1: Transcription repressor that is redox\nresponsive, found to be nitric oxide sensitive\n(n=32)", "Rv1828/SigH: iModulon that maps to the union of\nthe Rv1828 and Rv3223c regulons (n=61)", "Positive Regulation of Growth: iModulon enriched\nwith the Positive Regulation of Growth Gene\nOntology term (n=26)", "GroEL-GroES Complex: Gene Ontology Term specific\nto the bacterial chaperonin complex (n=8)", "DevR-2: Transcription factor that is one of the\nprimary regulators of hypoxia onset response\n(n=36)", "DevR-1: Transcription factor that is one of the\nprimary regulators of hypoxia onset response\n(n=22)", "Zur: Transcription factor that is activated by and\nregulates the metabolism of zinc (n=20)", "RicR: Copper responsive transcription factor that\nregulates metal metabolism genes (n=14)", "IdeR: Metal-Dependent DNA-Binding protein that\ncontrols genes related to iron metabolism (n=44)")
Sputum_data_6 <- Sputum_data_6 %>%
  mutate(Label_2 = factor(Label_2, levels = my_labels))

# STILL NOT GETTING THE ORDER RIGHT!

ForestPlot1 <- Sputum_data_6 %>% 
  forestplot(name = Label_2, estimate = Mean, pvalue = P_Value, 
                           se = CI2/1.96) + # Need to divide by 1.96 because I am giving it the true CI, not the SE.
  # scale_color_manual(values = my_fav_colors) + 
  labs(x = "LOG2FOLD change") + 
  # geom_text(aes(x = -3.55, label = paste0("n = ", N)), hjust = 1, size = 2.7) +
  facet_grid(rows = vars(iModulonCategory2), scales = "free_y", space = "free") + 
  my_plot_themes + facet_themes
ForestPlot1
ggsave(ForestPlot1,
       file = "ForestPlot_1.pdf",
       path = "Figures/ForestPlots/EllaNew/iModulons",
       width = 7.5, height = 6, units = "in")






