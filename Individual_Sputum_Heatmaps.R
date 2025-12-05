# Try heatmaps of all (some) of the genesets
# 12/5/25
# In response to reviewer comments

# Using Pheatmap
library(pheatmap)

### pheatmap save function (https://gist.github.com/mathzero/a2070a24a6b418740c44a5c023f5c01e)
# save_pheatmap <- function(x, filename, width=12, height=12){
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   if(grepl(".png",filename)){
#     png(filename, width=width, height=height, units = "in", res=300)
#     grid::grid.newpage()
#     grid::grid.draw(x$gtable)
#     dev.off()
#   }
#   else if(grepl(".pdf",filename)){
#     pdf(filename, width=width, height=height)
#     grid::grid.newpage()
#     grid::grid.draw(x$gtable)
#     dev.off()
#   }
#   else{
#     print("Filename did not contain '.png' or '.pdf'")
#   }
# }

source("Import_data.R")
source("Import_GeneSets.R")
source("Import_DEG_sets.R")


###########################################################
##################### ORGANIZE THE DATA ###################

# Just want my 9 sputum samples and the Rv samples (which will have to be averaged)
my_tpm <- All_tpm_f %>% select(all_of(Sputum_GoodSampleList), contains("H37Ra"))

# Grab the columns needed to give colors
my_pipeSummary <- BiolSamples_pipeSummary %>%
  filter(SampleID %in% colnames(my_tpm)) %>%
  select(SampleID, Type3) %>%
  column_to_rownames("SampleID")

# Define the colors
my_annotation_colors <- list(
  Type3 = c("Sputum L2" = "#99CCE8",
           "Sputum L4" = "#004C73",
           "Broth" = "#999999"))

###########################################################
######################## PHEATMAP #########################

# Rv0681
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList[["MTb.iModulons"]][["Rv0681: Transcription with previously unknown function, regulates core fatty acid response "]])
heatmap_plot <- pheatmap(myData,
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         annotation_col = my_pipeSummary, 
         annotation_colors = my_annotation_colors,
         scale = "row",
         # display_numbers = T,
         annotation_names_col = F,
         cluster_rows = F,
         cutree_cols = 1,
         main = "Rv0681",
         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
         treeheight_col = 10,
         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/Rv0681.png", width = 5, height = 6)

# PrpR
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`PrpR: Transcription factor that is involved in catabolism of short chain fatty acids via gloxylate and methylcitrate cycle `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "PrpR",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/PrpR.png", width = 5, height = 2.8)

# Peptidoglycan Biosynthesis
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`Peptidoglycan Biosynthesis: iModulon enriched with the Peptidoglycan biosyntehsis KEGG pathway `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "Peptidoglycan Biosynthesis",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/PeptidoglycanBiosynthesis.png", width = 5, height = 3.2)

# Mycofactocin Synthesis Pathway
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`Mycofactocin Synthesis Pathway: Manually annotated iModulon that captures multiple genes associated with mycofactocin metabolism `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "Mycofactocin Synthesis Pathway",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/MycofactocinSynthesisPathway.png", width = 5, height = 2.9)

# KstR2
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`KstR2: Transcription factor that regulates a small regulon related to cholesterol utilization `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "KstR2",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/KstR2.png", width = 5, height = 3.8)

# Fatty Acid Biosynthesis
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`Fatty Acid Biosynthesis: iModulon that maps to the intersection of the Rv1033c, Rv1776c, and Rv3681c regulons `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "Fatty Acid Biosynthesis",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/FattyAcidBiosynthesis.png", width = 5, height = 4.4)


# FasR
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`FasR: Transcription factor that controls the expression of fatty acid metabolism genes `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "FasR",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/FasR.png", width = 5, height = 3.7)


# BkaR
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`BkaR: Transcription factor that regulates genes related to branched-chain keto-acid metabolism `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "BkaR",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/BkaR.png", width = 5, height = 3.2)

# WhiB1
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`WhiB1: Transcription repressor that is redox responsive, found to be nitric oxide sensitive `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "WhiB1",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/WhiB1.png", width = 5, height = 6)

# Rv1828/SigH
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`Rv1828/SigH: iModulon that maps to the union of the Rv1828 and Rv3223c regulons `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "Rv1828/SigH",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/Rv1828SigH.png", width = 5, height = 8)

# Positive Regulation of Growth
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`Positive Regulation of Growth: iModulon enriched with the Positive Regulation of Growth Gene Ontology term `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "Positive Regulation of Growth",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/PositiveRegulationofGrowth.png", width = 5, height = 6)

# GroEL-GroES Complex
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`GroEL-GroES Complex: Gene Ontology Term specific to the bacterial chaperonin complex `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "GroEL-GroES Complex",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/GroELGroESComplex.png", width = 5, height = 3.2)

# DevR-2
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`DevR-2: Transcription factor that is one of the primary regulators of hypoxia onset response `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "DevR-2",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/DevR2.png", width = 5, height = 6.3)


# DevR-1
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`DevR-1: Transcription factor that is one of the primary regulators of hypoxia onset response `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "DevR-1",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/DevR1.png", width = 5, height = 4.5)

# Zur
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`Zur: Transcription factor that is activated by and regulates the metabolism of zinc `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "Zur",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/Zur.png", width = 5, height = 4.3)

# RicR
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`RicR: Copper responsive transcription factor that regulates metal metabolism genes `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "RicR",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/RicR.png", width = 5, height = 3.8)

# IdeR
myData <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList$MTb.iModulons$`IdeR: Metal-Dependent DNA-Binding protein that controls genes related to iron metabolism `)
heatmap_plot <- pheatmap(myData, 
                         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         annotation_col = my_pipeSummary, 
                         annotation_colors = my_annotation_colors,
                         scale = "row",
                         # display_numbers = T,
                         annotation_names_col = F,
                         cluster_rows = F,
                         cutree_cols = 1,
                         main = "IdeR",
                         fontsize = 9, fontsize_col = 9, fontsize_row = 9,
                         treeheight_col = 10,
                         border_color = "#3b3b3b") 
save_pheatmap(heatmap_plot, "Figures/Pheatmap/IdeR.png", width = 5, height = 7)






