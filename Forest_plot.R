# Forest plot from Bob
# 7/31/25
# E. Lamont

# From Bob:
# Main arguments are the first 2:
#   1) file - either the text string of the full pathname to a DE results file, such as "results/MetaResults/MTb.Test/W0.Mtb.Meta.JOINED.txt".   Or a data frame of DE results
# 2) geneSets - either a GeneSet as an R list object (each list element a vector of geneIDs), or a text string of the name of a predefined gene list.  Such as "MTb.Phenotypes"
# Optional arguments control various behaviors.  the "left.label" and "right.label" are the text we paste below the forest.  "xRange" gives more control over the X axis portion of the plot. Etc.
# Lastly, there is a Duffy function called "printPlot()" that turns whatever is in the current graphics window into a PDF file, with arguments to control plot shape.  Most pipeline tools use that to make the various PDF files the workflow creates.
# Totally forgot to mention-  the plotGeneSetForest function returns an (invisible) data frame of those plotted facts.  So when you call it, assign the result to an R object.  As in:
#   forestAns <- plotGeneSetForest( file, geneSets, ....)
#   print( forestAns)


###########################################################
###################### iMODULONS ##########################

plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
                  geneSets = allGeneSetList$Ella_GeneSets,
                  main = "Ella Gene Sets Forest Plot H37Ra vs W0 sputum",
                  left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=3)",
                  xRange = 4, # Changes how far out the log2fold change axis goes
                  text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 

# How to save the data:
# Save the plot
printPlot(filename = "ForestPlot_Figures/ForestPlot_EllaGeneSets_v1.pdf", width = 12, height = 8)


# Save the .csv
Forestplot_Ella <- plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
                                     geneSets = allGeneSetList$Ella_GeneSets,
                                     main = "Ella Gene Sets Forest Plot H37Ra vs W0 sputum",
                                     left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=3)",
                                     xRange = 4, # Changes how far out the log2fold change axis goes
                                     text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 
class(Forestplot_Ella)
write.csv(Forestplot_Ella,
          file = "ForestPlot_Figures/ForestPlot_EllaGeneSets_data.csv")


###########################################################
############ WALTER2015 GENE SETS FOREST PLOT #############
plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
                  geneSets = allGeneSetList$Walter2015GeneSets,
                  main = "Walter2015 Gene Sets Forest Plot H37Ra vs W0 sputum",
                  left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=3)",
                  xRange = 4, # Changes how far out the log2fold change axis goes
                  text.cex = 1, pt.cex = 1.25, lwd = 3.5) 


###########################################################
################## CLEARTB FOREST PLOT ####################

allGeneSetList$ClearTB_GeneSetList

colors <- c(
  rep("brown", 2),
  rep("#17becf", 4),
  rep("green4", 2)
)

plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
                  geneSets = allGeneSetList$ClearTB_GeneSetList,
                  main = "Forest Plot H37Ra vs W0 sputum",
                  # left.label = "Mtb broth (n=3)", right.label = "W0 sputum (n=3)",
                  xRange = 2.5, # Changes how far out the log2fold change axis goes
                  text.cex = 1.15, pt.cex = 1.25, lwd = 3.5,
                  min.genes.per.set = 5,
                  col = colors) 

# Save the plot
printPlot(filename = "ForestPlot_Figures/ForestPlot_ClearTB_v3.pdf", width = 12, height = 8)


###########################################################
################# iModulons FOREST PLOT ###################

plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
                  geneSets = allGeneSetList$MTb.iModulons,
                  # geneSets = allGeneSetList$MTb.iModulons[str_detect(names(allGeneSetList$MTb.iModulons), pattern)], # Messing around here
                  max.show = 40,
                  min.genes.per.set = 3, 
                  main = "iModulons Forest Plot H37Ra vs W0 sputum",
                  left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=3)",
                  xRange = 4, # Changes how far out the log2fold change axis goes
                  text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 
# Save the plot
printPlot(filename = "ForestPlot_Figures/ForestPlot_iModulons_Subset40.pdf", width = 18, height = 10)
printPlot(filename = "ForestPlot_Figures/ForestPlot_iModulons_Subset77.pdf", width = 18, height = 18)

# Save the .csv
Forestplot_Ella <- plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
                                     geneSets = allGeneSetList$MTb.iModulons,
                                     max.show = 80,
                                     main = "iModulons Forest Plot H37Ra vs W0 sputum",
                                     left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=3)",
                                     xRange = 4, # Changes how far out the log2fold change axis goes
                                     text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 
write.csv(Forestplot_Ella,
          file = "ForestPlot_Figures/ForestPlot_iModulons_data.csv")


###########################################################
################## Regulons FOREST PLOT ###################

plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
                  geneSets = allGeneSetList$MTb.Regulons,
                  max.show = 40,
                  min.genes.per.set = 3, 
                  main = "Regulons Forest Plot H37Ra vs W0 sputum",
                  left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=3)",
                  xRange = 4, # Changes how far out the log2fold change axis goes
                  text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 
# Save the plot
printPlot(filename = "ForestPlot_Figures/ForestPlot_Regulons_Subset40.pdf", width = 18, height = 10)


###########################################################
############### TFOE.Regulons FOREST PLOT #################

plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
                  geneSets = allGeneSetList$MTb.TFOE.Regulons,
                  max.show = 40,
                  min.genes.per.set = 3, 
                  main = "TFOE.Regulons Forest Plot H37Ra vs W0 sputum",
                  left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=3)",
                  xRange = 4, # Changes how far out the log2fold change axis goes
                  text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 
# Save the plot
printPlot(filename = "ForestPlot_Figures/ForestPlot_TFOE.Regulons_Subset40.pdf", width = 18, height = 10)


###########################################################
########## Tuberculist.GO.Ontology FOREST PLOT ############

plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
                  geneSets = allGeneSetList$MTb.Tuberculist.GO.Ontology,
                  max.show = 40,
                  min.genes.per.set = 3, 
                  main = "Tuberculist.GO.Ontology Forest Plot H37Ra vs W0 sputum",
                  left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=3)",
                  xRange = 4, # Changes how far out the log2fold change axis goes
                  text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 
# Save the plot
printPlot(filename = "ForestPlot_Figures/ForestPlot_Tuberculist.GO.Ontology_Subset40.pdf", width = 18, height = 10)


###########################################################
############# CellWallSynthesis FOREST PLOT ###############

plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
                  geneSets = allGeneSetList$CellWallSynthesis_GeneSets,
                  max.show = 40,
                  min.genes.per.set = 3, 
                  main = "Cell Wall Synthesis Forest Plot H37Ra vs W0 sputum",
                  left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=3)",
                  xRange = 4, # Changes how far out the log2fold change axis goes
                  text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 
# Save the plot
printPlot(filename = "ForestPlot_Figures/ForestPlot_W0vsBroth_CellWallSynthesis.pdf", width = 16, height = 8)

plotGeneSetForest(file = list_dfs$W2.ComparedTo.Broth,
                  geneSets = allGeneSetList$CellWallSynthesis_GeneSets,
                  max.show = 40,
                  min.genes.per.set = 3, 
                  main = "Cell Wall Synthesis Forest Plot H37Ra vs W2 sputum",
                  left.label = "H37Ra broth (n=3)", right.label = "W2 sputum (n=3)",
                  xRange = 4, # Changes how far out the log2fold change axis goes
                  text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 
# Save the plot
printPlot(filename = "ForestPlot_Figures/ForestPlot_W2vsBroth_CellWallSynthesis.pdf", width = 16, height = 8)

