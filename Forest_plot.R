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

source("Import_DEG_sets.R")


###########################################################
###################### iMODULONS ##########################

plotGeneSetForest(file = list_dfs_2$GoodSputumSubset.ComparedTo.Broth,
                  geneSets = allGeneSetList$MTb.iModulons,
                  main = "iModulons Forest Plot Good Sputum subset vs Broth",
                  left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=11)",
                  max.show = 30,
                  xRange = 4, # Changes how far out the log2fold change axis goes
                  text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 

# How to save the data:
# Save the plot
printPlot(filename = "Figures/ForestPlots/GoodSputumSamples_iModulons_v1.pdf", width = 18, height = 10)


# Save the .csv
# Forestplot_Ella <- plotGeneSetForest(file = list_dfs$W0.ComparedTo.Broth,
#                                      geneSets = allGeneSetList$Ella_GeneSets,
#                                      main = "Ella Gene Sets Forest Plot H37Ra vs W0 sputum",
#                                      left.label = "H37Ra broth (n=3)", right.label = "W0 sputum (n=3)",
#                                      xRange = 4, # Changes how far out the log2fold change axis goes
#                                      text.cex = 1.1, pt.cex = 1.25, lwd = 3.5) 
# class(Forestplot_Ella)
# write.csv(Forestplot_Ella,
#           file = "ForestPlot_Figures/ForestPlot_EllaGeneSets_data.csv")


