# Import Bob's differential gene expression data and the metagenesets data
# 7/31/25
# E. Lamont

source("Import_data.R")

# Data is coming from the Lenovo PredictTB_Run1


###########################################################
################### IMPORT BOB's DE DATA ##################

`GoodSputumSubset.ComparedTo.Broth` <- read.delim("Data/Differential_Expression/GoodSputumSubset_vs_Broth/W0.MTb.Meta.JOINED.txt")
`CaseumMimic.ComparedTo.Broth` <- read.delim("Data/Differential_Expression/CaseumMimic_vs_Broth/CaseumMimic.MTb.Meta.JOINED.txt")
`Rabbit.ComparedTo.Broth` <- read.delim("Data/Differential_Expression/Rabbit_vs_Broth/Rabbit.MTb.Meta.JOINED.txt")
`Marmoset.ComparedTo.Broth` <- read.delim("Data/Differential_Expression/Marmoset_vs_Broth/Marmoset.MTb.Meta.JOINED.txt")

# Extras
`GoodSputumSubset.ComparedTo.CaseumMimic` <- read.delim("Data/Differential_Expression/GoodSputumSubset_vs_CaseumMimic/W0.MTb.Meta.JOINED.txt")
`GoodSputumSubset.ComparedTo.Rabbit` <- read.delim("Data/Differential_Expression/GoodSputumSubset_vs_Rabbit/W0.MTb.Meta.JOINED.txt")
`GoodSputumSubset.ComparedTo.Marmoset` <- read.delim("Data/Differential_Expression/GoodSputumSubset_vs_Marmoset/W0.MTb.Meta.JOINED.txt")
`CaseumMimic.ComparedTo.Rabbit` <- read.delim("Data/Differential_Expression/CaseumMimic_vs_Rabbit/CaseumMimic.MTb.Meta.JOINED.txt")
`CaseumMimic.ComparedTo.Marmoset` <- read.delim("Data/Differential_Expression/CaseumMimic_vs_Marmoset/CaseumMimic.MTb.Meta.JOINED.txt")
`Rabbit.ComparedTo.Marmoset` <- read.delim("Data/Differential_Expression/Rabbit_vs_Marmoset/Rabbit.MTb.Meta.JOINED.txt")

# Lineage split
`GoodSputumSubset_L2.ComparedTo.Broth` <- read.delim("Data/Differential_Expression/GoodSputumSubset.L2_vs_Broth/L2.MTb.Meta.JOINED.txt")
`GoodSputumSubset_L4.ComparedTo.Broth` <- read.delim("Data/Differential_Expression/GoodSputumSubset.L4_vs_Broth/L4.MTb.Meta.JOINED.txt")
`GoodSputumSubset_L4.ComparedTo.L2` <- read.delim("Data/Differential_Expression/GoodSputumSubset.L2_vs_L4/L4.MTb.Meta.JOINED.txt")

###########################################################
################ MAKE A LIST OF ALL DFs ###################
list_dfs <- list(`GoodSputumSubset.ComparedTo.Broth`,
                 `CaseumMimic.ComparedTo.Broth`, 
                 `Rabbit.ComparedTo.Broth`, 
                 `Marmoset.ComparedTo.Broth`,
                 
                 `GoodSputumSubset.ComparedTo.CaseumMimic`,
                 `GoodSputumSubset.ComparedTo.Rabbit`,
                 `GoodSputumSubset.ComparedTo.Marmoset`,
                 `CaseumMimic.ComparedTo.Rabbit`,
                 `CaseumMimic.ComparedTo.Marmoset`,
                 `Rabbit.ComparedTo.Marmoset`,
                 
                 `GoodSputumSubset_L2.ComparedTo.Broth`,
                 `GoodSputumSubset_L4.ComparedTo.Broth`,
                 `GoodSputumSubset_L4.ComparedTo.L2`)

# Make a list of the names
df_names <- c("GoodSputumSubset.ComparedTo.Broth",
              "CaseumMimic.ComparedTo.Broth", 
              "Rabbit.ComparedTo.Broth", 
              "Marmoset.ComparedTo.Broth",
              
              "GoodSputumSubset.ComparedTo.CaseumMimic",
              "GoodSputumSubset.ComparedTo.Rabbit",
              "GoodSputumSubset.ComparedTo.Marmoset",
              "CaseumMimic.ComparedTo.Rabbit",
              "CaseumMimic.ComparedTo.Marmoset",
              "Rabbit.ComparedTo.Marmoset",
              
              "GoodSputumSubset_L2.ComparedTo.Broth",
              "GoodSputumSubset_L4.ComparedTo.Broth",
              "GoodSputumSubset_L4.ComparedTo.L2")

# Give the df list the correct df names
names(list_dfs) <- df_names



###########################################################
############### ADD COLUMNS OF DE VALUES ##################

# Make a new list to hold dataframes with extra columns
list_dfs_2 <- list()

ordered_DE <- c("significant down", "not significant", "significant up")

# Add extra DE columns to each dataframe
for (i in 1:length(list_dfs)) {
  
  current_df <- list_dfs[[i]]
  current_df_name <- df_names[i]
  
  # Make the column pointing out which ones are differentially expressed
  current_df$DE <- ifelse(current_df$LOG2FOLD < -1 & current_df$AVG_PVALUE < 0.05, "significant down",
                          ifelse(current_df$LOG2FOLD > 1 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE <- factor(current_df$DE, levels = ordered_DE)
  
  # Make another column where the threshold for DEG is 2 (above is 1)
  current_df$DE_2 <- ifelse(current_df$LOG2FOLD < -2 & current_df$AVG_PVALUE < 0.05, "significant down",
                          ifelse(current_df$LOG2FOLD > 2 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE_2 <- factor(current_df$DE_2, levels = ordered_DE)
  
  # Make the column with DE gene names for plotting on graph
  current_df$DE_labels <- ifelse(current_df$DE != "not significant", current_df$GENE_NAME, NA)
  current_df$DE_2_labels <- ifelse(current_df$DE_2 != "not significant", current_df$GENE_NAME, NA)
  
  # Columns for Log2Fold>abs(2) and FDR corrected p-values
  current_df$FDR_PVALUE <- p.adjust(current_df$AVG_PVALUE, method = "fdr")
  current_df$DE2_FDR <- ifelse(current_df$LOG2FOLD < -2 & current_df$FDR_PVALUE < 0.05, "significant down",
                               ifelse(current_df$LOG2FOLD > 2 & current_df$FDR_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE2_FDR <- factor(current_df$DE2_FDR, levels = ordered_DE)
  current_df$DE2_FDR_labels <- ifelse(current_df$DE2_FDR != "not significant", current_df$GENE_NAME, NA)
  
  
  list_dfs_2[[current_df_name]] <- current_df
}

###########################################################
################# REMOVE NON-CODING GENES #################
# 8/15/25: After talking to DRS, decided to remove all the non-coding genes and all the MT genes, and leave just the coding genes starting with Rv. So need to remove these at the raw read level and calculate new TPM
# The Pathcap people also had issues with ncRNAs: https://www.nature.com/articles/s41598-019-55633-6#Sec8

list_dfs_f <- lapply(list_dfs_2, function(df) {
  df %>% filter(str_detect(GENE_ID, "^Rv\\d+.*"))
})



###########################################################
############# IMPORT BOB's METAGENESETS DATA ##############


`MetaGeneSets_GoodSputumSubset.vs.Broth_UP` <- read.delim("Data/Differential_Expression/GoodSputumSubset_vs_Broth/W0.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_Marmoset.vs.Broth_UP` <- read.delim("Data/Differential_Expression/Marmoset_vs_Broth/Marmoset.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_Rabbit.vs.Broth_UP` <- read.delim("Data/Differential_Expression/Rabbit_vs_Broth/Rabbit.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_CaseumMimic.vs.Broth_UP` <- read.delim("Data/Differential_Expression/CaseumMimic_vs_Broth/CaseumMimic.MTb.MetaGeneSets.UP.txt")

`MetaGeneSets_GoodSputumSubset_L2.vs.Broth_UP` <- read.delim("Data/Differential_Expression/GoodSputumSubset.L2_vs_Broth/L2.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_GoodSputumSubset_L4.vs.Broth_UP` <- read.delim("Data/Differential_Expression/GoodSputumSubset.L4_vs_Broth/L4.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_GoodSputumSubset_L4.vs.L2_UP` <- read.delim("Data/Differential_Expression/GoodSputumSubset.L2_vs_L4/L4.MTb.MetaGeneSets.UP.txt")

