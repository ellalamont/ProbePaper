# Import Bob's differential gene expression data and the metagenesets data
# 7/31/25
# E. Lamont

source("Import_data.R")

# Data is coming from the Lenovo PredictTB_Run1


###########################################################
################### IMPORT BOB's DE DATA ##################

# Don't know why these are only working with the full pathname....
`GoodSputumSubset.ComparedTo.Broth` <- read.delim("Data/Differential_Expression/GoodSputumSubset_vs_Broth/W0.MTb.Meta.JOINED.txt")
`CaseumMimic.ComparedTo.Broth` <- read.delim("Data/Differential_Expression/CaseumMimic_vs_Broth/CaseumMimic.MTb.Meta.JOINED.txt")
`Rabbit.ComparedTo.Broth` <- read.delim("Data/Differential_Expression/Rabbit_vs_Broth/Rabbit.MTb.Meta.JOINED.txt")
`Marmoset.ComparedTo.Broth` <- read.delim("Data/Differential_Expression/Marmoset_vs_Broth/Marmoset.MTb.Meta.JOINED.txt")


###########################################################
################ MAKE A LIST OF ALL DFs ###################
list_dfs <- list(`GoodSputumSubset.ComparedTo.Broth`,
                 `CaseumMimic.ComparedTo.Broth`, 
                 `Rabbit.ComparedTo.Broth`, 
                 `Marmoset.ComparedTo.Broth`)

# Make a list of the names
df_names <- c("GoodSputumSubset.ComparedTo.Broth",
              "CaseumMimic.ComparedTo.Broth", 
              "Rabbit.ComparedTo.Broth", 
              "Marmoset.ComparedTo.Broth")

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
  
  # Make the column with DE gene names for plotting on graph
  current_df$DE_labels <- ifelse(current_df$DE != "not significant", current_df$GENE_NAME, NA)
  
  list_dfs_2[[current_df_name]] <- current_df
  
}
