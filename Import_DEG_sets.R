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

# With Indigo Rv or Lance Rv
`Indigo_Rv.ComparedTo.Ra` <- read.delim("Data/Differential_Expression/Ra_vs_Rv.Indigo/Rv.MTb.Meta.JOINED.txt")
`GoodSputumSubset.ComparedTo.Indigo_Rv` <- read.delim("Data/Differential_Expression/SputumSubset_vs_Rv.Indigo/W0.MTb.Meta.JOINED.txt")
`LancepH7_Rv.ComparedTo.Ra` <- read.delim("Data/Differential_Expression/Ra_vs_Rv.LancepH7/Rv_pH7.MTb.Meta.JOINED.txt")
`GoodSputumSubset.ComparedTo.LancepH7_Rv` <- read.delim("Data/Differential_Expression/SputumSubset_vs_Rv.LancepH7/W0.MTb.Meta.JOINED.txt")
`LancepH7_Rv.ComparedTo.Indigo_Rv` <- read.delim("Data/Differential_Expression/LancepH7.Rv_vs_Indigo_Rv/Rv_pH7.MTb.Meta.JOINED.txt")

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
                 `GoodSputumSubset_L4.ComparedTo.L2`,
                 
                 `Indigo_Rv.ComparedTo.Ra`,
                 `GoodSputumSubset.ComparedTo.Indigo_Rv`,
                 `LancepH7_Rv.ComparedTo.Ra`,
                 `GoodSputumSubset.ComparedTo.LancepH7_Rv`,
                 `LancepH7_Rv.ComparedTo.Indigo_Rv`)

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
              "GoodSputumSubset_L4.ComparedTo.L2",
              
              "Indigo_Rv.ComparedTo.Ra",
              "GoodSputumSubset.ComparedTo.Indigo_Rv",
              "LancepH7_Rv.ComparedTo.Ra",
              "GoodSputumSubset.ComparedTo.LancepH7_Rv",
              "LancepH7_Rv.ComparedTo.Indigo_Rv")

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
################# EXPORT THE LIST OF DFs ##################
# Just export the ones I use in the table

wb <- createWorkbook()

addWorksheet(wb, "Sputum.vs.Broth")
writeData(wb, "Sputum.vs.Broth", list_dfs_f$GoodSputumSubset.ComparedTo.Broth %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))
addWorksheet(wb, "CaseumMimic.vs.Broth")
writeData(wb, "CaseumMimic.vs.Broth", list_dfs_f$CaseumMimic.ComparedTo.Broth %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))
addWorksheet(wb, "Rabbit.vs.Broth")
writeData(wb, "Rabbit.vs.Broth", list_dfs_f$Rabbit.ComparedTo.Broth %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))
addWorksheet(wb, "Marmoset.vs.Broth")
writeData(wb, "Marmoset.vs.Broth", list_dfs_f$Marmoset.ComparedTo.Broth %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))
addWorksheet(wb, "Sputum.vs.CaseumMimic")
writeData(wb, "Sputum.vs.CaseumMimic", list_dfs_f$GoodSputumSubset.ComparedTo.CaseumMimic %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))
addWorksheet(wb, "Sputum.vs.Rabbit")
writeData(wb, "Sputum.vs.Rabbit", list_dfs_f$GoodSputumSubset.ComparedTo.Rabbit %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))
addWorksheet(wb, "Sputum.vs.Marmoset")
writeData(wb, "Sputum.vs.Marmoset", list_dfs_f$GoodSputumSubset.ComparedTo.Marmoset %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))
addWorksheet(wb, "CaseumMimic.vs.Rabbit")
writeData(wb, "CaseumMimic.vs.Rabbit", list_dfs_f$CaseumMimic.ComparedTo.Rabbit %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))
addWorksheet(wb, "CaseumMimic.vs.Marmoset")
writeData(wb, "CaseumMimic.vs.Marmoset", list_dfs_f$CaseumMimic.ComparedTo.Marmoset %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))
addWorksheet(wb, "Rabbit.vs.Marmoset")
writeData(wb, "Rabbit.vs.Marmoset", list_dfs_f$Rabbit.ComparedTo.Marmoset %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))
addWorksheet(wb, "L4.vs.L2")
writeData(wb, "L4.vs.L2", list_dfs_f$GoodSputumSubset_L4.ComparedTo.L2 %>% select(GENE_NAME, GENE_ID, PRODUCT, LOG2FOLD, FDR_PVALUE, DE2_FDR, DE2_FDR_labels) %>% rename(DE = DE2_FDR, Labels = DE2_FDR_labels))

# saveWorkbook(wb, "DEG.xlsx", overwrite = TRUE)

###########################################################
############# IMPORT BOB's METAGENESETS DATA ##############


`MetaGeneSets_GoodSputumSubset.vs.Broth_UP` <- read.delim("Data/Differential_Expression/GoodSputumSubset_vs_Broth/W0.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_Marmoset.vs.Broth_UP` <- read.delim("Data/Differential_Expression/Marmoset_vs_Broth/Marmoset.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_Rabbit.vs.Broth_UP` <- read.delim("Data/Differential_Expression/Rabbit_vs_Broth/Rabbit.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_CaseumMimic.vs.Broth_UP` <- read.delim("Data/Differential_Expression/CaseumMimic_vs_Broth/CaseumMimic.MTb.MetaGeneSets.UP.txt")

`MetaGeneSets_GoodSputumSubset_L2.vs.Broth_UP` <- read.delim("Data/Differential_Expression/GoodSputumSubset.L2_vs_Broth/L2.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_GoodSputumSubset_L4.vs.Broth_UP` <- read.delim("Data/Differential_Expression/GoodSputumSubset.L4_vs_Broth/L4.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_GoodSputumSubset_L4.vs.L2_UP` <- read.delim("Data/Differential_Expression/GoodSputumSubset.L2_vs_L4/L4.MTb.MetaGeneSets.UP.txt")

`MetaGeneSets_Indigo_Rv.vs.Ra_UP` <- read.delim("Data/Differential_Expression/Ra_vs_Rv.Indigo/Rv.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_GoodSputumSubset.vs.Indigo_Rv_UP` <- read.delim("Data/Differential_Expression/SputumSubset_vs_Rv.Indigo/W0.MTb.MetaGeneSets.UP.txt")
`MetaGeneSetsLancepH7_Rv.vs.Ra_UP` <- read.delim("Data/Differential_Expression/Ra_vs_Rv.LancepH7/Rv_pH7.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_GoodSputumSubset.vs.LancepH7_Rv_UP` <- read.delim("Data/Differential_Expression/SputumSubset_vs_Rv.LancepH7/W0.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_LancepH7_Rv.vs.Indigo_Rv_UP` <- read.delim("Data/Differential_Expression/LancepH7.Rv_vs_Indigo_Rv/Rv_pH7.MTb.MetaGeneSets.UP.txt")



###########################################################
#################### SAVING FOR PAPER #####################

SputumVsBroth_iModulons <- MetaGeneSets_GoodSputumSubset.vs.Broth_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% # Because Bob's pipeline doesn't actually do this...
  # rename("Sputum_LOG2FOLD" = "LOG2FOLD", "Sputum_AVG.PVALUE" = "AVG_PVALUE", "Sputum_AVG.RANK" = "AVG_RANK") %>%
  # mutate()
  select(PathName, N_Genes, LOG2FOLD, FDR.pvalue)


# Save the one being used in the paper
addWorksheet(wb, "Sputum.vs.Broth_iModulons")
writeData(wb, "Sputum.vs.Broth_iModulons", SputumVsBroth_iModulons)
saveWorkbook(wb, "DEG.xlsx", overwrite = TRUE)
