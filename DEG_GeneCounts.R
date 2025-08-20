# Count the number of up and downregulated genes between each condition
# E. Lamont
# 8/5/25

source("Import_DEG_sets.R")

###########################################################
################## LOG2FOLD >=2 TABLE #####################
# Make an empty list to hold the count summaries
DEG2count_list <- list()

# Loop through list of dataframes
for (list_name in names(list_dfs_2)) {
  
  # Get the current dataframe
  df <- list_dfs_2[[list_name]]
  
  # Count the DE_2 values
  GeneCounts <- df %>%
    count(DE_2, name = "Count") %>% # Does the actual counting 
    pivot_wider(names_from = DE_2, values_from = Count, values_fill = 0)
  
  # Add back the name of the current list
  GeneCounts <- GeneCounts %>%
    mutate(DEG_Group = list_name) %>%
    select(DEG_Group, everything()) # Puts list name first
  
  # Add the results to my list
  DEG2count_list[[list_name]] <- GeneCounts
    
}

# Convert the list of dataframes into one dataframe
DEG2count_df <- bind_rows(DEG2count_list)

DEG2count_df <- DEG2count_df %>% mutate(All_Signifiant_DEG2 = `significant down` + `significant up`)


### Move to excel and make the rest by hand!!

###########################################################
############# LOG2FOLD >=2 AND FDR PVALUE TABLE ###########
# Make an empty list to hold the count summaries
DEG2_FDR_count_list <- list()

# Loop through list of dataframes
for (list_name in names(list_dfs_2)) {
  
  # Get the current dataframe
  df <- list_dfs_2[[list_name]]
  
  # Count the DE2_FDR values
  GeneCounts <- df %>%
    count(DE2_FDR, name = "Count") %>% # Does the actual counting 
    pivot_wider(names_from = DE2_FDR, values_from = Count, values_fill = 0)
  
  # Add back the name of the current list
  GeneCounts <- GeneCounts %>%
    mutate(DEG_Group = list_name) %>%
    select(DEG_Group, everything()) # Puts list name first
  
  # Add the results to my list
  DEG2_FDR_count_list[[list_name]] <- GeneCounts
  
}

# Convert the list of dataframes into one dataframe
DEG2_FDR_count_df <- bind_rows(DEG2_FDR_count_list)

DEG2_FDR_count_df <- DEG2_FDR_count_df %>% mutate(All_Signifiant_DEG2_FDR = `significant down` + `significant up`)



###########################################################
####### FILTERED LOG2FOLD >=2 AND FDR PVALUE TABLE ########
# Make an empty list to hold the count summaries
DEG2_FDR_count_list_f <- list()

# Loop through list of dataframes
for (list_name in names(list_dfs_f)) {
  
  # Get the current dataframe
  df <- list_dfs_f[[list_name]]
  
  # Count the DE2_FDR values
  GeneCounts <- df %>%
    count(DE2_FDR, name = "Count") %>% # Does the actual counting 
    pivot_wider(names_from = DE2_FDR, values_from = Count, values_fill = 0)
  
  # Add back the name of the current list
  GeneCounts <- GeneCounts %>%
    mutate(DEG_Group = list_name) %>%
    select(DEG_Group, everything()) # Puts list name first
  
  # Add the results to my list
  DEG2_FDR_count_list_f[[list_name]] <- GeneCounts
  
}

# Convert the list of dataframes into one dataframe
DEG2_FDR_count_df_f <- bind_rows(DEG2_FDR_count_list_f)

DEG2_FDR_count_df_f <- DEG2_FDR_count_df_f %>% mutate(DEG2_FDR_count_df_f = `significant down` + `significant up`)


### Move to excel and make the rest by hand!!



