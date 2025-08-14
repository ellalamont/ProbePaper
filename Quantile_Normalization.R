# Quantile Normalization on my W0 TPM values
# 4/16/25

# Following Coppola et al. 2021
# https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.763364/full#h13

source("Import_data.R")

# GoodSputumSubset_RawReads
# Won't work if starting with averages! Needs multiple columns

# Convert the X to rownames
GoodSputumSubset_RawReads_2 <- GoodSputumSubset_RawReads %>% column_to_rownames(("X"))

################################################
###### FUNCTION FROM COPPOLA ET AL. 2021 #######

quantile_normalisation <- function(df)
  
{df_rank <- apply (df,2, rank, ties.method="min")

df_sorted <- data.frame (apply (df, 2, sort))

df_mean <- apply (df_sorted, 1, mean)

index_to_mean <- function (my_index, my_mean) {return (my_mean [my_index])}

df_final <- apply (df_rank, 2, index_to_mean, my_mean = df_mean)

rownames (df_final) <- rownames(df)

return(df_final)}

test <- as.data.frame(quantile_normalisation(GoodSputumSubset_RawReads_2))
class(test)


#########################################################
###### DO THIS ALL WITH THE GOOD SPUTUM RAW READS #######
# GoodSputumSubset_RawReads

GoodSputumSubset_RawReads_QuantileNormalization <- as.data.frame(quantile_normalisation(GoodSputumSubset_RawReads_2))

# ChatGPT made this to change the Quantile Normalization to a 0-100 scale
rank_cols_to_percentile <- function(df) {
  apply(df, 2, function(col) {
    ranks <- rank(col, ties.method = "min")
    percentiles <- (ranks - 1) / (length(ranks) - 1) * 100
    return(percentiles)
  }) |> as.data.frame()
}

GoodSputumSubset_RawReads_QuantileNormalization_0.100 <- rank_cols_to_percentile(GoodSputumSubset_RawReads_QuantileNormalization)

# Get the average for my 9 sputum samples
RANK_Average_RawReads_GoodSputumSubset <- GoodSputumSubset_RawReads_QuantileNormalization_0.100 %>% 
  mutate(RANK_Average = rowMeans(across(where(is.numeric)))) %>%
  select(RANK_Average)

