# Look at gene length and if that is a cause of missing genes
# E. Lamont
 # 7/31/25

# This is after I saw that Rvnc0009 (mpr6) is not present in any of my sputum samples
# But I know it is being captured because it is showing up in other samples!


source("Import_data.R") # for BiolSamples_tpm

load("Data/MTb.MapSet.rda") # This is named mapSet


# I think this has all the lengths:
mapSet[["geneMap"]]

my_geneLengths <- mapSet[["geneMap"]] %>% select(GENE_ID, NAME, N_EXON_BASES)


# Merge the two dfs

my_tpm <- BiolSamples_tpm %>% rename(GENE_ID = X)
merge_lengths <- full_join(my_geneLengths, my_tpm, by = "GENE_ID")

# Lets just put in excel for filtering
write.csv(merge_lengths, "TPM_w_GeneLength.csv")
