# Calculate TPM from Raw Reads
# E. Lamont
# 8/1/25

# Want my own way to convert raw reads to TPM in case I want to do any filtering or batch correction on the raw reads

# https://support.bioconductor.org/p/91218/#91256
# https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

# source("Import_data.R") # To get All_RawReads
# 
# # Get all the gene lengths
# load("Data/MTb.MapSet.rda") # This is named mapSet
# # This has all the lengths:
# mapSet[["geneMap"]]
# my_geneLengths <- mapSet[["geneMap"]] %>% select(GENE_ID, NAME, N_EXON_BASES)
# 
# All_RawReads_2 <- All_RawReads %>% column_to_rownames("X")
# 
# # Order the geneLengths so they are in the same order as the RawReads
# my_geneLengths_ordered <- my_geneLengths[match(All_RawReads$X, my_geneLengths$GENE_ID), ]
# my_geneLengths_ordered <- my_geneLengths_ordered %>% mutate(Kilobases = N_EXON_BASES/1000)
# 
# # Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
# All_RPK <- All_RawReads_2 / my_geneLengths_ordered$Kilobases
# 
# # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
# ScalingFactor <- colSums(All_RPK) / 1e6
# 
# # Divide the RPK values by the “per million” scaling factor. This gives you TPM.
# All_tpm_byHand <- sweep(All_RPK, 2, ScalingFactor, FUN = "/")
# # This isn't exactly identical to the All_tpm from Bob but it is really close
# 
# # This below get the same thing but is a little more confusing
# # tpm.mat <- t( t(All_RPK) * 1e6 / colSums(All_RPK) )

###########################################################
#################### MAKE TPM FUNCTION ####################

CalculateTPM <- function(Raw_reads) {
  ## This takes the raw reads from Bob's pipeline and turns them into TPM
  ## Requires the MTb.MapSet.rda to be in the R environment
  
  # 1. Extract and organize the gene lengths 
  load("Data/MTb.MapSet.rda")
  my_geneLengths <- mapSet[["geneMap"]] %>% select(GENE_ID, NAME, N_EXON_BASES)
  my_geneLengths_ordered <- my_geneLengths[match(All_RawReads$X, my_geneLengths$GENE_ID), ]
  my_geneLengths_ordered <- my_geneLengths_ordered %>% mutate(Kilobases = N_EXON_BASES/1000)
  
  # 2. Convert column to rowname in the raw data
  Raw_reads_2 <- Raw_reads %>% column_to_rownames("X")
  
  # 3. Divide the raw reads by the gene lengths in KB
  All_RPK <- All_RawReads_2 / my_geneLengths_ordered$Kilobases
  
  # 4. Generate the "per million" scaling factor (sum of RPK values / 1e6 for each column)
  ScalingFactor <- colSums(All_RPK) / 1e6
  
  # 5. Divide the RPK by the Scaling Factor
  All_tpm <- sweep(All_RPK, 2, ScalingFactor, FUN = "/")
  
  return(All_tpm)
  
}

# testing <- CalculateTPM(All_RawReads)



