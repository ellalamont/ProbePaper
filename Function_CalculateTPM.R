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