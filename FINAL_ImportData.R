# Import Data

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr) 
library(ggrepel) 
library(ggcorrplot)
library(openxlsx)
library(readxl)


###########################################################
##################### IMPORT METADATA #####################

excel_sheets("metadata.xlsx")

LimitofDetect_pipeSummary <- read_excel("metadata.xlsx", sheet = "LimitofDetection")
CapturedVsNot_pipeSummary <- read_excel("metadata.xlsx", sheet = "CapturedVsNot")
BiolSamples_pipeSummary <- read_excel("metadata.xlsx", sheet = "BiolSamples")

# Rearrange
ordered_Probe <- c("Uncaptured", "Captured")
CapturedVsNot_pipeSummary$Probe <- factor(CapturedVsNot_pipeSummary$Probe, levels = ordered_Probe)


###########################################################
##################### IMPORT RAW READS ####################

excel_sheets("rawreads.xlsx")

All_RawReads <- read_excel("rawreads.xlsx", sheet = "All_RawReads")


###########################################################
################ CALCULATE TXN COVERAGE  ##################

# Count, for each column (sample), how many genes have >= 10 reads
NumGoodReads <- colSums(All_RawReads >= 10)

# Add as new column in pipeSummary, matching by SampleID
LimitofDetect_pipeSummary$AtLeast.10.Reads <- NumGoodReads[LimitofDetect_pipeSummary$SampleID]
CapturedVsNot_pipeSummary$AtLeast.10.Reads <- NumGoodReads[CapturedVsNot_pipeSummary$SampleID]
BiolSamples_pipeSummary$AtLeast.10.Reads <- NumGoodReads[BiolSamples_pipeSummary$SampleID]

# Add transcriptional coverage
LimitofDetect_pipeSummary <- LimitofDetect_pipeSummary %>% mutate(Txn_Coverage = round(AtLeast.10.Reads/4030*100, 2))
CapturedVsNot_pipeSummary <- CapturedVsNot_pipeSummary %>% mutate(Txn_Coverage = round(AtLeast.10.Reads/4030*100, 2))
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>% mutate(Txn_Coverage = round(AtLeast.10.Reads/4030*100, 2))

###########################################################
##################### CONVERT TO TPM ######################

CalculateTPM_RvOnly <- function(Raw_reads) {
  ## This takes the raw reads from Bob's pipeline and turns them into TPM
  ## Requires the MTb.MapSet.rda to be in the R environment
  ## The genes are filtered here to only include protein coding Rv genes
  ## Gene names in column X
  
  # 1. Extract and organize the gene lengths 
  load("Data/MTb.MapSet.rda")
  my_geneLengths <- mapSet[["geneMap"]] %>% select(GENE_ID, NAME, N_EXON_BASES)
  
  my_geneLengths_f <- my_geneLengths %>%
    filter(grepl("^Rv[0-9]+[A-Za-z]?$", GENE_ID))
  
  my_geneLengths_ordered <- my_geneLengths_f[match(Raw_reads$X, my_geneLengths_f$GENE_ID), ]
  my_geneLengths_ordered <- my_geneLengths_ordered %>% mutate(Kilobases = N_EXON_BASES/1000)
  
  if (any(is.na(my_geneLengths_ordered$Kilobases))) {
    warning("Some genes in Raw_reads did not match gene lengths")
  }
  
  # 2. Convert column to rowname in the raw data
  Raw_reads_2 <- Raw_reads %>% column_to_rownames("X")
  
  # 3. Divide the raw reads by the gene lengths in KB
  All_RPK <- Raw_reads_2 / my_geneLengths_ordered$Kilobases
  
  # 4. Generate the "per million" scaling factor (sum of RPK values / 1e6 for each column)
  ScalingFactor <- colSums(All_RPK) / 1e6
  
  # 5. Divide the RPK by the Scaling Factor
  All_tpm <- sweep(All_RPK, 2, ScalingFactor, FUN = "/")
  
  return(All_tpm)
  
}

All_tpm <- CalculateTPM_RvOnly(All_RawReads)


