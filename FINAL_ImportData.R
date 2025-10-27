# Import Data

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(ggpmisc)
library(tidyverse)
library(ggpubr) 
library(ggrepel) 
library(ggcorrplot)
library(openxlsx)
library(readxl)
library(cowplot)

###########################################################
##################### IMPORT METADATA #####################

excel_sheets("metadata.xlsx")

LimitofDetect_pipeSummary <- read_excel("metadata.xlsx", sheet = "LimitofDetection")
CapturedVsNot_pipeSummary <- read_excel("metadata.xlsx", sheet = "CapturedVsNot")
BiolSamples_pipeSummary <- read_excel("metadata.xlsx", sheet = "BiolSamples")
ProbeTests_pipeSummary <- read_excel("metadata.xlsx", sheet = "ProbeTests")

# Rearrange
ordered_Probe <- c("Uncaptured", "Captured")
CapturedVsNot_pipeSummary$Probe <- factor(CapturedVsNot_pipeSummary$Probe, levels = ordered_Probe)
ordered_BiolSample <- c("Caseum mimic", "Rabbit", "Marmoset", "Sputum")
BiolSamples_pipeSummary$Type <- factor(BiolSamples_pipeSummary$Type, levels = ordered_BiolSample)

###########################################################
##################### IMPORT RAW READS ####################

excel_sheets("rawreads.xlsx")

All_RawReads <- read_excel("rawreads.xlsx", sheet = "All_RawReads")

###########################################################
############ SUM TOTAL READS ALIGNING TO Rv  ##############

NumRawReads <- round(colSums(All_RawReads %>% column_to_rownames("X")))

# Add N_Genomic
LimitofDetect_pipeSummary$N_Genomic_Rv <- NumRawReads[LimitofDetect_pipeSummary$SampleID]
CapturedVsNot_pipeSummary$N_Genomic_Rv <- NumRawReads[CapturedVsNot_pipeSummary$SampleID]
BiolSamples_pipeSummary$N_Genomic_Rv <- NumRawReads[BiolSamples_pipeSummary$SampleID]
ProbeTests_pipeSummary$N_Genomic_Rv <- NumRawReads[ProbeTests_pipeSummary$SampleID]

# Add P_Genomic
LimitofDetect_pipeSummary <- LimitofDetect_pipeSummary %>% mutate(P_Genomic_Rv = round((N_Genomic_Rv/RawReads)*100, 2))
CapturedVsNot_pipeSummary <- CapturedVsNot_pipeSummary %>% mutate(P_Genomic_Rv = round((N_Genomic_Rv/RawReads)*100, 2))
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>% mutate(P_Genomic_Rv = round((N_Genomic_Rv/RawReads)*100, 2))
ProbeTests_pipeSummary <- ProbeTests_pipeSummary %>% mutate(P_Genomic_Rv = round((N_Genomic_Rv/RawReads)*100, 2))

# N_NonCodingMtbRNA
LimitofDetect_pipeSummary <- LimitofDetect_pipeSummary %>% mutate(N_NonCodingMtbRNA = RawReads - N_Genomic_Rv - N_RiboClear - N_NoHit)
CapturedVsNot_pipeSummary <- CapturedVsNot_pipeSummary %>% mutate(N_NonCodingMtbRNA = RawReads - N_Genomic_Rv - N_RiboClear - N_NoHit)
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>% mutate(N_NonCodingMtbRNA = RawReads - N_Genomic_Rv - N_RiboClear - N_NoHit)
ProbeTests_pipeSummary <- ProbeTests_pipeSummary %>% mutate(N_NonCodingMtbRNA = RawReads - N_Genomic_Rv - N_RiboClear - N_NoHit)

LimitofDetect_pipeSummary <- LimitofDetect_pipeSummary %>% mutate(P_NonCodingMtbRNA = round((N_NonCodingMtbRNA/RawReads)*100, 2))
CapturedVsNot_pipeSummary <- CapturedVsNot_pipeSummary %>% mutate(P_NonCodingMtbRNA = round((N_NonCodingMtbRNA/RawReads)*100, 2))
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>% mutate(P_NonCodingMtbRNA = round((N_NonCodingMtbRNA/RawReads)*100, 2))
ProbeTests_pipeSummary <- ProbeTests_pipeSummary %>% mutate(P_NonCodingMtbRNA = round((N_NonCodingMtbRNA/RawReads)*100, 2))


###########################################################
################ CALCULATE TXN COVERAGE  ##################

# Count, for each column (sample), how many genes have >= 10 reads
NumGoodReads <- colSums(All_RawReads >= 10)

# Add as new column in pipeSummary, matching by SampleID
LimitofDetect_pipeSummary$AtLeast.10.Reads <- NumGoodReads[LimitofDetect_pipeSummary$SampleID]
CapturedVsNot_pipeSummary$AtLeast.10.Reads <- NumGoodReads[CapturedVsNot_pipeSummary$SampleID]
BiolSamples_pipeSummary$AtLeast.10.Reads <- NumGoodReads[BiolSamples_pipeSummary$SampleID]
ProbeTests_pipeSummary$AtLeast.10.Reads <- NumGoodReads[ProbeTests_pipeSummary$SampleID]

# Add transcriptional coverage
LimitofDetect_pipeSummary <- LimitofDetect_pipeSummary %>% mutate(Txn_Coverage = round(AtLeast.10.Reads/4030*100, 2))
CapturedVsNot_pipeSummary <- CapturedVsNot_pipeSummary %>% mutate(Txn_Coverage = round(AtLeast.10.Reads/4030*100, 2))
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>% mutate(Txn_Coverage = round(AtLeast.10.Reads/4030*100, 2))
ProbeTests_pipeSummary <- ProbeTests_pipeSummary %>% mutate(Txn_Coverage = round(AtLeast.10.Reads/4030*100, 2))

###########################################################
############ EXPORT COMPLETE PIPE SUMMARIES  ##############

pipeSummary_wb <- createWorkbook()
addWorksheet(pipeSummary_wb, "CapturedVsNot")
writeData(pipeSummary_wb, "CapturedVsNot", CapturedVsNot_pipeSummary)
addWorksheet(pipeSummary_wb, "LimitofDetection")
writeData(pipeSummary_wb, "LimitofDetection", LimitofDetect_pipeSummary)
addWorksheet(pipeSummary_wb, "BiolSamples")
writeData(pipeSummary_wb, "BiolSamples", BiolSamples_pipeSummary)
addWorksheet(pipeSummary_wb, "ProbeTests")
writeData(pipeSummary_wb, "ProbeTests", ProbeTests_pipeSummary)

# Save the workbook
saveWorkbook(pipeSummary_wb, "Metadata_2025.10.27_v2.xlsx", overwrite = FALSE)






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


###########################################################
####################### IMPORT DEG ########################

file <- "Data.xlsx"

# Get all sheet names
sheet_names <- getSheetNames(file)

# Read each sheet into a list
DEG_dfs <- lapply(sheet_names, function(s) read.xlsx(file, sheet = s))

# Name the list elements
names(DEG_dfs) <- sheet_names



###########################################################
############ iMODULONS: MAKE LISTS OF GROUPS ##############

CentralCarbon_iModulons <- c("Peptidoglycan Biosynthesis", "Central Carbon Metabolism", "Fumarate Reductase", "PrpR", "BkaR", "Nicotinate Metabolism")
CentralCarbon_iModulons_pattern <- str_c(CentralCarbon_iModulons, collapse = "|") #
AminoAcid_iModulons <- c("GroEL-GroES Complex", "Leucine Related", "LysG", "ArgR") 
AminoAcid_iModulons_pattern <- str_c(AminoAcid_iModulons, collapse = "|") 
NucleicAcid_iModulons <- c("PyrR", "Rv0135\\+Rv1019", "Nucleic Acid Hydrolysis") 
NucleicAcid_iModulons_pattern <- str_c(NucleicAcid_iModulons, collapse = "|") 
FattyAcid.Cholesterol_iModulons <- c("Fatty Acid Biosynthesis", "KstR2", "Mycofactocin Synthesis Pathway", "FasR", "Polyketide Synthase Complex", "Rv0681") 
FattyAcid.Cholesterol_iModulons_pattern <- str_c(FattyAcid.Cholesterol_iModulons, collapse = "|") 
Metal_iModulons <- c("RicR", "IdeR", "M-box", "Zur", "Hpt-2b Induced") 
Metal_iModulons_pattern <- str_c(Metal_iModulons, collapse = "|") 
SulfurMetabolism_iModulons <- c("Sulfur Metabolism") 
SulfurMetabolism_iModulons_pattern <- str_c(SulfurMetabolism_iModulons, collapse = "|") 
Growth_iModulons <- c("Positive Regulation of Growth") 
Growth_iModulons_pattern <- str_c(Growth_iModulons, collapse = "|") 
Redox_iModulons <- c("DevR-1", "WhiB4", "DevR-2", "WhiB1", "WhiB4/IdeR", "Rv1828/SigH", "Rv1776c\\+WhiB4", "VirS", "WhiB6") 
Redox_iModulons_pattern <- str_c(Redox_iModulons, collapse = "|") 
AcidStress_iModulons <- c("MarR") 
AcidStress_iModulons_pattern <- str_c(AcidStress_iModulons, collapse = "|") 
Antibiotic_iModulons <- c("Lsr2", "Blal", "Rv0078\\+Rv2034", "WhiB7", "IniR") 
Antibiotic_iModulons_pattern <- str_c(Antibiotic_iModulons, collapse = "|") 
Virulence.Persistence_iModulons <- c("Rv0576", "Mce1R", "SigH", "PhoP", "Mce3R", "MprA", "PDIM\\;PGL Synthesis", "Rv2488c", "SigC", "SigD", "MbcA\\+Rv3249c\\+Rv3066", "SigK") 
Virulence.Persistence_iModulons_pattern <- str_c(Virulence.Persistence_iModulons, collapse = "|") 






