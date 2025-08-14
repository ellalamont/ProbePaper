# Import data for the probe MS
# 7/31/25

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)
# library(ggprism) # for add_pvalue()
library(rstatix) # for adjust_pvalue
library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
library(pheatmap)
# library(dendextend) # May need this for looking at pheatmap clustering
library(ggplotify) # To convert pheatmaps to ggplots
library(corrplot)
library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly
library(edgeR) # for cpm
library(sva) # For ComBat_seq batch correction
# devtools::install_github("NightingaleHealth/ggforestplot")
library(ggforestplot)
library(ggforce)

# DuffyTools
library(devtools)
# install_github("robertdouglasmorrison/DuffyTools")
library(DuffyTools)
# install_github("robertdouglasmorrison/DuffyNGS")
# BiocManager::install("robertdouglasmorrison/DuffyTools")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biobase") 


###########################################################
########################## COLORS #########################

# Labelled Colors
my_fav_colors <- c(`Sputum` = "#0072B2", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00")

# Labelled Shapes
my_fav_shapes <- c(`Sputum` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette4 <- c("#56B4E9", "#009E73", "#F0E442","#CC79A7")
c25 <- c(
  "dodgerblue2", "#E31A1C", "green4",
  "#6A3D9A","#FF7F00","black", "gold1",
  "skyblue2", "#FB9A99","palegreen2","#CAB2D6",
  "#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown"
)
c12 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "palegreen2", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4")
c16 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black","gold1", "#FB9A99", "#CAB2D6", "palegreen2", "gray70", "maroon", "orchid1", "blue1", "darkturquoise", "darkorange4")


# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default


###########################################################
############### IMPORT PIPELINE SUMMARY DATA ##############

## PREDICTTB_RUN1
# This has been edited to include more metadata!
Run1_pipeSummary <- read.csv("Data/PredictTB_Run1/Pipeline.Summary.Details.csv") 
# Just get the samples I want 
Run1_pipeSummary <- Run1_pipeSummary %>% filter(Type != "TBAIT")
Run1_pipeSummary <- Run1_pipeSummary %>% filter(Type != 'NA')
# Add a column for the Run
Run1_pipeSummary <- Run1_pipeSummary %>% mutate(Run = "PredictTB_Run1")

# Randomize and subset to get 12 sputum samples
# First remove all the relapse 
sputumOnly_pipeSummary <- Run1_pipeSummary %>% filter(str_detect(Type, "sputum")) %>% filter(Outcome == "Cure" & Week == "Week 0")
# set.seed(2)
# set.seed(5)
set.seed(4)
# SputumSubset_pipeSummary <- slice_sample(sputumOnly_pipeSummary, n = 12)
# SputumSubset_list <- SputumSubset_pipeSummary %>% pull(SampleID)
# Save this list just in case it keeps changing
# SputumSubset_list <- c("W0_13004_S53", "W0_12034_S29", "W0_12009_S4", "W0_15081_S49", "W0_14113_S41", "W0_12012_S6", "W0_12038_S31", "W0_14044_S36", "W0_15035_S25", "W0_12029_S10", "W0_12028_S9", "W0_11011_S16") # set seed 2
# SputumSubset_list <- c("W0_12028_S9", "W0_11047_S22",  "W0_12010_S5", "W0_14048_S37", "W0_15045_S26", "W0_11012_S17", "W0_13027_S15", "W0_11011_S16", "W0_12083_S39", "W0_14005_S55", "W0_15065_S47", "W0_14113_S41") # set seed 4
# Actually now thinking I should only inlucde samples with TOTAL reads > 1 million (everything else is failed)
sputumOnly_pipeSummary2 <- sputumOnly_pipeSummary %>% filter(RawReads >= 1e6)
set.seed(42)
SputumSubset_pipeSummary <- slice_sample(sputumOnly_pipeSummary2, n = 12)
SputumSubset_list <- SputumSubset_pipeSummary %>% pull(SampleID)
SputumSubset_list <- c("W0_15081_S49", "W0_11011_S16", "W0_13045_S34", "W0_12024_S8", "W0_12043_S32", "W0_13027_S15", "W0_12010_S5", "W0_12007_S2", "W0_15089_S51", "W0_12008_S3", "W0_12032_S28", "W0_12083_S39")



### PREDICTTB RUN 1 MIMIC and RABBIT DATA ###
MimicRabbit_pipeSummary <- Run1_pipeSummary %>% filter(Type %in% c("Caseum mimic", "Rabbit")) %>% filter(RawReads >=1e6)

### MARMOSET DATA FROM PROBETEST 3 ###
# Pull out the marmoset information from ProbeTest3
ProbeTest3_pipeSummary <- read.csv("Data/ProbeTest3/ProbeTest3_Pipeline.Summary.Details.csv") 
Marmoset_pipeSummary <- ProbeTest3_pipeSummary %>% filter(SampleID %in% c("BQ12_10_Probe_3A_S29", "BQ12_3_Probe_4A_50_S27", "BQ12_8_Probe_4A_50_S28")) %>% filter(RawReads >= 1e6)
Marmoset_pipeSummary <- Marmoset_pipeSummary %>% mutate(Run = "ProbeTest3")

### BROTH DATA FROM PROBETEST 5 ###
ProbeTest5_pipeSummary <- read.csv("Data/ProbeTest5/ProbeTest5_Pipeline.Summary.Details_moreTrim.csv") 
Broth_pipeSummary <- ProbeTest5_pipeSummary %>% filter(SampleID %in% c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9"))

### CAPTURED VS NOT DATA FROM PROBETEST 5 ###
CapturedVsNot_pipeSummary <- ProbeTest5_pipeSummary %>% 
  filter(SampleID %in% c("THP1_1e6_1a_S28", "THP1_1e6_1b_S29", "THP1_1e6_2a_S30", "THP1_1e6_2b_S31", "THP1_1e6_3a_S32", "THP1_1e6_3b_S33")) %>%
  mutate(Replicates = c("R1","R1","R2","R2","R3","R3"))
ordered_Probe <- c("None", "JA2")
CapturedVsNot_pipeSummary$Probe <- factor(CapturedVsNot_pipeSummary$Probe, levels = ordered_Probe)

### LIMIT OF DETECTION DATA FROM PROBETEST 5 ###
LimitofDetect_pipeSummary <- ProbeTest5_pipeSummary %>% 
  filter(Sample_Type == "THP1") %>% 
  filter(Ra_cells != "none") %>% 
  filter(Probe != "None") %>%
  mutate(Ra_cells2 = case_when(
    Ra_cells == "one_e_2" ~ "1e2",
    Ra_cells == "one_e_3" ~ "1e3",
    Ra_cells == "one_e_4" ~ "1e4",
    Ra_cells == "one_e_5" ~ "1e5",
    Ra_cells == "one_e_6" ~ "1e6",
    Ra_cells == "one_e_8" ~ "1e8",
    TRUE ~ "default_value"  # this is optional for values that don't meet any condition
  ))
# Should maybe go back and remove all the 1e6 cells that are not specifically for the limit of detection, but not doing that right now.
# Or what might be easier is just include the samples from ProbeTest5, where I have at least 3 replicates each! 


### MERGE the pipeSummaries for Biological Samples ###
BiolSamples_pipeSummary <- merge(SputumSubset_pipeSummary, Marmoset_pipeSummary, all = T)
BiolSamples_pipeSummary <- merge(BiolSamples_pipeSummary, MimicRabbit_pipeSummary, all = T)
BiolSamples_pipeSummary <- merge(BiolSamples_pipeSummary, Broth_pipeSummary, all = T)

# Merge two columns
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>% mutate(Type = coalesce(Type, Sample_Type)) %>%
  mutate(Type2 = coalesce(Type2, Sample_Type)) %>% 
  select(-Sample_Type)
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>% mutate(CFU_per_g.or.mL = coalesce(CFU_per_g, CFU_per_mL))
BiolSamples_pipeSummary$X <- NULL
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>%
  mutate(Type = case_when(Type == "Week 0 sputum" ~ "Sputum",
                          TRUE ~ Type))

# Reorder things
ordered_Type <- c("Caseum mimic", "Rabbit", "Marmoset", "Sputum", "Broth")
BiolSamples_pipeSummary$Type <- factor(BiolSamples_pipeSummary$Type, levels = ordered_Type)



# Already did the below in excel
# All_pipeSummary$SampleID <- gsub(x = All_pipeSummary$SampleID, pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# All_pipeSummary <- All_pipeSummary %>% mutate(Sputum_Number = str_extract(SampleID, "S_[0-9]+"))

###########################################################
########### LIST OF SAMPLES PASSING INSPECTION ############

GoodSampleList <- BiolSamples_pipeSummary %>%
  filter(N_Genomic >= 1000000 & AtLeast.10.Reads >= (4499*0.8)) %>% 
  pull(SampleID)

Sputum_GoodSampleList <- GoodSampleList[grep("^W", GoodSampleList)]

Marm_GoodSampleList <- GoodSampleList[grep("^BQ", GoodSampleList)]

Mimic_GoodSampleList <- GoodSampleList[grep("^HN878", GoodSampleList)]

Rabbit_GoodSampleList <- GoodSampleList[grep("Cav", GoodSampleList)]

###########################################################
############### IMPORT AND PROCESS TPM VALUES #############

Run1_tpm <- read.csv("Data/PredictTB_Run1/Mtb.Expression.Gene.Data.TPM.csv")
Run1_tpm <- Run1_tpm %>% select(-contains("TBAIT"))

Run1_SputumSubset_tpm <- Run1_tpm %>% select(X, all_of(SputumSubset_list))

# Just pull the tpm of the THP1 spiked from another run: THP1 1e6_1 (Predict rack 2 box 1 I04)
# Need THP1 1e6_1a from the Januaray run. Also need 
ProbeTest5_tpm <- read.csv("Data/ProbeTest5/ProbeTest5_Mtb.Expression.Gene.Data.TPM_moreTrim.csv") 
ProbeTest5_tpm_subset <- ProbeTest5_tpm %>% select(X, THP1_1e6_1a_S28)
ProbeTest5_tpm_Broth <- ProbeTest5_tpm %>% select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9)
ProbeTest5_tpm_CapturedVsNot_wBroth <- ProbeTest5_tpm %>% select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9, THP1_1e6_1a_S28, THP1_1e6_1b_S29, THP1_1e6_2a_S30, THP1_1e6_2b_S31, THP1_1e6_3a_S32, THP1_1e6_3b_S33)

# Get the Marmoset TPM
ProbeTest3_tpm <- read.csv("Data/ProbeTest3/ProbeTest3_Mtb.Expression.Gene.Data.TPM.csv")
ProbeTest3_tpm_marm <- ProbeTest3_tpm %>% select(X, BQ12_10_Probe_3A_S29, BQ12_3_Probe_4A_50_S27, BQ12_8_Probe_4A_50_S28)

# PROBETEST 4
ProbeTest4_tpm <- read.csv("Data/ProbeTest4/Mtb.Expression.Gene.Data.TPM.csv")



# Adjust the names so they are slightly shorter
# names(All_tpm) <- gsub(x = names(All_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# rownames(All_tpm) <- All_tpm[,1] # add the rownames

# Need to make sure there is a Gene column (gets lots)
# All_tpm <- All_tpm %>% 
#   rename(Gene = X) 


###########################################################
###### MAKE TPM WITH ALL CLINICAL AND ANIMAL MODELS #######

# Merge the tpms I collected above
All_tpm <- merge(Run1_tpm, ProbeTest3_tpm_marm, all = T)
All_tpm <- merge(All_tpm, ProbeTest5_tpm_Broth)

# Just keep the samples passing filter
GoodBiolSamples_tpm <- All_tpm %>% select("X", all_of(GoodSampleList))



###########################################################
########### FILTER TPM TO REMOVE NON CODING RNA ###########

# Not sure if this is right but will do it for now
# Removing the Rvncs (non coding RNAs) because I don't think the current probes contain them! (Previous probes have)
## I think this because they are missing in PredictTB Run 1 THP1 spiked sample but present in that same sample from ProbeTest5
# The Pathcap people also had issues with ncRNAs: https://www.nature.com/articles/s41598-019-55633-6#Sec8

GoodBiolSamples_tpmF <- GoodBiolSamples_tpm %>% filter(!str_detect(X, regex("Rvnc", ignore_case = T)))


###########################################################
############### IMPORT AND PROCESS RAW READS ##############

Run1_RawReads <- read.csv("Data/PredictTB_Run1/Mtb.Expression.Gene.Data.readsM.csv")
Run1_RawReads <- Run1_RawReads %>% select(-contains("TBAIT"))

ProbeTest5_RawReads <- read.csv("Data/ProbeTest5/ProbeTest5_Mtb.Expression.Gene.Data.readsM_moreTrim.csv") 
ProbeTest5_RawReads_Broth <- ProbeTest5_RawReads %>% select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9)

ProbeTest3_RawReads <- read.csv("Data/ProbeTest3/Mtb.Expression.Gene.Data.readsM_old.csv")
ProbeTest3_RawReads_marm <- ProbeTest3_RawReads %>% select(X, BQ12_10_Probe_3A_S29, BQ12_3_Probe_4A_50_S27, BQ12_8_Probe_4A_50_S28)

# Merge the RawReads I collected above
All_RawReads <- merge(Run1_RawReads, ProbeTest3_RawReads_marm, all = T)
All_RawReads <- merge(All_RawReads, ProbeTest5_RawReads_Broth)

# Just keep the samples passing filter
GoodBiolSamples_RawReads <- All_RawReads %>% select("X", all_of(GoodSampleList), "H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9")

# Just keep the sputum samples
GoodSputumSubset_RawReads <- GoodBiolSamples_RawReads %>% select("X", any_of(SputumSubset_list))
