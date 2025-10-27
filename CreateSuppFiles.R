# Combine all the files to make clean supplemental files for import



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

# THP1 CONTROL FROM PROBETEST 5
THP1_Control_pipeSummary <- ProbeTest5_pipeSummary %>% 
  filter(str_detect(SampleID, "Control"))

# TWO PROBE PREPS 
ProbeA_pipeSummary <- ProbeTest3_pipeSummary %>% filter(SampleID %in% c("THP1_1e6_1_Probe_3D_100_S1", "THP1_1e6_2_Probe_3D_50_S2", "THP1_1e6_3_Probe_3D_25_S3", "THP1_1e6_4_Probe_3D_10_S4", "THP1_1e6_5_Probe_1_S5")) %>% 
  mutate(Run = "ProbeTest3")
ProbeTest4_pipeSummary <- read.csv("Data/ProbeTest4/Pipeline.Summary.Details.csv") 
ProbeB_pipeSummary <- ProbeTest4_pipeSummary %>% filter(SampleID %in% c("THP1_1e6_1_S41", "THP1_1e6_2_S42", "THP1_1e6_3_S43", "THP1_1e6_4_S44", "THP1_1e6_5_S45")) %>% 
  mutate(Run = "ProbeTest4")
ProbeTestAB_pipeSummary <- merge(ProbeA_pipeSummary, ProbeB_pipeSummary, all = T)


### MERGE the pipeSummaries for Biological Samples ###
BiolSamples_pipeSummary <- merge(SputumSubset_pipeSummary, Marmoset_pipeSummary, all = T)
BiolSamples_pipeSummary <- merge(BiolSamples_pipeSummary, MimicRabbit_pipeSummary, all = T)
BiolSamples_pipeSummary <- merge(BiolSamples_pipeSummary, Broth_pipeSummary, all = T)

# Merge two columns
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>% 
  mutate(Type = coalesce(Type, Sample_Type)) %>%
  mutate(Type2 = coalesce(Type2, Sample_Type)) %>% 
  mutate(Type3 = coalesce(Type3, Sample_Type)) %>%
  select(-Sample_Type)
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>% mutate(CFU_per_g.or.mL = coalesce(CFU_per_g, CFU_per_mL))
BiolSamples_pipeSummary$X <- NULL
BiolSamples_pipeSummary <- BiolSamples_pipeSummary %>%
  mutate(Type = case_when(Type == "Week 0 sputum" ~ "Sputum",
                          TRUE ~ Type))

# Reorder things
ordered_Type <- c("Caseum mimic", "Rabbit", "Marmoset", "Sputum", "Broth")
BiolSamples_pipeSummary$Type <- factor(BiolSamples_pipeSummary$Type, levels = ordered_Type)

# Export all the pipeSummary as excel in different tabs
install.packages("openxlsx")
library(openxlsx)

pipeSummary_wb <- createWorkbook()
addWorksheet(pipeSummary_wb, "CapturedVsNot")
writeData(pipeSummary_wb, "CapturedVsNot", CapturedVsNot_pipeSummary)
addWorksheet(pipeSummary_wb, "LimitofDetection")
writeData(pipeSummary_wb, "LimitofDetection", LimitofDetect_pipeSummary)
addWorksheet(pipeSummary_wb, "BiolSamples")
writeData(pipeSummary_wb, "BiolSamples", BiolSamples_pipeSummary)
addWorksheet(pipeSummary_wb, "ProbeTests")
writeData(pipeSummary_wb, "ProbeTests", ProbeTestAB_pipeSummary)

# Save the workbook
saveWorkbook(pipeSummary_wb, "Metadata_2025.10.27.xlsx", overwrite = FALSE)
# This has then been edited in excel


###########################################################
############### IMPORT AND PROCESS RAW READS ##############

Run1_RawReads <- read.csv("Data/PredictTB_Run1/Mtb.Expression.Gene.Data.readsM.csv")
Run1_RawReads <- Run1_RawReads %>% 
  select(X, HN878_mimic_D14_R1_S63, HN878_mimic_D14_R2_S64, HN878_mimic_D28_R1_S65, HN878_mimic_D28_R2_S66, HN878_mimic_D7_R2_S62, LLL_Cav_L2_18_S58, LLL_Cav_L2_21_S59, LU_Cav_L1_12_S57, RLL_Cav_Cluester_R4_23_S60, RUL_Cab_R1_26_S61,
         W0_15081_S49, W0_11011_S16, W0_13045_S34, W0_12024_S8, W0_12043_S32, W0_13027_S15,
         W0_12010_S5, W0_12007_S2, W0_15089_S51, W0_12008_S3, W0_12032_S28, W0_12083_S39)


ProbeTest5_RawReads <- read.csv("Data/ProbeTest5/ProbeTest5_Mtb.Expression.Gene.Data.readsM_moreTrim.csv") 
ProbeTest5_RawReads_CapturedVsNot_Broth_LimitofDetect <- ProbeTest5_RawReads %>% 
  select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9, 
         THP1_1e6_1a_S28, THP1_1e6_1b_S29, THP1_1e6_2a_S30, THP1_1e6_2b_S31, 
         THP1_1e6_3a_S32, THP1_1e6_3b_S33, THP1_1e6_4_S34,
         THP1_1e2_1_S14, THP1_1e2_2_S15, THP1_1e2_3_S16, THP1_1e2_4_S17, THP1_1e3_1_S18, THP1_1e3_2_S19, THP1_1e3_3_S20, THP1_1e3_4_S21, THP1_1e4_1_S22, THP1_1e4_2_S23, THP1_1e4_3_S24, THP1_1e5_1_S25, THP1_1e5_2_S26, THP1_1e5_3_S27)

ProbeTest3_RawReads <- read.csv("Data/ProbeTest3/Mtb.Expression.Gene.Data.readsM_old.csv")
ProbeTest3_RawReads <- ProbeTest3_RawReads %>% 
  select(X, BQ12_10_Probe_3A_S29, BQ12_3_Probe_4A_50_S27, BQ12_8_Probe_4A_50_S28, THP1_1e6_1_Probe_3D_100_S1, THP1_1e6_2_Probe_3D_50_S2, THP1_1e6_3_Probe_3D_25_S3, THP1_1e6_4_Probe_3D_10_S4, THP1_1e6_5_Probe_1_S5) %>%
  rename(ProbeA_THP1_1e6_1 = THP1_1e6_1_Probe_3D_100_S1,
         ProbeA_THP1_1e6_2 = THP1_1e6_2_Probe_3D_50_S2,
         ProbeA_THP1_1e6_3 = THP1_1e6_3_Probe_3D_25_S3,
         ProbeA_THP1_1e6_4 = THP1_1e6_4_Probe_3D_10_S4,
         ProbeA_THP1_1e6_5 = THP1_1e6_5_Probe_1_S5)

ProbeTest4_RawReads <- read.csv("Data/ProbeTest4/Mtb.Expression.Gene.Data.readsM.csv")
ProbeTest4_RawReads <- ProbeTest4_RawReads %>%
  select(X, THP1_1e6_1_S41, THP1_1e6_2_S42, THP1_1e6_3_S43, THP1_1e6_4_S44, THP1_1e6_5_S45) %>%
  rename(ProbeB_THP1_1e6_1 = THP1_1e6_1_S41,
         ProbeB_THP1_1e6_2 = THP1_1e6_2_S42,
         ProbeB_THP1_1e6_3 = THP1_1e6_3_S43,
         ProbeB_THP1_1e6_4 = THP1_1e6_4_S44,
         ProbeB_THP1_1e6_5 = THP1_1e6_5_S45)

# Merge the RawReads I collected above
All_RawReads <- merge(Run1_RawReads, ProbeTest3_RawReads, all = T)
All_RawReads <- merge(All_RawReads, ProbeTest5_RawReads_CapturedVsNot_Broth_LimitofDetect)
All_RawReads <- merge(All_RawReads, ProbeTest4_RawReads)

# Keep only Rv protein coding genes
All_RawReads_f <- All_RawReads %>% 
  filter(str_detect(X, "^Rv\\d+.*"))

RawReads_wb <- createWorkbook()
addWorksheet(RawReads_wb, "All_RawReads")
writeData(RawReads_wb, "All_RawReads", All_RawReads_f)
saveWorkbook(RawReads_wb, "RawReads.xlsx", overwrite = TRUE)
