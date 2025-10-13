

###################################################################
############################## SETUP ##############################

library( DuffyNGS)
nCores <- multicore.setup( max.cores=6)
annT <- readAnnotationTable( "Annotation.txt", sep = "\t")
anot = read.delim('Annotation.txt',as.is=T, sep = "\t")
setCurrentSpecies( "MT_H37")


###################################################################
########################## ALIGNMENT ##############################

for (a in anot$SampleID){
  pipeline(a, optionsFile = "Options.txt")
}

s = extractPipelineSummaryDetails(anot$SampleID, optionsFile = "Options.txt")
write.csv(s,'Pipeline.Summary.Details.csv')

sids <- anot$SampleID
fset <- file.path( "results/transcript/", paste( sids, "MTb.Transcript.txt",sep="."))

n <- expressionFileSetToMatrix(fset, sids, intensityColumn="READS_M")
write.csv(n, 'Mtb.Expression.Gene.Data.ReadsM.csv')

cleanupBAMfiles()


###################################################################
############# DIFFERENTIAL GENE EXPRESSION - SETUP ################

# Samples passing filter
GoodSampleList <- c("BQ12_10_Probe_3A_S29", "HN878_mimic_D14_R1_S63", "HN878_mimic_D14_R2_S64", "BQ12_8_Probe_4A_50_S28", "HN878_mimic_D28_R1_S65",  "H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "HN878_mimic_D7_R2_S62", "H37Ra_Broth_6_S9", "LLL_Cav_L2_18_S58", "LLL_Cav_L2_21_S59", "LU_Cav_L1_12_S57", "W0_11011_S16", "W0_12007_S2", "W0_12008_S3", "W0_12010_S5", "W0_12032_S28", "W0_12083_S39", "W0_13045_S34", "W0_15081_S49", "W0_15089_S51")

# Keep just the samples that are in my GoodSampleList
annT_GoodSamples <- annT %>% filter(SampleID %in% GoodSampleList)

# Change the group and color for the new annotation file
annT_GoodSamples <- annT_GoodSamples %>% 
  mutate(Group = case_when(SampleID %in% c("HN878_mimic_D14_R1_S63", "HN878_mimic_D14_R2_S64", "HN878_mimic_D28_R1_S65",  "HN878_mimic_D7_R2_S62") ~ "CaseumMimic", TRUE ~ Group)) %>%
  mutate(Group = case_when(SampleID %in% c("BQ12_10_Probe_3A_S29", "BQ12_8_Probe_4A_50_S28") ~ "Marmoset", TRUE ~ Group)) %>%
  mutate(Group = case_when(SampleID %in% c("LLL_Cav_L2_18_S58", "LLL_Cav_L2_21_S59", "LU_Cav_L1_12_S57") ~ "Rabbit", TRUE ~ Group))

annT_GoodSamples <- annT_GoodSamples %>% 
  mutate(Color = case_when(SampleID %in% c("HN878_mimic_D14_R1_S63", "HN878_mimic_D14_R2_S64", "HN878_mimic_D28_R1_S65",  "HN878_mimic_D7_R2_S62") ~ "green", TRUE ~ Color)) %>%
  mutate(Color = case_when(SampleID %in% c("BQ12_10_Probe_3A_S29", "BQ12_8_Probe_4A_50_S28") ~ "purple", TRUE ~ Color)) %>%
  mutate(Color = case_when(SampleID %in% c("LLL_Cav_L2_18_S58", "LLL_Cav_L2_21_S59", "LU_Cav_L1_12_S57") ~ "orange", TRUE ~ Color)) %>%
  mutate(Color = case_when(SampleID %in% c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9") ~ "black", TRUE ~ Color))


# Need to re-save this because the code below finds the actual txt file
write.table(annT_GoodSamples, "Annotation_GoodSampleList.txt", sep = "\t", row.names = F, quote = F)

my_results.path <- "results"

###################################################################
###################### COMPARE SPUTUM L2 to L4 ####################
L2_Names <- c("W0_12083_S39", "W0_13045_S34", "W0_15081_S49", "W0_15089_S51")
L4_Names <- c("W0_11011_S16", "W0_12007_S2", "W0_12008_S3", "W0_12010_S5", "W0_12032_S28")
pipe.MetaResults(c(L2_Names, L4_Names), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "MetaResults.GoodSamples_L2.vs.L4", 
                 groupColumn ="Group2",
                 colorColumn ="Color", 
                 results.path = my_results.path, 
                 PLOT.FUN=NA)


###################################################################
###################### COMPARE SPUTUM TO BROTH ####################
SputumNames <- c("W0_11011_S16", "W0_12007_S2", "W0_12008_S3", "W0_12010_S5", "W0_12032_S28", "W0_12083_S39", "W0_13045_S34", "W0_15081_S49", "W0_15089_S51")
BrothNames <- c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9")
pipe.MetaResults(c(SputumNames, BrothNames), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "iModulons.GoodSamples_Sputum.vs.Broth", 
                 groupColumn ="Group",
                 colorColumn ="Color", 
                 results.path = my_results.path, 
                 PLOT.FUN=NA)
pipe.MetaGeneSets(c(SputumNames, BrothNames), 
                  annotationFile = "Annotation_GoodSampleList.txt",
                  folderName = "iModulons.GoodSamples_Sputum.vs.Broth", 
                  groupColumn="Group", 
                  colorColumn="Color", 
                  geneSets = defaultGeneSets(speciesID = "MT_H37")[c(18)], # This is iModulons
                  results.path= my_results.path)

###################################################################
################### COMPARE CASEUM MIMIC TO BROTH #################
MimicNames <- c("HN878_mimic_D14_R1_S63", "HN878_mimic_D14_R2_S64", "HN878_mimic_D28_R1_S65",  "HN878_mimic_D7_R2_S62")
BrothNames <- c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9")
pipe.MetaResults(c(MimicNames, BrothNames), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "MetaResults.GoodSamples_CaseumMimic.vs.Broth", 
                 groupColumn ="Group", 
                 colorColumn ="Color", 
                 results.path = my_results.path,
                 PLOT.FUN=NA)

###################################################################
##################### COMPARE RABBIT TO BROTH #####################
RabbitNames <- c("LLL_Cav_L2_18_S58", "LLL_Cav_L2_21_S59", "LU_Cav_L1_12_S57")
BrothNames <- c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9")
pipe.MetaResults(c(RabbitNames, BrothNames), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "MetaResults.GoodSamples_Rabbit.vs.Broth", 
                 groupColumn ="Group",
                 colorColumn ="Color", 
                 results.path = my_results.path,
                 PLOT.FUN=NA)

###################################################################
#################### COMPARE MARMOSET TO BROTH ####################
MarmosetNames <- c("BQ12_10_Probe_3A_S29", "BQ12_8_Probe_4A_50_S28")
BrothNames <- c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9")
pipe.MetaResults(c(MarmosetNames, BrothNames), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "MetaResults.GoodSamples_Marmoset.vs.Broth", 
                 groupColumn ="Group", 
                 colorColumn ="Color", 
                 results.path = my_results.path,
                 PLOT.FUN=NA)

###################################################################
#################### COMPARE SPUTUM TO MARMOSET ###################
SputumNames <- c("W0_11011_S16", "W0_12007_S2", "W0_12008_S3", "W0_12010_S5", "W0_12032_S28", "W0_12083_S39", "W0_13045_S34", "W0_15081_S49", "W0_15089_S51")
MarmosetNames <- c("BQ12_10_Probe_3A_S29", "BQ12_8_Probe_4A_50_S28")
pipe.MetaResults(c(SputumNames, MarmosetNames), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "MetaResults.GoodSamples_Sputum.vs.Marmoset", 
                 groupColumn ="Group",
                 colorColumn ="Color", 
                 results.path = my_results.path,
                 PLOT.FUN=NA)

###################################################################
#################### COMPARE SPUTUM TO RABBIT #####################
SputumNames <- c("W0_11011_S16", "W0_12007_S2", "W0_12008_S3", "W0_12010_S5", "W0_12032_S28", "W0_12083_S39", "W0_13045_S34", "W0_15081_S49", "W0_15089_S51")
RabbitNames <- c("LLL_Cav_L2_18_S58", "LLL_Cav_L2_21_S59", "LU_Cav_L1_12_S57")
pipe.MetaResults(c(SputumNames, RabbitNames), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "MetaResults.GoodSamples_Sputum.vs.Rabbit", 
                 groupColumn ="Group", 
                 colorColumn ="Color", 
                 results.path = my_results.path, 
                 PLOT.FUN=NA)

###################################################################
################## COMPARE SPUTUM TO CASEUM MIMIC #################
SputumNames <- c("W0_11011_S16", "W0_12007_S2", "W0_12008_S3", "W0_12010_S5", "W0_12032_S28", "W0_12083_S39", "W0_13045_S34", "W0_15081_S49", "W0_15089_S51")
MimicNames <- c("HN878_mimic_D14_R1_S63", "HN878_mimic_D14_R2_S64", "HN878_mimic_D28_R1_S65",  "HN878_mimic_D7_R2_S62")
pipe.MetaResults(c(SputumNames, MimicNames), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "MetaResults.GoodSamples_Sputum.vs.CaseumMimic", 
                 groupColumn ="Group",
                 colorColumn ="Color", 
                 results.path = my_results.path,
                 PLOT.FUN=NA)

###################################################################
################ COMPARE CASEUM MIMIC TO MARMOSET #################
MimicNames <- c("HN878_mimic_D14_R1_S63", "HN878_mimic_D14_R2_S64", "HN878_mimic_D28_R1_S65",  "HN878_mimic_D7_R2_S62")
MarmosetNames <- c("BQ12_10_Probe_3A_S29", "BQ12_8_Probe_4A_50_S28")
pipe.MetaResults(c(MimicNames, MarmosetNames), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "MetaResults.GoodSamples_CaseumMimic.vs.Marmoset", 
                 groupColumn ="Group", 
                 colorColumn ="Color", 
                 results.path = my_results.path, 
                 PLOT.FUN=NA)

###################################################################
################## COMPARE CASEUM MIMIC TO RABBIT #################
MimicNames <- c("HN878_mimic_D14_R1_S63", "HN878_mimic_D14_R2_S64", "HN878_mimic_D28_R1_S65",  "HN878_mimic_D7_R2_S62")
RabbitNames <- c("LLL_Cav_L2_18_S58", "LLL_Cav_L2_21_S59", "LU_Cav_L1_12_S57")
pipe.MetaResults(c(MimicNames, RabbitNames), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "MetaResults.GoodSamples_CaseumMimic.vs.Rabbit", 
                 groupColumn ="Group",
                 colorColumn ="Color", 
                 results.path = my_results.path, 
                 PLOT.FUN=NA)

###################################################################
##################### COMPARE RABBIT TO MARMOSET ##################
RabbitNames <- c("LLL_Cav_L2_18_S58", "LLL_Cav_L2_21_S59", "LU_Cav_L1_12_S57")
MarmosetNames <- c("BQ12_10_Probe_3A_S29", "BQ12_8_Probe_4A_50_S28")
pipe.MetaResults(c(RabbitNames, MarmosetNames), 
                 annotationFile = "Annotation_GoodSampleList.txt",
                 folderName = "MetaResults.GoodSamples_Rabbit.vs.Marmoset", 
                 groupColumn ="Group", 
                 colorColumn ="Color", 
                 results.path = my_results.path, 
                 PLOT.FUN=NA)

