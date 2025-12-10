# Import Literature Data
# 4/16/25

# Want to import the literature data and process it in some way so I can compare....
# Maybe something like Coppola et al. (2021)? https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.763364/full#h13

# 4/25/25
# For DEG, set the log2fold thershold to 2.5, because that's what Garton2008 did. Honeyborne2016 is at 2. But I don't think I can do less than 2.5 because thats all I have for Garton2008


#############################################################
################ GENES ALL FROM COPPOLA 2021 ################
# Just take from the supplemental of Coppola2021 for all so it will be more consistent! 

Coppola2021_all <- read.csv("Data/LiteratureData/Coppola2021_EL.csv")
Coppola2021_Walter2015 <- Coppola2021_all %>% 
  select(Walter2015_Gene, Walter2015_MedianRelativeScoreRank) %>% 
  drop_na()
Coppola2021_Garcia2016 <- Coppola2021_all %>% 
  select(Garcia2016_Gene, Garcia2016_PercentileRelativeScoreRank) %>% 
  drop_na()
Coppola2021_Sharma2017 <- Coppola2021_all %>% 
  select(Sharma2017_Gene, Sharma2017_MedianRelativeScoreRank) %>% 
  drop_na()
Coppola2021_Lai2021 <- Coppola2021_all %>% 
  select(Lai2021_Gene, Lai2021_MedianRelativeScoreRank) %>% 
  drop_na()

###########################################################
############# PUT EVERTHING IN ONE DATAFRAME ##############

combined_lit_df <- full_join(Coppola2021_Walter2015 %>% rename(Gene = Walter2015_Gene),
                         Coppola2021_Garcia2016 %>% rename(Gene = Garcia2016_Gene),
                         by = "Gene")
combined_lit_df <- full_join(combined_lit_df,
                         Coppola2021_Sharma2017 %>% rename(Gene = Sharma2017_Gene),
                         by = "Gene")
combined_lit_df <- full_join(combined_lit_df,
                         Coppola2021_Lai2021 %>% rename(Gene = Lai2021_Gene),
                         by = "Gene")


###########################################################
##################### SAVE FOR PAPER ######################
# Adding to the wb

addWorksheet(wb, "Literature.Data")
writeData(wb, "Literature.Data", combined_lit_df)
saveWorkbook(wb, "DEG.xlsx", overwrite = TRUE)



###########################################################
###################### LAI2021 TPM ########################
# 12/08/25

Lai2021_RawReads <- read.csv("Data/LiteratureData/Lai2021_BobsPipeline/Mtb.Expression.Gene.Data.ReadsM_Lai2021.csv")

source("Function_CalculateTPM.R")

# Keep only the Rv#* coding genes 
Lai2021_RawReads_f <- Lai2021_RawReads %>% 
  filter(str_detect(X, "^Rv\\d+.*"))
# Convert all to TPM
Lai2021_tpm_f <- CalculateTPM_RvOnly(Lai2021_RawReads_f)

