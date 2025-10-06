# Probe Prep testing 5
# Probes made by Jessica, run included THP1 test samples and sputum
# See Test_Sputum_Jan2025_RNALibraryPrep
# 3/11/25

# Importing data on SNPs for the phoP gene (Rv0757)
# H37Ra should have a SNP at this gene

# App is showing:
# Mutation from C -> T at 852263
# Since it starts at 851608, this is at position 655
# "one single nucleotide polymorphism affecting codon 219 (TCG → TTG) in phoP and changing a serine to a leucine. S219 in PhoP of M. tuberculosis is equivalent to R200 in PhoB of E. coli" - https://pmc.ncbi.nlm.nih.gov/articles/PMC2238218/


###########################################################
############### IMPORT phoP SEQUENCE DATA #################

SNPs_df <- read.csv("Data/ProbeTest5/EL_phoP_All.Base.Depth.Calls.csv")

levels(as.factor(SNPs_df$SampleID))
# This only contains the 3 1e6 H37Ra spiked THP1 samples that were captured

SNPs_df$SampleID <- gsub("_S(.*)", "", as.character(SNPs_df$SampleID))

# Save in the wb
# addWorksheet(wb, "phoP.SNPs")
# writeData(wb, "phoP.SNPs", SNPs_df)
# saveWorkbook(wb, "DEG.xlsx", overwrite = TRUE)

###########################################################
################# ISOLATE RRDR REGION #####################

# Will do this if need to to remove ends

# class(df_3$POSITION)
# df_4 <- df_3 %>% filter(between(POSITION, 761085, 761171))


###########################################################
######################## GATHER ###########################

# Change the name of the column containing the number of reads aligning to reference from X. to REF
names(SNPs_df)[names(SNPs_df) == "X."] <- "REF"

# Gather (pivot_longer) ####
SNPs_df_2 <- SNPs_df %>% pivot_longer(cols = c("A","C", "G", "T", "N", "Indel"),
                              names_to = "base_call",
                              values_to = "raw_read_number")

# Right now, the reads for the reference are stored in column X. and the base this is just has 0 for the raw_read_number

# Need to give the Ref base position the correct read numbers!

# Using mutate and case_when: case_when(condition ~ output-value)
# Condition must be true or false. The TRUE in the second line forces case_when to output the else-output-value if none of the previous conditions were true
# If REF_BASE is the same as base call, then change raw_read_number to be the value of REF, if this is false, then go to the next line where everything is TRUE and raw_read_number is changed to just be raw_read_number again
# This way is better than ifelse when there are multiple conditions 
SNPs_df_3 <- SNPs_df_2 %>% mutate(raw_read_number = case_when(REF_BASE == base_call ~ REF, TRUE ~ raw_read_number))

# Change correct base call to be "genome" ####
# Doing this so it can all plot as one color on the graph

SNPs_df_4 <- SNPs_df_3 %>% mutate(base_call = case_when(REF_BASE == base_call ~ "Genome", TRUE ~ base_call))

# Order the values ####
ordered_bases <- c("Genome", "A", "C", "G", "T", "N", "Indel")
SNPs_df_4$base_call <- factor(SNPs_df_4$base_call, levels = ordered_bases)


###########################################################
############### ADD PROPORTION COLUMN #####################

SNPs_df_4 <- SNPs_df_4 %>% mutate(Proportion = raw_read_number / DEPTH)


