# Do a batch correction and new PCA with Lance's Ra
# 8/31/25
# E. Lamont


source("Import_data.R") # To get GoodBiolSamples_wLance_RawReads
source("TPM_from_RawReads.R") # To get the calculate TPM function


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

###########################################################
###################### PROCESS DATA #######################

pipeSummary2 <- BiolSamples_pipeSummary_2

# Add Batch designation to the metadata:
pipeSummary2 <- pipeSummary2 %>%
  mutate(Batch = case_when(
    Run == "ProbeTest3" ~ 1,
    Run == "ProbeTest5" ~ 2,
    Run == "PredictTB_Run1" ~ 3,
    Run == "LanceRun" ~ 4
  ))

All_RawReads_2 <- GoodBiolSamples_wLance_RawReads %>%
  column_to_rownames(var = "X") %>%
  as.matrix()

# Need to make sure the order of the SampleIDs matches
pipeSummary2 <- pipeSummary2 %>%
  arrange(match(SampleID, colnames(All_RawReads_2))) %>%
  slice(1:(n() - 4)) # Remove the last 4 rows which are not passing filter


###########################################################
#################### BATCH CORRECTION #####################
# This basically all taken from Mark and I'm not sure what's happening
# https://academic.oup.com/nargab/article/2/3/lqaa078/5909519?login=true

batch <- pipeSummary2$Batch # Check this and make sure it is in the correct order!!
counts_corrected <- ComBat_seq(All_RawReads_2, batch = batch)


###########################################################
######################## PCA w CPM ########################
# Just do the PCA the way I have been doing it

my_fav_colors3 <- c(`Sputum L4` = "#004C73", `Sputum L2` = "#99CCE8", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00", `Rv7` = "#FC8D59", `Rv8.3` = "#FEE090", `Rv5.7` = "#D73027") 
# Labelled Shapes
my_fav_shapes3 <- c(`Sputum L4` = 21, `Sputum L2` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25, `Rv7` = 23, `Rv8.3` = 23, `Rv5.7` = 23)

# Convert the batch corrected counts to counts per million (cpm)
my_cpm <- cpm(counts_corrected)

# transform the data 
my_cpm_t <- as.data.frame(t(my_cpm))

# Remove columns that are all zero so the scale works for prcomp
my_cpm_t <- my_cpm_t %>% select_if(colSums(.) != 0) # Down to 4494 genes 

# Make the actual PCA
my_PCA_cpm <- prcomp(my_cpm_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA_cpm)
summary_PCA_cpm <- format(round(as.data.frame(summary(my_PCA_cpm)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA_cpm[1,1] # PC1 explains 27.5% of variance
summary_PCA_cpm[2,1] # PC2 explains 11.5% of variance
summary_PCA_cpm[3,1] # PC3 explains 10.4% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_cpm_df <- as.data.frame(my_PCA_cpm$x[, 1:3]) # Extract the first 3 PCs
my_PCA_cpm_df <- data.frame(SampleID = row.names(my_PCA_cpm_df), my_PCA_cpm_df)
my_PCA_cpm_df <- merge(my_PCA_cpm_df, BiolSamples_pipeSummary_2, by = "SampleID", )

PCA_BatchCorrected <- my_PCA_cpm_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type3, shape = Type3)) + 
  geom_point(aes(fill = Type3, shape = Type3), size = 5, alpha = 0.7, stroke = 0.8) +
  # geom_text_repel(aes(label = SampleID), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors3) +  
  scale_shape_manual(values = my_fav_shapes3) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "Batch Corrected CPM",
       x = paste0("PC1: ", summary_PCA_cpm[1,1], "%"),
       y = paste0("PC2: ", summary_PCA_cpm[2,1], "%")) +
  my_plot_themes
PCA_BatchCorrected
# ggsave(PCA_BatchCorrected,
#        file = paste0("BatchCorrected_CPM_wRv.pdf"),
#        path = "Figures_preNonCodingRemoval/PCA/wRv",
#        width = 8, height = 5, units = "in")


###########################################################
######################## PCA w TPM ########################
# Just do the PCA the way I have been doing it

#### NEED TO CHECK THE TPM CONVERSION HAPPENING HERE!!!

my_fav_colors3 <- c(`Sputum L4` = "#004C73", `Sputum L2` = "#99CCE8", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00", `Rv7` = "#FC8D59", `Rv8.3` = "#FEE090", `Rv5.7` = "#D73027") 
# Labelled Shapes
my_fav_shapes3 <- c(`Sputum L4` = 21, `Sputum L2` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25, `Rv7` = 23, `Rv8.3` = 23, `Rv5.7` = 23)

# Convert the batch corrected counts to counts per million (cpm)
counts_corrected_df <- as.data.frame(counts_corrected) %>% rownames_to_column("X")
my_tpm_bc <- CalculateTPM(counts_corrected_df)

# transform the data 
my_tpm_bc_t <- as.data.frame(t(my_tpm_bc))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_bc_t <- my_tpm_bc_t %>% select_if(colSums(.) != 0) # Down to 4494 genes 

# Make the actual PCA
my_PCA_tpm_bc <- prcomp(my_tpm_bc_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA_tpm_bc)
summary_PCA_tpm_bc <- format(round(as.data.frame(summary(my_PCA_tpm)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA_tpm_bc[1,1] # PC1 explains 32.0% of variance
summary_PCA_tpm_bc[2,1] # PC2 explains 13.8% of variance
summary_PCA_tpm_bc[3,1] # PC3 explains 8.5% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_tpm_bc_df <- as.data.frame(my_PCA_tpm_bc$x[, 1:3]) # Extract the first 3 PCs
my_PCA_tpm_bc_df <- data.frame(SampleID = row.names(my_PCA_tpm_bc_df), my_PCA_tpm_bc_df)
my_PCA_tpm_bc_df <- merge(my_PCA_tpm_bc_df, BiolSamples_pipeSummary_2, by = "SampleID", )

PCA_BatchCorrected <- my_PCA_tpm_bc_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type3, shape = Type3)) + 
  geom_point(aes(fill = Type3, shape = Type3), size = 5, alpha = 0.7, stroke = 0.8) +
  # geom_text_repel(aes(label = SampleID), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors3) +  
  scale_shape_manual(values = my_fav_shapes3) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "Batch Corrected Converted to TPM",
       x = paste0("PC1: ", summary_PCA_tpm_bc[1,1], "%"),
       y = paste0("PC2: ", summary_PCA_tpm_bc[2,1], "%")) +
  my_plot_themes
PCA_BatchCorrected



###########################################################
############## PCA w CPM NOT BATCH CORRECTED ##############
# Just do the PCA the way I have been doing it

my_fav_colors3 <- c(`Sputum L4` = "#004C73", `Sputum L2` = "#99CCE8", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00", `Rv7` = "#FC8D59", `Rv8.3` = "#FEE090", `Rv5.7` = "#D73027") 
# Labelled Shapes
my_fav_shapes3 <- c(`Sputum L4` = 21, `Sputum L2` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25, `Rv7` = 23, `Rv8.3` = 23, `Rv5.7` = 23)

# Convert the batch corrected counts to counts per million (cpm)
my_cpm <- cpm(All_RawReads_2)

# transform the data 
my_cpm_t <- as.data.frame(t(my_cpm))

# Remove columns that are all zero so the scale works for prcomp
my_cpm_t <- my_cpm_t %>% select_if(colSums(.) != 0) # Down to 4494 genes 

# Make the actual PCA
my_PCA_cpm <- prcomp(my_cpm_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA_cpm)
summary_PCA_cpm <- format(round(as.data.frame(summary(my_PCA_cpm)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA_cpm[1,1] # PC1 explains 29.9% of variance
summary_PCA_cpm[2,1] # PC2 explains 14.0% of variance
summary_PCA_cpm[3,1] # PC3 explains 9.0% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_cpm_df <- as.data.frame(my_PCA_cpm$x[, 1:3]) # Extract the first 3 PCs
my_PCA_cpm_df <- data.frame(SampleID = row.names(my_PCA_cpm_df), my_PCA_cpm_df)
my_PCA_cpm_df <- merge(my_PCA_cpm_df, BiolSamples_pipeSummary_2, by = "SampleID", )

PCA_cpm <- my_PCA_cpm_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type3, shape = Type3)) + 
  geom_point(aes(fill = Type3, shape = Type3), size = 5, alpha = 0.7, stroke = 0.8) +
  # geom_text_repel(aes(label = SampleID), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors3) +  
  scale_shape_manual(values = my_fav_shapes3) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "NOT Batch Corrected CPM",
       x = paste0("PC1: ", summary_PCA_cpm[1,1], "%"),
       y = paste0("PC2: ", summary_PCA_cpm[2,1], "%")) +
  my_plot_themes
PCA_cpm
# Doesn't look like CPM is causing much difference, so what is going on with the TPM converted batch corrected samples???
