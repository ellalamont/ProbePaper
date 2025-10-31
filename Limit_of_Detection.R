# Limit of Detection with the THP1 spiked samples graphs
# E. Lamont
# 7/31/25

# 9/8/25: Looking at correlation between the lower burden samples

source("Import_data.R") # to get LimitofDetect_pipeSummary


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
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
################ CELL NUMBER VS N_GENOMIC #################

# ggerrorplot
LimitofDetect_NumReads_Fig1 <- LimitofDetect_pipeSummary %>% 
  ggerrorplot(x = "Ra_cells2", y = "N_Genomic", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4,  # Size of error bars
              add.params = list(size = 0.4)) +  # Size of mean points
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # annotate("text", x = 5.6, y = 1000000*0.8, label = "1 million", 
  #          hjust = 1.1, vjust = -0.5, color = "black") + 
  # scale_y_continuous(limits = c(0,19000000), breaks = seq(0, 19000000, 2000000)) + 
  scale_y_continuous(limits = c(0,20000000), breaks = seq(0, 20000000, 4000000)) + 
  labs(# title = "ProbeTest5 THP1 cells spiked with H37Ra", 
       # subtitle = "Mean with standard deviation", 
       x = "# spiked in H37Ra cells", 
       y = "# reads aligning to Mtb") + 
  my_plot_themes
LimitofDetect_NumReads_Fig1
# ggsave(LimitofDetect_NumReads_Fig1,
#        file = "N_Genomic_Limit.of.Detection_fig1.pdf",
#        path = "Figures/Limit_of_Detection",
#        width = 7, height = 5, units = "in")


###########################################################
################ CELL NUMBER VS P_GENOMIC #################

# ggerrorplot
LimitofDetect_PercentReads_Fig1 <- LimitofDetect_pipeSummary %>% 
  ggerrorplot(x = "Ra_cells2", y = "P_Genomic", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4,  # Size of error bars
              add.params = list(size = 0.4)) +  # Size of mean points
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  labs(# title = "ProbeTest5 THP1 cells spiked with H37Ra", 
       # subtitle = "Mean with standard deviation", 
       x = "# spiked in H37Ra cells", 
       y = "% reads aligning to Mtb") + 
  my_plot_themes
LimitofDetect_PercentReads_Fig1
# ggsave(LimitofDetect_PercentReads_Fig1,
#        file = "P_Genomic_Limit.of.Detection_fig1.pdf",
#        path = "Figures/Limit_of_Detection",
#        width = 7, height = 5, units = "in")


###########################################################
############# CELL NUMBER VS AtLeast.10.Reads #############

# ggerrorplot
LimitofDetect_10Reads_Fig1 <- LimitofDetect_pipeSummary %>% 
  mutate(Txn_Coverage = (AtLeast.10.Reads/4499)*100) %>%
  ggerrorplot(x = "Ra_cells2", y = "Txn_Coverage", desc_stat = "mean_sd", error.plot = "errorbar", add = "mean", color = "black", size = 0.4,  # Size of error bars
              add.params = list(size = 0.4)) +  # Size of mean points
  geom_point(alpha = 0.7, position = position_jitter(width = 0.1, seed = 42), size = 1) + 
  geom_text_repel(aes(label = format(SampleID, big.mark = ",")), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  geom_hline(yintercept = 80, linetype = "dashed", alpha = 0.5) + 
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) + 
  labs(x = "# spiked in H37Ra cells", 
       y = "% transcriptional coverage") + 
  my_plot_themes
LimitofDetect_10Reads_Fig1
# ggsave(LimitofDetect_10Reads_Fig1,
#        file = "TenReads_Limit.of.Detection_fig1.pdf",
#        path = "Figures/Limit_of_Detection",
#        width = 7, height = 5, units = "in")


###########################################################
##################### CORRELATIONS ########################

# Start with the TPM ProbeTest5_tpm
# Grab just the samples I want
LimitofDetect_tpm <- ProbeTest5_tpm %>%
  select("X", all_of(LimitofDetect_pipeSummary$SampleID)) %>%
  rename("Gene" = "X")

names(LimitofDetect_tpm) <- gsub(x = names(LimitofDetect_tpm), pattern = "_S.*", replacement = "")

### Log10 transform ###
LimitofDetect_tpm_Log10 <- LimitofDetect_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

### Add averages columns ###
LimitofDetect_tpm_Log10 <- LimitofDetect_tpm_Log10 %>% 
  mutate(
    AVERAGE_1e6 = rowMeans(select(., c(THP1_1e6_1a, THP1_1e6_2b, THP1_1e6_3a)), na.rm = TRUE),
    AVERAGE_1e5 = rowMeans(select(., c(THP1_1e5_1, THP1_1e5_2, THP1_1e5_3)), na.rm = TRUE),
    AVERAGE_1e4 = rowMeans(select(., c(THP1_1e4_1, THP1_1e4_2, THP1_1e4_3)), na.rm = TRUE),
    AVERAGE_1e3 = rowMeans(select(., c(THP1_1e3_1, THP1_1e3_2, THP1_1e3_3, THP1_1e3_4)), na.rm = TRUE),
    AVERAGE_1e2 = rowMeans(select(., c(THP1_1e2_1, THP1_1e2_2, THP1_1e2_3, THP1_1e2_4)), na.rm = TRUE),
  )

###################################################################
############################## 1e2 VS 1e6 #########################

Sample1 <- "AVERAGE_1e2"
Sample2 <- "AVERAGE_1e6"
ScatterCorr <- LimitofDetect_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 0.5, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(# title = paste0("Samples AVERAGED: ", Sample1, " vs ", Sample2),
    # subtitle = "Pearson correlation; 1e6 Ra THP1 spiked captured VS Broth Not captured (Not scaled)",
    x = paste0("Log10(TPM+1)\n1e2 Mtb cells"),
    y = paste0("Log10(TPM+1)\n1e6 Mtb cells"), ) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
LimitofDetect_ScatterCorr <- ScatterCorr
# ggsave(ScatterCorr,
#        file = "1e2.vs.1e6_Correlation.pdf",
#        path = "Figures_preNonCodingRemoval/Limit_of_Detection",
#        width = 7, height = 5, units = "in")


###################################################################
###################### COMPARE SINGLE SAMPLES #####################

LimitofDetect_pipeSummary <- LimitofDetect_pipeSummary %>% mutate(Coverage = round(AtLeast.10.Reads/4499*100))

Sample1 <- "THP1_1e4_3" # 4_1, 4_2, 4_3, 3_3, 3_2, 2_3
Sample2 <- "THP1_1e6_1a" # 3a, 1a, 2b
ScatterCorr <- LimitofDetect_tpm_Log10 %>% 
  # select(Gene, .data[[Sample1]], .data[[Sample2]]) %>%
  # filter(if_all(everything(), ~ . >= 1)) %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 0.5, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(x = paste0("Log10(TPM+1)\n", Sample1, " ", LimitofDetect_pipeSummary %>% filter(str_detect(SampleID, Sample1)) %>% select(Coverage) %>% pull(), "% coverage"),
    y = paste0("Log10(TPM+1)\n", Sample2, " ", LimitofDetect_pipeSummary %>% filter(str_detect(SampleID, Sample2)) %>% select(Coverage) %>% pull(), "% coverage"), ) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0(Sample1, ".vs.", Sample2, ".pdf"),
       path = "Figures_preNonCodingRemoval/Limit_of_Detection/InterestingCorrelations",
       width = 7, height = 5, units = "in")


###################################################################
########## LOOP THROUGH SINGLE SAMPLE CORRELATION PLOTS ###########

# Not quite doing all of them...

# my_path <- "Figures_preNonCodingRemoval/Limit_of_Detection/AllCorrelations"
# 
# for (i in 2:(length(colnames(LimitofDetect_tpm_Log10)) -1 -1)) {
#   for (j in (i + 1):(length(colnames(LimitofDetect_tpm_Log10))-1)) {
#     if (i != j) { # Avoid comparing the same samples or repeating comparisons
#       # Access the samples
#       Sample1 <- colnames(LimitofDetect_tpm_Log10[i])
#       Sample2 <- colnames(LimitofDetect_tpm_Log10[j])
#       # cat("Comparing:", Sample1, "with", Sample2, "\n")
#       filename <- paste0(Sample1, "_ComparedTo_", Sample2, ".pdf")
#       
#       ScatterCorr <- LimitofDetect_tpm_Log10 %>%
#         ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) +
#         geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
#         geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
#         labs(title = paste0(Sample1, " vs ", Sample2, " TPM (4499 genes)"),
#              subtitle = "Pearson correlation",
#              x = paste0(Sample1, " Log10(TPM)"), y = paste0(Sample2, " Log10(TPM)")) +
#         stat_cor(method="pearson") + # add a correlation to the plot
#         my_plot_themes
#       
#       ggsave(ScatterCorr,
#              file = filename,
#              path = my_path,
#              width = 7, height = 5, units = "in")
#     }
#   }
# }





###########################################################
################# PEARSON LOG10 GGCORRPLOT ################

# Make the correlation
LimitofDetect_tpm_Log10_AVG <- LimitofDetect_tpm_Log10 %>% 
  select("Gene", contains("AVERAGE")) %>%
  column_to_rownames("Gene")
corr <- cor(LimitofDetect_tpm_Log10_AVG, method = "pearson")

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(LimitofDetect_tpm_Log10_AVG)
# head(p.mat[, 1:4])

min(corr) # 0.1863461
my_min <- round(min(corr), 1)

# Plot pearson
ggcorrplot_PearsonLog10 <- corr %>% 
  ggcorrplot(hc.order = F, 
             method = "square", 
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  # scale_fill_gradient2(limit = c(my_min,1), low = "blue", high =  "red", mid = "white", midpoint = (((1-my_min)/2)+my_min)) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Pearson Correlation Log10(TPM+1)", 
       subtitle = NULL, 
       fill = "Correlation")
ggcorrplot_PearsonLog10
# ggsave(ggcorrplot_PearsonLog10,
#        file = "LimitofDetect_PearsonCorrelation.pdf",
#        path = "Figures_preNonCodingRemoval/Limit_of_Detection",
#        width = 7, height = 5, units = "in")
