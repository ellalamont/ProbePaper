# Trying to combine figures so they are all consistent
# E. Lamont
# 8/8/25

library(cowplot)

###########################################################
######################## FIGURE 2 #########################

source("Captured.vs.Not_Correlations.R") # ScatterCorr_2D, ScatterCorr_2E
source("Comparing_across_Runs.R") # ScatterCorr_2F
source("Captured.vs.Not_ReadCounts.R") # PieChart_Averages_fig1, CapturedVsNot_N.Genomic_fig1, CapturedvsNot_10Reads_Fig1
source("Limit_of_Detection.R") # LimitofDetect_NumReads_Fig1, LimitofDetect_PercentReads_Fig1, LimitofDetect_10Reads_Fig1, LimitofDetect_ScatterCorr

# combined <- plot_grid(PieChart_Averages_fig1, CapturedVsNot_N.Genomic_fig1, CapturedvsNot_10Reads_Fig1, ScatterCorr_2D, ScatterCorr_2E, ScatterCorr_2F, LimitofDetect_NumReads_Fig1, LimitofDetect_PercentReads_Fig1, LimitofDetect_10Reads_Fig1, ncol = 3, nrow = 3, axis = "tblr", align = "hv", rel_widths = c(1, 1, 1), rel_heights = c(1, 1))

combined <- plot_grid(PieChart_Averages_fig1, CapturedVsNot_N.Genomic_fig1, CapturedvsNot_10Reads_Fig1,
                      LimitofDetect_NumReads_Fig1, LimitofDetect_PercentReads_Fig1, LimitofDetect_10Reads_Fig1,
                      LimitofDetect_ScatterCorr, ScatterCorr_2E, ScatterCorr_2F, 
                      ncol = 3, nrow = 3, axis = "tblr", align = "hv", rel_widths = c(1, 1, 1), rel_heights = c(1, 1))

combined
# ggsave(combined,
#        file = paste0("Figure2_v2.pdf"),
#        path = "Figures/CombinedFigures",
#        width = 15, height = 10, units = "in")
# ggsave(combined,
#        file = paste0("Figure2_v3.pdf"),
#        path = "Figures/CombinedFigures",
#        width = 15, height = 12, units = "in")
ggsave(combined,
       file = paste0("Figure2_v5.pdf"),
       path = "Figures_preNonCodingRemoval/CombinedFigures",
       width = 15, height = 12, units = "in")

###########################################################
######################## FIGURE 3 #########################

source("BiolSamples_GenomicReads.R")
source("Reads_w_Metadata.R")

combined4 <- plot_grid(N_Genomic_box1, P_Genomic_box1, TenReads_box1, CFU.Percent, ctVsPercent, ttdVsPercent, align = "hv", axis = "tblr", nrow = 2, ncol = 3, rel_widths = c(1, 1, 1), rel_heights = c(1, 1))
combined4
ggsave(combined4,
       file = paste0("Figure3_v3.pdf"),
       path = "Figures/CombinedFigures",
       width = 15, height = 10, units = "in")
# ggsave(combined4,
#        file = paste0("Figure3_v1.png"),
#        path = "Figures/CombinedFigures",
#        width = 15, height = 10, units = "in")


