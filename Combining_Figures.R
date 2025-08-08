# Trying to combine figures so they are all consistent
# E. Lamont
# 8/8/25

library(cowplot)

###########################################################
######################## FIGURE 2 #########################

source("Captured.vs.Not_Correlations.R")
source("Comparing_across_Runs.R")
source("Captured.vs.Not_ReadCounts.R")

combined <- plot_grid(PieChart_Averages_fig1, ScatterCorr_2D, CapturedVsNot_N.Genomic_fig1, ScatterCorr_2E, CapturedvsNot_10Reads_Fig1, ScatterCorr_2F, align = "h", axis = "tblr", ncol = 2, nrow = 3)
combined

###########################################################
######################## FIGURE 3 #########################

source("BiolSamples_GenomicReads.R")

combined2 <- plot_grid(N_Genomic_box1, P_Genomic_box1, TenReads_box1, align = "h", axis = "tblr", labels = "AUTO", nrow = 1)
combined2

source("Reads_w_Metadata.R")

combined3 <- plot_grid(CFU.Percent, ctVsPercent, ttdVsPercent, axis = "tblr", labels = "AUTO", nrow = 1)
combined3

combined4 <- plot_grid(N_Genomic_box1, P_Genomic_box1, TenReads_box1, CFU.Percent, ctVsPercent, ttdVsPercent, align = "hv", axis = "tblr", nrow = 2, ncol = 3, rel_widths = c(1, 1, 1), rel_heights = c(1, 1))
combined4
ggsave(combined4,
       file = paste0("Figure3_v1.pdf"),
       path = "Figures/CombinedFigures",
       width = 15, height = 10, units = "in")
ggsave(combined4,
       file = paste0("Figure3_v1.png"),
       path = "Figures/CombinedFigures",
       width = 15, height = 10, units = "in")


