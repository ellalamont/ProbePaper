
# Can I get all the way through the way people normally do?


# https://github.com/JulioLeonIncio/Tutorial-from-DEseq2-to-GSEA-in-R/blob/main/DEseq2_to_GSEA_JL.Rmd
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("DESeq2","sva","fgsea","clusterProfiler","GSEABase","tidyverse","pheatmap")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask=FALSE)
}
library(DESeq2)
library(sva) # ComBat_seq
library(fgsea)
library(clusterProfiler)
library(GSEABase)
library(tidyverse)
library(pheatmap)
library(fgsea)
library(dplyr)


# Start with GoodBiolSamples_wRv_RawReads
my_RawCounts <- GoodBiolSamples_wRv_RawReads %>% # rename("GENE_ID" = "X")
  column_to_rownames("X") %>%
  mutate_if(is.numeric, ~ as.integer(round(., 0)))

# Start with BiolSamples_pipeSummary_2
my_metadata <- BiolSamples_pipeSummary_2 %>% 
  filter(SampleID != "INDIGO.NoDrug.Control.E1") %>% 
  filter(N_Genomic >= 1000000 & AtLeast.10.Reads >= (4499*0.8) | Run == "LanceRun") %>% # To only keep the good sputum
  mutate(batch = case_when(
    Run == "ProbeTest3" ~ 1,
    Run == "ProbeTest5" ~ 2,
    Run == "PredictTB_Run1" ~ 3,
    Run == "LanceRun" ~ 4)) %>%
  dplyr::select(SampleID, Type3, batch) %>%
  mutate(Type3 = if_else(Type3 %in% c("Sputum L4", "Sputum L2"), "Sputum", Type3)) %>% 
  rename("sample" = "SampleID", "condition" = "Type3") %>%
  column_to_rownames("sample")

# Ensure sample names line up
stopifnot(all(colnames(my_RawCounts) %in% rownames(my_metadata)))
my_metadata <- my_metadata[colnames(my_RawCounts), ]
stopifnot(all(colnames(my_RawCounts) == rownames(my_metadata)))
# Basic sanity
head(my_metadata)
summary(colSums(my_RawCounts))    # library sizes


# PCA to visualize data before correction
summary(colSums(my_RawCounts))                   # total reads per sample
hist(log10(colSums(my_RawCounts)+1))             # distribution of library sizes

dds_raw <- DESeqDataSetFromMatrix(countData = my_RawCounts,
                                  colData = my_metadata,
                                  design = ~ condition)

dds_raw <- dds_raw[rowSums(counts(dds_raw)) > 10, ]
vsd_raw <- vst(dds_raw, blind = TRUE)

plotPCA(vsd_raw, intgroup="condition") + ggtitle("PCA by condition (raw)")
plotPCA(vsd_raw, intgroup="batch") + ggtitle("PCA by batch (raw)")


combat_counts <- ComBat_seq(counts = as.matrix(my_RawCounts),
                            batch = my_metadata$batch)
# Build DESeq on ComBat-corrected counts (if you choose this route)
dds_combat <- DESeqDataSetFromMatrix(countData = combat_counts,
                                     colData = my_metadata,
                                     design = ~ condition)   # now omit batch because already corrected
dds_combat <- DESeq(dds_combat)
res_sputum_vs_Ra <- results(dds_combat, contrast = c("condition", "Sputum", "Broth"))
res_sputum_vs_Rv7 <- results(dds_combat, contrast = c("condition", "Sputum", "Rv7"))


