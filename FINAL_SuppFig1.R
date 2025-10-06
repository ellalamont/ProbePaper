# Supplemental figure 1


source("FINAL_ImportData.R")

# Mutation from C -> T at 852263
# Since it starts at 851608, this is at position 655
# "one single nucleotide polymorphism affecting codon 219 (TCG → TTG) in phoP and changing a serine to a leucine. S219 in PhoP of M. tuberculosis is equivalent to R200 in PhoB of E. coli" - https://pmc.ncbi.nlm.nih.gov/articles/PMC2238218/

SNPs_df <- DEG_dfs$phoP.SNPs
levels(as.factor(SNPs_df$SampleID)) # This only contains the 3 1e6 H37Ra spiked THP1 samples that were captured
SNPs_df$SampleID <- gsub("_S(.*)", "", as.character(SNPs_df$SampleID))


###########################################################
##################### PROCESS DATA ########################

names(SNPs_df)[names(SNPs_df) == "X."] <- "REF"
SNPs_df_2 <- SNPs_df %>% pivot_longer(cols = c("A","C", "G", "T", "N", "Indel"),
                                      names_to = "base_call",
                                      values_to = "raw_read_number")
SNPs_df_3 <- SNPs_df_2 %>% mutate(raw_read_number = case_when(REF_BASE == base_call ~ REF, TRUE ~ raw_read_number))
SNPs_df_4 <- SNPs_df_3 %>% mutate(base_call = case_when(REF_BASE == base_call ~ "Genome", TRUE ~ base_call))
ordered_bases <- c("Genome", "A", "C", "G", "T", "N", "Indel")
SNPs_df_4$base_call <- factor(SNPs_df_4$base_call, levels = ordered_bases)
SNPs_df_4 <- SNPs_df_4 %>% mutate(Proportion = raw_read_number / DEPTH)
SNPs_df_5 <- SNPs_df_4 %>% filter(between(POSITION, 852256, 852270))



###########################################################
####################### MAKE GRAPHs #######################

my_colors <- c("grey", "green", "orange", "blue", "red", "black", "purple")

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=6, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9))


# my_sample <- "THP1_1e6_1a"
# my_sample <- "THP1_1e6_2b"
my_sample <- "THP1_1e6_3a"
SuppFig1 <- SNPs_df_5 %>%
  filter(SampleID == my_sample) %>%
  mutate(POSITION = as.character(POSITION)) %>% 
  ggplot(aes(x=POSITION, y=Proportion, fill=base_call, text = paste0("Reference Base: ", REF_BASE))) + 
  geom_bar(stat="identity", color = "black") + 
  geom_text(aes(label = ifelse(Proportion > 0.5, CALL_BASE, ""), y = 1.05),
            size = 4) +
  scale_fill_manual(values=my_colors) + 
  scale_y_continuous(limits = c(0,1.07), breaks = seq(0, 1.07, 0.2)) +
  labs(title = paste0(my_sample, " phoP (Rv0757) region with SNP"),
       subtitle = "Ra has mutation from C -> T at 852263, codon 219 (TCG -> TTG)",
       x = "Base position",
       y = "Proportion of reads") +
  my_plot_themes
SuppFig1

