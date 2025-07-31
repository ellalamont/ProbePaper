# Make a graph zoomed in on the H37Ra SNP
# 7/31/25


# App from ProbeTest5 is showing:
# Mutation from C -> T at 852263
# Since it starts at 851608, this is at position 655
# "one single nucleotide polymorphism affecting codon 219 (TCG → TTG) in phoP and changing a serine to a leucine. S219 in PhoP of M. tuberculosis is equivalent to R200 in PhoB of E. coli" - https://pmc.ncbi.nlm.nih.gov/articles/PMC2238218/


source("Import_Data_SNPs.R") # SNPs_df_4

SNPs_df_4$SampleID %>% unique()
# [1] "THP1_1e6_1a" "THP1_1e6_2b" "THP1_1e6_3a"


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
        plot.subtitle = element_text(size=9), 
        # plot.margin = margin(10, 10, 10, 20),
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )


###########################################################
############ SUBSET TO REGION OF INTEREST #################

SNPs_df_5 <- SNPs_df_4 %>% filter(between(POSITION, 852256, 852270))



###########################################################
####################### MAKE GRAPHs #######################

# my_sample <- "THP1_1e6_1a"
# my_sample <- "THP1_1e6_2b"
my_sample <- "THP1_1e6_3a"
my_barplot <- SNPs_df_5 %>%
  filter(SampleID == my_sample) %>%
  mutate(POSITION = as.character(POSITION)) %>% 
  ggplot(aes(x=POSITION, y=Proportion, fill=base_call, text = paste0("Reference Base: ", REF_BASE))) + 
  geom_bar(stat="identity", color = "black") + # position = stack for raw reads, position = fill for proportion
  geom_text(aes(label = ifelse(Proportion > 0.5, CALL_BASE, ""), y = 1.05),
            size = 4) +
  scale_fill_manual(values=my_colors) + 
  scale_y_continuous(limits = c(0,1.07), breaks = seq(0, 1.07, 0.2)) +
  labs(title = paste0(my_sample, " phoP (Rv0757) region with SNP"),
       subtitle = "Ra has mutation from C -> T at 852263, codon 219 (TCG -> TTG)",
       x = "Base position",
       y = "Proportion of reads") +
  my_plot_themes
# scale_x_continuous(breaks = seq(from = input$range[1], to = input$range[2], by = 2))
my_barplot
ggsave(my_barplot,
       file = paste0(my_sample, "_phoP_SNP_calls.pdf"),
       path = "Figures/PhoP_SNP",
       width = 7, height = 4, units = "in")




