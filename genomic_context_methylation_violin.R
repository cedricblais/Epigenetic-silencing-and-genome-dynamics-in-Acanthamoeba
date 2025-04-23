library(ggplot2)
library(rstatix)
library(tidyverse)
library(ggprism)


#### We load the necessary files #####

#Curated viral regions
#Neff: neff_viral_regions_curated.gff
#C3: c3_viral_regions_curated.gff
viral_table <- read.table("neff_viral_regions_curated.gff")

#All concatenated divergent regions
#Neff: neff_c3_mummer_concatenated_divergent_regions.gff
#C3: c3_neff_mummer_concatenated_divergent_regions.gff
divergent_table <- read.table("c3_neff_mummer_concatenated_divergent_regions.gff")

#Methylation
#Neff: Neff_methylation_frequency_split.tsv
#C3: C3_methylation_frequency_split.tsv
methyl_table <- read.table("C3_methylation_frequency_split.tsv", header=FALSE)

#We get mobile elements, ignoring the header lines
#Neff: Neff_repeatmasker_complex_repeats_formatted.tsv
#C3: C3_repeatmasker_complex_repeats_formatted.tsv
transposable_elements <- read.table("C3_repeatmasker_complex_repeats_formatted.tsv",sep="\t", 
                                    fill = TRUE, skip = 3)

#Gene bed file
#Neff: Neff_all_genes_info_virus_conservation.bed
#C3: C3_all_genes_info_virus_conservation.bed
genes <- read.table(text = gsub(";","\t", readLines("C3_all_genes_info_virus_conservation.bed")))

#We make it 1 based
genes$V2 <- genes$V2 +1

#We remove everything after "/"
transposable_elements$V11 <- gsub("\\/.*", "", transposable_elements$V11)


#We add a 17th column to the methylation, which specified whether the CpG site is in a viral gap, non-viral gap, or is conserved
methyl_table <- methyl_table %>%
  rowwise() %>%
  mutate(V9 = ifelse (any(V2 >= viral_table$V4 & V2 <= viral_table$V5 & V1 == viral_table$V1), "viral_gap", ifelse(any(V2 >= divergent_table$V4 & V2 <= divergent_table$V5 & V1 == divergent_table$V1), "non_viral_gap", 'conserved')))

methyl_viral <- methyl_table

methyl_viral <- methyl_viral %>%
  rowwise() %>%
  mutate(V10 = ifelse (any(V2 >= genes$V2 & V2 <= genes$V3 & V1 == genes$V1 & genes$V4 == "virus_gene"), "virus_gene", ifelse(any(V2 >= genes$V2 & V2 <= genes$V3 & V1 == genes$V1 & genes$V4 == "gene"), "gene", ifelse(any(V2 >= transposable_elements$V6 & V2 <= transposable_elements$V7 & V1 == transposable_elements$V5), "mobile element", "intergenic"))))

methyl_viral2 <- methyl_viral
methyl_viral$V11 <- paste(methyl_viral$V9,methyl_viral$V10)

methylation_df <- as.data.frame(methyl_viral)

label_list <- c("genes", "intergenic\nregions", "complex\nrepeats", "viral\ngenes")

p <- ggplot(data=methylation_df, mapping=aes(x=V10, y=V7)) + geom_violin(aes(fill = V10,), show.legend = FALSE,scale = 'width') +
  stat_summary(fun.y=mean, geom="segment", mapping=aes(xend=..x.. - 0.15, yend=..y..), color='red') +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.15, yend=..y..), color='red') +
  stat_summary(fun.y=median, geom="segment", mapping=aes(xend=..x.. - 0.15, yend=..y..)) +
  stat_summary(fun="median", geom="segment", mapping=aes(xend=..x.. + 0.15, yend=..y..)) +
  xlab(NULL)+ylab("Methylation") + #geom_boxplot(width = 0.05, outliers = FALSE)+ 
  theme_minimal() + theme(strip.text = element_text(size=15), axis.text.x = element_text(size = 12)) + facet_wrap(~ V9, scales = "free", labeller = labeller(V9 = 
                                                                                                                                                               c("viral_gap" = "Viral regions",
                                                                                                                                                                 "non_viral_gap" = "Divergent regions",
                                                                                                                                                                 "conserved" = "Conserved regions"
                                                                                                                                                               ))) +
  scale_y_continuous(breaks = pretty(c(0, 1), n = 5)) +
  scale_x_discrete(label = label_list);
#p
p
df_p_val1 <- methylation_df %>%
  rstatix::group_by(V9) %>%
  rstatix::t_test(V7 ~ V10)

#df_p_val1

p2 <- p + add_pvalue(df_p_val1, label.size = 3,
                     y.position = c(1.05, 1.1, 1.15, 1.05, 1.2, 1.05, 1.05, 1.1, 1.15, 1.05, 1.2, 1.05, 1.05, 1.1, 1.15, 1.05, 1.2, 1.05), bracket.shorten = c(0.025, 0, 0, 0.025, 0, 0.025, 0.025, 0, 0, 0.025, 0, 0.025, 0.025, 0, 0, 0.025, 0, 0.025))
p2
