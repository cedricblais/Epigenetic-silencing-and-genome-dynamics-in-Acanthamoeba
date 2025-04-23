library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggcorrplot)
library(ggplotify)
library(grid)


#We load the file, ignoring the header
#Neff: Neff_repeatmasker_complex_repeats_formatted.tsv
#C3: C3_repeatmasker_complex_repeats_formatted.tsv
transposable_elements <- read.table("Neff_repeatmasker_complex_repeats_formatted.tsv",sep="\t", 
                                    fill = TRUE, skip = 3)

#We remove everything after "/" in the family column
transposable_elements$V11 <- gsub("\\/.*", "", transposable_elements$V11)

#Divergent regions
#Neff: neff_c3_mummer_concatenated_divergent_regions.gff
#C3: c3_neff_mummer_concatenated_divergent_regions.gff
divergent_regions <- read.table("neff_c3_mummer_concatenated_divergent_regions.gff")

#Curated viral regions
#Neff: neff_viral_regions_curated.gff
#C3: c3_viral_regions_curated.gff
viral_table <- read.table("neff_viral_regions_curated.gff")




#We add a 17th column, which specified whether the element is in a viral gap, non-viral gap, or is conserved
transposable_elements <- transposable_elements %>%
  rowwise() %>%
  mutate(V17 = ifelse (any(((V6 >= viral_table$V4 & V6 <= viral_table$V5)|(V7 >= viral_table$V4 & V7 <= viral_table$V5)|(V6 < viral_table$V4 & V7 > viral_table$V5)) & V5 == viral_table$V1), "viral_gap", ifelse (any(((V6 >= divergent_regions$V4 & V6 <= divergent_regions$V5)|(V7 >= divergent_regions$V4 & V7 <= divergent_regions$V5)|(V6 < divergent_regions$V4 & V7 > divergent_regions$V5)) & V5 == divergent_regions$V1), 'non_viral_gap', 'conserved')))

#We get the total number of predicted transposable elements
length <- nrow(transposable_elements)

#We get counts for each pair of mobile element and type of region
counts <- transposable_elements %>%
  group_by(V17, V11) %>%
  summarise(N = n(), .groups = "drop")

freq_table2 <- table(transposable_elements$V17, transposable_elements$V11)

# Chi-square test
chi_square_test <- chisq.test(freq_table2, simulate.p.value = TRUE)

residuals_matrix <- chi_square_test$stdres

# Convert to tidy format
residuals_df <- as.data.frame(as.table(residuals_matrix)) %>%
  rename(Row = Var1, Column = Var2, Residual = Freq)

#label list 
label_list = c("conserved region", "divergent region", "viral region")

# Bubble plot
bubblesneff <- ggplot(residuals_df, aes(x = Column, y = Row)) +
  geom_point(aes(size = abs(Residual), fill = Residual), shape = 21, color = "black") +
  scale_size_continuous(range = c(2, 15), guide = "none") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, name = "Standardized\nResidual") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "Enrichment of mobile elements based on Chi-Square test",
       x = "Repeatmasker category", y = "", subtitle = paste("simulated pval =", round(chi_square_test$p.value, digits = 4)))+theme(plot.subtitle = element_text(hjust = 0.5, vjust = 3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_y_discrete(label = label_list) 


# Bar plot
totals <- counts %>%
  group_by(V17) %>%
  summarise(total = sum(N))

library(RColorBrewer)

gneff <- ggplot(counts, aes(x = V17, y = N, fill = V11)) +
  geom_bar(position="fill", stat = "identity", width = 0.5) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, vjust = 3), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 15), axis.title.y = element_text(vjust=2, size = 10))+
  xlab(NULL) +
  labs(title = "Counts of Col1 Values by Col2 Group", fill = 'Mobile element') +
  geom_text(data = totals,
            aes(x = V17, y = 1.05, label = paste("n =" ,total)),  # y = 1.05 places text just above the top (1.0)
            inherit.aes = FALSE) + ylab("Proportion of mobile element") + 
  ggtitle('Mobile elements in conserved, divergent, \nand viral regions in Neff') + 
  scale_fill_brewer(palette="Accent") + scale_x_discrete(label = label_list)
gneff
bubblesneff
ggarrange(gneff, bubblesneff, nrow = 2)


