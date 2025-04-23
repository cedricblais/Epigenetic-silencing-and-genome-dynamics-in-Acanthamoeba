library(ggplot2)
library(ggprism)
library(patchwork)
library(magrittr)
library(dplyr)
library(tidyverse)

expression_table_1 <- read.table("C3_TPM_with_info.tsv", header=TRUE)

expression_table_1 <- expression_table_1 %>% add_column(Organism = "C3", .after = 'type.conservation')


expression_table_2 <- read.table("Neff_TPM_with_info.tsv", header=TRUE)

expression_table_2 <- expression_table_2 %>% add_column(Organism = "Neff", .after = 'type.conservation')



expression_table <- rbind(expression_table_1[1:7], expression_table_2[1:7])


#We integrated degraded with conserved:

expression_table$type.conservation[expression_table$type.conservation == 'gene,DEGRADED'] <- 'gene,CONSERVED'

expression_table$type.conservation[expression_table$type.conservation == 'virus_gene,DEGRADED'] <- 'virus_gene,CONSERVED'

table(expression_table$type.conservation)



expression_table <- expression_table %>% mutate(viral = ifelse(str_detect(expression_table$type.conservation, "virus"), "Viral", "Non-viral"))


label_list <- c("Conserved\nnon-viral genes", "Unconserved\nnon-viral genes", "Conserved\nviral genes", "Unconserved\nviral genes")

my_colors= c('#33adff' ,'#0066ff','#47d147', 'red')

p <- ggplot(data=expression_table, mapping=aes(x=type.conservation, y=log2(tpm + 1))) + geom_violin(aes(fill = type.conservation), show.legend = FALSE,scale = 'width') + 
  scale_fill_manual(values = c('#88bdfa', '#2477da', '#f27c8d', '#B1172C')) +xlab(NULL)+ylab("TPM") + geom_boxplot(width = 0.075, outliers = FALSE)+ 
  scale_x_discrete(label = label_list) + facet_wrap(~ Organism, scales = "free") + theme_minimal() + theme(strip.text = element_text(size=15), axis.text.x = element_text(size = 12));

df_p_val <- expression_table %>%
  rstatix::group_by(Organism) %>%
  rstatix::t_test(tpm ~ type.conservation)

p2 <- p + add_pvalue(df_p_val, label.size = 5,
               y.position = c(18, 16.8, 16, 18, 15, 18, 18, 16.8, 16, 18, 15, 18), bracket.shorten = c(0.025, 0, 0, 0.025, 0, 0.025, 0.025, 0, 0, 0.025, 0, 0.025))

p2

