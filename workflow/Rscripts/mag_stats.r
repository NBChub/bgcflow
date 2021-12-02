library(dplyr)
library(readr)

old_stats <- read_tsv("MAG_statistics_STABLEX_20200911_RS95.tsv")
bac202 <- read_tsv("gtdbtk.bac120.summary.tsv")
ar202 <- read_tsv("gtdbtk.ar122.summary.tsv")
midas4.8 <- read_tsv("MiDAS4.8_reformatted_correct.txt", col_names = c("MAG", "midas4.8tax", "midas_identity")) %>%
  mutate(MAG = str_sub(MAG, end = -4)) %>%
  select(1:3)

MAGs <- old_stats %>%
  mutate(MAG = str_sub(MAG, end = -4)) %>%
  left_join(
    bind_rows(bac202[,c("user_genome", "classification")],
             ar202[,c("user_genome", "classification")]), 
    by = c('MAG' = 'user_genome')) %>%
  select(MAG, TotBP, gtdb202tax = classification) %>%
  left_join(midas4.8, by = "MAG")

write.csv(MAGs, "mag_stats.csv")
