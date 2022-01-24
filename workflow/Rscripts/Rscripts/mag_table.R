suppressPackageStartupMessages({
library("argparser")
library("tidyverse")
library("readxl")
})

# Parse command line arguments
p <- arg_parser(hide.opts = TRUE, 
                paste0("Take a table produced by bgcs_table.R ",
                       "and a GTDB-Tk classify output directory ", 
                       "to produce a csv table where each observation is a MAG"))

p <- add_argument(p, "--bgcs_table",
                  short = "-b",
                  help = "Table with bgcs produced by bgcs_table.R ")
p <- add_argument(p, "--gtdbtk_dir",
                  short = "-g",
                  help = "Directory with GTDB-Tk summary.tsv files")
p <- add_argument(p, "--supplementary_file", 
                  short = "-s", 
                  help = "Filepath with supplementary info from Singleton et al 2021")
p <- add_argument(p, "--output", 
                  short = "-o", 
                  help = "Filepath for output table ")


argv <- parse_args(p)

gtdbtk_files <- list.files(argv$gtdbtk_dir, 
           pattern = ".summary.tsv", 
           recursive = TRUE, 
           full.names = TRUE)

gtdb_tax <- distinct(
  bind_rows(
    lapply(
      gtdbtk_files, 
      read_tsv,  
      col_select = 1:2, 
      col_types = "cc")))

colnames(gtdb_tax) <- c("MAG", "GTDB_tax")

size <- read_xlsx(
  argv$supplementary_file, 
  col_names = TRUE,
  skip = 1)[,c(1,4)] %>%
  mutate(MAG = str_replace(MAG, ".fa", ""))

bgc_count <- read_csv(argv$bgcs_table, 
                      col_types = "cccccddcl") %>%
  group_by(MAG) %>%
  summarise(bgcs = n())

mags <- left_join(gtdb_tax, bgc_count, by = "MAG") %>%
  mutate(bgcs = replace_na(bgcs, 0)) %>%
  left_join(size, by = "MAG")
  
write_csv(mags, file = argv$output)

  

