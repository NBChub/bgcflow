library(tidyverse, quietly=TRUE)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Joins GTDBtk summary tables for bacteria and archaea")

# Add command line arguments

p <- add_argument(p, "--output_directory", short = "-o", help = "output_directory")
p <- add_argument(p, "--archaea", short = "-a", help = "Summary file for archaea")
p <- add_argument(p, "--bacteria", short = "-b", help = "Summary file for bacteria")


# Parse the command line arguments
argv <- parse_args(p)

bac202 <- read_tsv(argv$bacteria)
ar202 <- read_tsv(argv$archaea)

MAGs <- bind_rows(bac202[,c("user_genome", "classification")],
                  ar202[,c("user_genome", "classification")])

write.csv(MAGs, file.path(argv$output_directory, "mag_stats.csv"))