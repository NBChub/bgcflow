#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gggenomes))
suppressPackageStartupMessages(library(jsonlite))

region_feat <- read_feats(commandArgs(trailingOnly = TRUE)[1])
cogs <- read_csv(commandArgs(trailingOnly = TRUE)[2])

cogs <- cogs %>%
  mutate(
    cluster_label = paste0(cluster_id, " (", cluster_n, ")")
  )

merged_feat <- merge(region_feat, cogs, by.x = "feat_id", by.y = "feature_id")

# Output results or further processing...
output_file <- commandArgs(trailingOnly = TRUE)[3]
write_json(merged_feat, output_file)