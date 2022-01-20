suppressPackageStartupMessages({
  library("tidyverse")
  library("argparser") 
  library("parallel") 
  #library("GenomicRanges"), 
  library("jsonlite") 
  library("pbapply")
  })

# Parse command line arguments
p <- arg_parser(hide.opts = TRUE, 
  paste0("Take an antiSMASH and BiG-SCAPE output directories, ",
         "and search recursively in them ", 
         "to produce a csv table where each observation is a BGCc."))

p <- add_argument(p, "--antismash_dir",
                  short = "-a",
                  help = "Directory with antiSMASH-generated json files")
p <- add_argument(p, "--bigscape_dir",
                  short = "-b",
                  help = "Directory with BiG-SCAPE output")
p <- add_argument(p, "--output", 
                  short = "-o", 
                  help = "Filepath of output table ")
p <- add_argument(p, "--threads",
                  short = "-t",
                  help = "Number of threads to use in parallel", 
                  default = detectCores())

argv <- parse_args(p)


#Extracts info from antismash .json files

get_region_rows <-   function(MAG){
    json_file <- str_replace(MAG, '$','.json')
    filepath <- paste0(argv$antismash_dir, '/', json_file)
    MAG_json <- fromJSON(filepath)
    records_df <- MAG_json$records
    contig_names <- records_df$id
    features <- records_df$features
    names(features) <- contig_names
    features_df <- bind_rows(features, .id = 'contig')
    regions_df <- features_df[features_df$type == 'region',]
    ends <- str_sub(regions_df$location, start = 2, end = -2) %>%
      str_split(":", n =2, simplify = T)
    start <- as.numeric(ends[,1])
    end <- as.numeric(ends[,2])
    length <- as.numeric(ends[,2]) - as.numeric(ends[,1])
    product <- sapply(regions_df$qualifiers$product, paste0, collapse = ";") %>%
      as.character()
    contig_edge <- as.logical(unlist(regions_df$qualifiers$contig_edge))
    contig <- as.character(regions_df$contig)
    regions <- suppressMessages(bind_cols(
    contig, start, end, product, contig_edge))
    colnames(regions) <- c(
     "contig", "start", "end", "product", "contig_edge")
    regions
   
#    modules_df <- features_df[features_df$type == 'aSModule',]
#    modules_df <- modules_df[modules_df$qualifiers$type == 'pks',]
#    module_ends <- str_sub(modules_df$location, start = 2, end = -2) %>%
#      str_split(":", n =2, simplify = T)
#    module_start <- as.numeric(module_ends[,1])
#    module_end <- as.numeric(module_ends[,2])
#    module_contig <- as.character(modules_df$contig)
#    
#    regions_ranges <-  GRanges(
#      seqnames = regions$contig, 
#      ranges = IRanges(
#        start = regions$start,
#        end = regions$end, 
#        names = rownames(regions)))
#    
#    modules_ranges <- GRanges(
#      seqnames = module_contig,
#      ranges = IRanges(
#        start = module_start, 
#        end = module_end))
#    
#    module_count <-  table(subjectHits(findOverlaps(
#      modules_ranges, regions_ranges)))
#    
#    regions$modules <- module_count[
#      match(rownames(regions), names(module_count))]
#    
#  mutate(regions, modules = as.integer(modules))
  }

MAGs <- list.files(argv$antismash_dir, recursive = T, pattern = '*.json') %>%
  str_sub(end = -6)

# Parallelizes json parsing, shows a progress bar on the terminal

cluster <- makeForkCluster(nnodes = argv$threads)
pbo  <- pboptions(type="timer")
regions_list <- pbsapply(MAGs, get_region_rows, simplify = FALSE, cl  = cluster)

regions <- bind_rows(regions_list, .id = 'MAG') %>%
  group_by(contig) %>%
  mutate(bgc_id = paste0(
    contig, 
    '.region', 
    str_pad(1:n(), width = 3, pad = '0')))



annotation_files <- list.files(path = argv$bigscape_dir,
                               pattern = "Network_Annotations",
                               recursive = TRUE, 
                               full.names = TRUE)

bigscape_annotations <- distinct(
  bind_rows(
  lapply(annotation_files, 
         read_tsv,
         col_select = c(1,5),
         col_names =  TRUE,
         col_types = "cc"))) 

colnames(bigscape_annotations) <- c("bgc_id", "class")


clustering_files <- list.files(path = argv$bigscape_dir,
                               pattern = "clustering",
                               recursive = TRUE, 
                               full.names = TRUE)

bigscape_clustering <- distinct(
  bind_rows(
  lapply(clustering_files, 
         read_tsv, 
         col_types = 'cc')))

colnames(bigscape_clustering) <- c("bgc_id", "GCF")

bgcs_dataframe <- bigscape_clustering %>% 
  group_by(bgc_id) %>% 
  summarise(GCF = paste0(unique(GCF), collapse =';')) %>%
  left_join(bigscape_annotations, by = "bgc_id") %>%
  left_join(regions, by = 'bgc_id')


write_csv(bgcs_dataframe, 
          file = argv$output, 
          num_threads = argv$threads)




