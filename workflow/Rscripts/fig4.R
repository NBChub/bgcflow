suppressPackageStartupMessages({
  library("argparser")
  library("tidyverse")
  library("treeio")
  library("ggtree")
  library("ape")
  library("ggtreeExtra")
})

p <- arg_parser(hide.opts = TRUE, 
                paste0("Generates figure 3. ", 
                "Subsets the tree from Fig2 to the selected taxon, ", 
                "labels species according to MiDAS 4.8 taxonomy ", 
                "and displays a matrix with absence/presence of GCFs"))
p <- add_argument(p, "--taxon",
                  help = "Taxon to analyse, ej: g__Nitrospira for genus Nitrospira", 
                  short = "-t")
p <- add_argument(p, "--midas_taxonomy",
                  short = "-m",
                  help = "File with MiDAS taxonomy")
p <- add_argument(p, "--bgcs_table",
                  short = "-b",
                  help = "Table with bgcs produced by bgcs_table.R")
p <- add_argument(p, "--gtdb_tree",
                  short = "-g",
                  help = "GTDB-tk tree file")
p <- add_argument(p, "--output", 
                  short = "-o", 
                  help = "Filepath for output png figure")

argv <- parse_args(p)

bigscape_classes <- c("Terpene",
                      "RiPPs", 
                      "NRPS",
                      "PKS-NRP_Hybrids",
                      "PKSI",
                      "PKSother", 
                      "Saccharides",
                      "Others")                             

my_classes <- c("Terpene", 
                "RiPP", 
                "NRPS", 
                "NRPS-PKS", 
                "PKS", 
                "PKS", 
                "Other", 
                "Other")


my_colors <- c("#ff8800",
               "#20c200",
               "#ff0000",
               "#dd00ff",
               "#0d00ff",
               "#757575")

midas <- read_tsv(argv$midas_taxonomy, 
                  col_select = 1:2, 
                  col_types = "cc", 
                  col_names = FALSE)%>%
  mutate(MAG = str_sub(X1, end = -4), 
         tax = str_replace_all(X2, ":", '__'), 
         tax = str_replace_all(tax, ',', ';'), .keep = "none")

taxon_mags <- midas %>%
  filter(str_detect(tax, argv$taxon)) %>%
  pull(MAG)


tree<- gsub(';', '_', readLines(argv$gtdb_tree))
tree <- gsub('$', ';', tree)
tree <- read.newick(text = tree)
tree$edge.length <- replace_na(tree$edge.length, 0)
reference_tips <- tree$tip.label[str_detect(tree$tip.label, '18-Q3') == FALSE]
tree <- drop.tip(tree, reference_tips)
drop_tips <- tree$tip.label[tree$tip.label %in% taxon_mags == FALSE]
tree <- drop.tip(tree, drop_tips)
tree <- ape::drop.tip(tree, drop_tips)
tree1 <- ggtree(tree)

midas_species_ds <- midas %>%
  filter(MAG %in% taxon_mags) %>%
  mutate(MAG = MAG, 
         midas_species  = str_extract(tax, "s__.*$"),
         midas_species = str_sub(midas_species, start = 4, end = -2),
         .keep = "none") %>%
  group_by(midas_species) %>%
  mutate(node = MRCA(tree, MAG))


tree2 <- tree1 +
  geom_highlight(data = distinct(midas_species_ds, node, midas_species), 
                 mapping = aes(
                   node = node),
                 fill = "grey", 
                 alpha = 0.3) +
  geom_cladelab(data = distinct(midas_species_ds, node, midas_species), 
                mapping = aes(
                  node = node, 
                  label = midas_species), 
                geom = 'shadowtext',
                hjust = 1,
                bg.colour = 'white',
                align = TRUE,
                offset = -0.01,
                barsize = NA) 


matrix_ds <- read_csv(argv$bgcs_table, 
                      col_select = c("MAG", "class", "GCF", "product"), 
                      col_types = "cccc")%>%
  filter(MAG %in% taxon_mags) %>%
  mutate(class = factor(class, levels = bigscape_classes), 
         class = my_classes[as.numeric(class)], 
         class = factor(class, levels = unique(my_classes))) %>%
  arrange(class, product) %>%
  mutate(GCF = factor(GCF, levels = unique(GCF))) %>%
  arrange(GCF) %>%
  mutate(string = paste0(GCF, ": ", product), 
         string = factor(string, levels = unique(string))) %>%
  select(MAG, string, class)

matrix_ds <- expand_grid(matrix_ds$string, matrix_ds$MAG) %>%
  setNames(c("string", "MAG")) %>%
  left_join(matrix_ds, by = c("string", "MAG")) 


figure <- tree2 + geom_fruit(
  geom = "geom_point", 
  data = matrix_ds, 
  mapping = aes(
    fill = class,
    x = string, 
    y = MAG), 
  size = 3.5, 
  shape = 21, 
  pwidth = 2.5,
  axis.params = list(axis = "x", 
                     text.size = 3, 
                     text.angle = 45, 
                     hjust =1),
  grid.params = list(vline = TRUE, 
                     size = 1, 
                     color = "gray95"),
  color = "gray95") +
  scale_fill_manual(values = my_colors, 
                    na.value = "gray95")  +
  guides(fill=guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal") +
  scale_y_continuous(expand = c(0, 4, 0, 0))


ggsave(figure,
       filename = argv$output,
       device = "png", 
       height= 6.9,
       width = 9, 
       units = "in")

