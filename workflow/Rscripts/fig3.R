suppressPackageStartupMessages({
  library("argparser")
  library("tidyverse")
})

p <- arg_parser(hide.opts = TRUE, 
                paste0("Generates figure 3"))
p <- add_argument(p, "--midas_taxonomy",
                  short = "-g",
                  help = "Tfile with MiDAS taxonomy")
p <- add_argument(p, "--bgcs_table",
                  short = "-b",
                  help = "Table with bgcs produced by bgcs_table.R ")
p <- add_argument(p, "--mags_table",
                  short = "-m",
                  help = "Table with mags produced by mags_table.R")
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

gtdbtax <- read_csv(argv$mags_table,
                    col_select = c("MAG", "GTDB_tax"), 
                    col_types = "cc")


filamentous <- c(
  "Ca_Microthrix",
  "Ca_Amarolinea",
  "Ca_Promineofilum", 
  "Ca_Villigracilis")

nitrifiers <- c(
  "Nitrosomonas",
  "Nitrospira",
  "Nitrotoga")

denitrifiers <- c(
  "Rhodoferax",
  "Acidovorax",
  "Zoogloea")

PAOs <- c(
  "Tetrasphaera",
  "Dechloromonas",
  "Ca_Accumulibacter")

GAOs <- c(
  "Ca_Competibacter",
  "Micropruina",
  "Defluviicoccus")

pets <- list(
  filamentous, 
  nitrifiers, 
  denitrifiers, 
  PAOs, 
  GAOs)

guilds <- c(
  "Filamentous", 
  "Nitrifiers", 
  "Denitrifiers", 
  "PAOs", 
  "GAOs")

names(pets) <- guilds

guilds_df <- bind_rows(lapply(pets, tibble), .id = "guild")
colnames(guilds_df) <- c("guild", "genus")

midas <- read_tsv(argv$midas_taxonomy, 
                  col_select = 1:2, 
                  col_types = "cc", 
                  col_names = F)

colnames(midas) <- c("MAG", "midas_tax")

midas <-  midas %>%
  mutate(genus = str_extract(midas_tax, 'g:.*?,'), 
  genus = str_sub(genus, end = -2, start = 3),  
  MAG = str_sub(MAG, end = -4), 
  .keep = "unused")

ds <- read_csv(argv$bgcs_table, 
               col_select = c(3,4), 
               col_types = "cc") %>%
  mutate(class = factor(class, levels = bigscape_classes), 
         class = my_classes[as.numeric(class)], 
         class = factor(class, levels = unique(my_classes))) %>%
  group_by(MAG, class) %>%
  tally()


genera_list <- unlist(pets)

pets_ds <- 
  expand.grid(MAG = unique(ds$MAG), 
              class = unique(my_classes))%>%
  left_join(midas, by = "MAG") %>%
  left_join(ds, by = c("MAG", "class")) %>%
  mutate(n = replace_na(n, 0))%>%
  left_join(guilds_df, by = "genus") %>% 
  filter(is.na(guild) == FALSE)  %>%
  mutate(guild = factor(guild, levels = guilds)) %>%
  group_by(genus) %>%
  mutate(mags = length(unique(MAG)), 
         string = paste0(genus, ' (', mags,')'))


figure <- ggplot(pets_ds) +
  geom_boxplot(aes(x = n, y = string, color = class), 
               fatten = 2.5) +
  facet_grid(guild ~ class, scales = "free", space = "free") +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(breaks = seq(0, 8, by =2)) +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.box = "horizontal",
        panel.background =  element_rect(color = "gray"), 
        strip.text.x = element_blank()) +
  guides(color = guide_legend(nrow = 1)) +
  xlab("BGCs") +
  ylab("Genus") 
 
ggsave(figure,
       filename = argv$output,
       height= 6.9,
       width = 9, 
       units = "in")

