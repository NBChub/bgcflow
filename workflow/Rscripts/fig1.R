suppressPackageStartupMessages({
  library("argparser")
  library("tidyverse")
  library("randomcoloR")
  library("cowplot")
  library("shadowtext")
})

p <- arg_parser(hide.opts = TRUE, 
                paste0("Generates figure 1"))
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

ds <- read_csv(argv$bgcs_table, 
               col_select = c(1,3,4,8), 
               col_types = "cccc") %>%
  mutate(class = factor(class, levels = bigscape_classes), 
         class = my_classes[as.numeric(class)], 
         class = factor(class, levels = unique(my_classes))) 

get_color <- function(class) {
  #Randomizes colors for plot
  hue <-  if(class == 'Terpene'){"orange"}
          else if(class == 'RiPP'){"green"} 
          else if(class == 'NRPS'){"red"}
          else if(class == 'NRPS-PKS'){"purple"}
          else if(class == 'PKS'){"blue"} 
          else if(class == 'Other'){"monochrome"}
          else{NULL}
    color <- randomcoloR::randomColor(hue = hue)
  color <- ifelse(class == "Terpene", "#ff8800", color)
  color
  }

get_colors <- Vectorize(get_color)

barplot_ds <- ds %>%
  group_by(class, product) %>% 
  tally() %>%
  mutate(product = str_replace_all(product, ";", " / "),
         product = str_replace_all(product, "_", " ")) %>%
  arrange(class, desc(n)) %>%
  mutate(product = factor(product, levels = rev(product))) %>% 
  mutate(ymax = cumsum(n), 
         ymin = ymax - n,
         xmin = as.integer(class),
         xmax = xmin + 0.95) %>%
  ungroup() %>% 
  mutate(fill = get_colors(class)) %>%
  group_by(class) %>%
  mutate(label = paste0(class, "\n(", sum(n), ")"))

 barplot <- ggplot(barplot_ds, 
       aes(fill = product, y = n, x = 0)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_shadowtext(data = ~ filter(.x, n > 25), 
                  aes(label = product, 
                      y = ymin +(ymax-ymin)/2, 
                      x = 0), 
                  size = 3,
                  color = "black",
                  bg.colour = "white", 
                  bg.alpha = 0.5,
                  bg.r = 0.1) +
  theme(legend.position = "none") +
  facet_grid( . ~ class) +
  scale_fill_manual(values = barplot_ds$fill) +
  scale_y_continuous(breaks = seq(0, 1200, by = 200)) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        strip.placement = "inside") +
  ylab("BGC count") +
  xlab("BGC class")
  
boxplot_ds <- read_csv(argv$mags_table, col_types = "ccd") %>%
  mutate(TotBP = TotBP / 1000000, 
         phylum = str_extract(GTDB_tax, 'p__.*?;'), ) %>%
  group_by(phylum) %>%
  mutate(n = n(), 
         phylum = if_else(n > 7, phylum, "   other phyla ")) %>%
  group_by(phylum) %>%
  mutate(n = n(), 
         string = paste0(
    str_sub(phylum, start = 4, end = -2),
    " (", n, " )")) %>%
  arrange(desc(n))

boxplots <- ggplot(boxplot_ds, 
                   aes(x = fct_infreq(string), 
                       y = bgcs, 
                       color = TotBP)) + 
  geom_boxplot(fill = NA,
               outlier.shape = NA, 
               fatten = 7,
               width =0.9, 
               size = 0.2) +
  geom_jitter(height = 0.1, width = 0.25) +
  scale_color_gradient2(low = "steelblue",
                        mid = "lemonchiffon",
                        high = "red3", 
                        midpoint = 6, 
                        limits = c(0,12)) +
  scale_x_discrete(expand = c(0.1,0, 0, 0)) +
  scale_y_continuous(breaks = seq(0, 22, by = 4)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   size =12, 
                                   vjust = 0.5), 
        legend.position = "top", 
        legend.box = "horizontal") +
  ylab('BGCs per MAG') +
  xlab('Phylum') +
  guides(color = guide_colorbar(title = "Genome size (Mbp)"))


figure <- cowplot::plot_grid(barplot,
                   boxplots, 
                   ncol =2,
                   align = "h",
                   axis = "l", 
                   labels = c("A", "B")) 

ggsave(argv$output, 
       figure, 
       device = "png",
       height = 6.9,
       width = 9, 
       units = "in")
