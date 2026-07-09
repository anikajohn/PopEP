library(data.table)
library(tidyverse)
library(ggplot2)
library(pals)
library(yaml)
require(optparse)

option_list = list(
  make_option(c("-i", "--index_file"), action="store", default=NA, type='character',
              help="Path to coverage index_file.tsv file "),
  make_option(c("-k", "--count_file"), action="store", default=NA, type='character',
              help="Path to combined_amp_read_versions_counts.tsv file"),
  make_option(c("-f", "--output_file"), action="store", 
              default="index_file_processed.tsv", type='character', 
              help="Path to file used to generate plots"),
  make_option(c("-w", "--output_plot"), action="store", 
              default="plot_index.png", type='character',
              help="Path output plot (subclades by week)"),
  make_option(c("-a", "--amplicon_only"), action="store_true", default=FALSE,
              help="If TRUE, facet only by amplicon (default: FALSE, facet by location and amplicon)"),
  make_option(c("-c", "--colour_palette"), action="store", default=NA, type='character',
              help="path to json file with colour palette")
)



opt = parse_args(OptionParser(option_list=option_list))



######INPUT########

index_file = opt$index_file

count_file = opt$count_file

palette_file = opt$colour_palette



######OUTPUT########
#this would need to be part of rule all

output_plot = opt$output_plot

output_file = opt$output_file

amplicon_only = opt$amplicon_only


######PROCESSING######


df_count = fread(count_file, sep="\t", header=TRUE)
df_count = df_count[df_count$read_name != "read_name", ]
df_count$count = as.numeric(df_count$count)


if (nrow(df_count) == 0L){ 
  
  print("No reads available, empty plot outputted.")
  
  empty_plot <- ggplot() + theme_void()
  ggsave(output_plot, empty_plot, height = 6, width = 9)
  write.table(NULL, output_file, row.names = FALSE, col.names = FALSE)
  
}else{

  dt = fread(index_file)
  
  dt = merge(dt, df_count, by.x = c("seqName"), 
             by.y = c("read_name"))
  

  
  ##### PLOTTING ######

  get_colors <- function(n) {
    # polychrome is capped at 36; interpolate beyond that
    if (n <= 36) {
      as.vector(pals::polychrome(max(2L, n)))[seq_len(n)]
    } else {
      colorRampPalette(as.vector(pals::polychrome(36)))(n)
    }
  }

  if (is.na(palette_file)){
    all_subclades <- unique(dt$subclade)
    custom_palette <- setNames(get_colors(length(all_subclades)), all_subclades)

  }else{

    config <- read_yaml(palette_file)
    defined_palette <- unlist(config)
    all_subclades <- unique( dt$subclade)
    missing_subclades <- setdiff(all_subclades, names(defined_palette))

    if (length(missing_subclades) > 0) {
      random_palette <- setNames(get_colors(length(missing_subclades)), missing_subclades)
      custom_palette<- c(defined_palette, random_palette)

      message("Warning: Assigned random colors to missing subclades: ",
              paste(missing_subclades, collapse = ", "))
    } else {
      custom_palette <- defined_palette
    }
  }
  
  
  
  dt <- dt %>% mutate(date = lubridate::date_decimal(time))

  p <- ggplot(dt, aes(x = date, y = indicies, colour = subclade)) +
    geom_point(alpha = 0.4) +
    geom_smooth(aes(x = date, y = indicies, colour = subclade, weight = count),
                method = "gam", se = TRUE) +
    scale_colour_manual(values = custom_palette) +
    facet_grid(cols = vars(time_window),
                rows = vars(timescale),
                labeller = label_both) +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1))

  
  ##### OUTPUTTING ######
  
  
  ggsave(output_plot,p,height = 6,width = 9)
  write.table(dt,output_file, sep = "\t", row.names = FALSE)
  
}
