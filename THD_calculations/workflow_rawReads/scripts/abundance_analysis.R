library(data.table)
library(tidyverse)
library(ggplot2)
library(pals)
library(yaml)
require(optparse)

option_list = list(
  make_option(c("-m", "--nextclade_meta_file"), action="store", default=NA, type='character',
              help="Path to coverage nextclade_metadata.tsv file (see rule merge_nextclade_meta_ww)"),
  make_option(c("-f", "--output_file"), action="store", 
              default="nextclade_metadata_processed.tsv", type='character', 
              help="Path to file used to generate plots"),
  make_option(c("-d", "--output_plot_byDay"), action="store", 
              default="plot_subclade_frequency_byDay.png", type='character', 
              help="Path output plot (subclades by date)"),
  make_option(c("-w", "--output_plot_byWeek"), action="store", 
              default="plot_subclade_frequency_byWeek.png", type='character',
              help="Path output plot (subclades by week)"),
  make_option(c("-a", "--amplicon_only"), action="store_true", default=FALSE,
              help="If TRUE, facet only by amplicon (default: FALSE, facet by location and amplicon)"),
  make_option(c("-k", "--count_file"), action="store", default=NA, type='character',
              help="Path to combined_amp_read_versions_counts.tsv file"),
  make_option(c("-c", "--colour_palette"), action="store", default=NA, type='character',
              help="path to json file with colour palette")
)



opt = parse_args(OptionParser(option_list=option_list))



######INPUT########

nextclade_meta_file = opt$nextclade_meta_file

count_file = opt$count_file

palette_file = opt$colour_palette



######OUTPUT########
#this would need to be part of rule all

output_plot_day = opt$output_plot_byDay

output_plot_week = opt$output_plot_byWeek

output_file = opt$output_file

amplicon_only = opt$amplicon_only


######PROCESSING######


df_count = fread(count_file, sep="\t", header=TRUE)
df_count = df_count[df_count$read_name != "read_name", ]
df_count$count = as.numeric(df_count$count)


if (nrow(df_count) == 0L){ 
  
  print("No reads available, empty plot outputted.")
  
  empty_plot <- ggplot() + theme_void()
  ggsave(output_plot_day, empty_plot, height = 6, width = 9)
  ggsave(output_plot_week, empty_plot, height = 6, width = 9)
  write.table(NULL, output_file, row.names = FALSE, col.names = FALSE)
  
}else{

  dt = fread(nextclade_meta_file, 
             select = c("strain","legacy-clade","subclade", "date", "amplicon", "location")) %>% 
    mutate(`legacy-clade` = as.character(`legacy-clade`),
           week_date = floor_date(as.Date(date), unit = "week"),
           month_date = floor_date(as.Date(date), unit = "month"))
  
  dt = merge(dt, df_count, by.x = c("strain", "amplicon", "location"), 
             by.y = c("read_name", "amplicon", "location"))
  

  
  ##### PLOTTING ######
  
  if (is.na(palette_file)){
    all_subclades <- unique(dt$subclade)
    random_colors <- as.vector(pals::polychrome(length(all_subclades)))
    custom_palette <- setNames(random_colors, all_subclades)
    
  }else{
    
    config <- read_yaml(palette_file)
    defined_palette <- unlist(config)
    all_subclades <- unique( dt$subclade)
    missing_subclades <- setdiff(all_subclades, names(defined_palette))
    
    if (length(missing_subclades) > 0) {
      random_colors <- as.vector(pals::polychrome(length(missing_subclades)))
      random_palette <- setNames(random_colors, missing_subclades)
      custom_palette<- c(defined_palette, random_palette)
      
      message("Warning: Assigned random colors to missing subclades: ", 
              paste(missing_subclades, collapse = ", "))
    } else {
      custom_palette <- defined_palette
    }
  }
  
  
  
  dt_day <- dt %>%
    group_by(date, subclade, amplicon, location) %>%
    summarise(count_sum = sum(count, na.rm = TRUE), .groups = "drop")
  
  dt_week <- dt %>%
    group_by(week_date, subclade, amplicon, location) %>%
    summarise(count_sum = sum(count, na.rm = TRUE), .groups = "drop")
  
  by_day = ggplot(dt_day, aes(fill=subclade, x=date, y=count_sum)) + 
    geom_col(position="stack")+
    scale_fill_manual(values = custom_palette) +
    {if(amplicon_only) facet_grid(cols = vars(amplicon)) else facet_grid(vars(location), vars(amplicon))} +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  by_week = ggplot(dt_week, aes(fill=subclade, x=week_date, y=count_sum)) + 
    geom_col(position="stack")+
    scale_fill_manual(values = custom_palette) +
    {if(amplicon_only) facet_grid(cols = vars(amplicon)) else facet_grid(vars(location), vars(amplicon))} +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ##### OUTPUTTING ######
  
  
  ggsave(output_plot_day,by_day,height = 6,width = 9)
  ggsave(output_plot_week,by_week,height = 6,width = 9)
  write.table(dt,output_file, sep = "\t", row.names = FALSE)
  
}
