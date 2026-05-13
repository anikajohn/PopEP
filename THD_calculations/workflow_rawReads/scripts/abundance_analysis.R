library(data.table)
library(tidyverse)
library(ggplot2)
library(pals)
library(yaml)
require(optparse)

option_list = list(
  make_option(c("-m", "--nextclade_meta_file"), action="store", default=NA, type='character',
              help="Path to coverage nextclade_metadata.tsv file (see rule merge_nextclade_meta_ww)"),
  make_option(c("-i", "--THD_IDs"), action="store", default="usedSeq_IDs.tsv", 
              type='character', help="Path to IDs of sequences used in THD_analysis.R"),
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
  make_option(c("-c", "--colour_palette"), action="store", default=NA, type='character',
              help="path to json file with colour palette")
)



opt = parse_args(OptionParser(option_list=option_list))



######INPUT########

nextclade_meta_file = opt$nextclade_meta_file

readIDs_file = opt$THD_IDs

palette_file = opt$colour_palette



######OUTPUT########
#this would need to be part of rule all

output_plot_day = opt$output_plot_byDay

output_plot_week = opt$output_plot_byWeek

output_file = opt$output_file

amplicon_only = opt$amplicon_only


######PROCESSING######

if (file.size(readIDs_file) == 0L){ 
  
  print("No reads available, empty files outputted.")
  
  write.table(NULL,output_file, row.names = FALSE, col.names = FALSE)
  file.create(output_plot_day)
  file.create(output_plot_week)
  
}else{

  readIDs = fread(readIDs_file)
  
  dt = fread(nextclade_meta_file) %>% 
    filter(seqName %in% as.list(readIDs)$usedIDs) %>% 
    mutate(`legacy-clade` = as.character(`legacy-clade`),
           week_date = floor_date(as.Date(date), unit = "week"),
           month_date = floor_date(as.Date(date), unit = "month"))
  
  
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
  
  
  
  by_day = ggplot(dt, aes(fill=subclade, x=date )) + 
    geom_bar(position="stack")+
    scale_fill_manual(values = custom_palette) +
    {if(amplicon_only) facet_grid(cols = vars(amplicon)) else facet_grid(vars(location), vars(amplicon))} +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  by_week = ggplot(dt, aes(fill=subclade, x=week_date )) + 
    geom_bar(position="stack")+
    scale_fill_manual(values = custom_palette) +
    {if(amplicon_only) facet_grid(cols = vars(amplicon)) else facet_grid(vars(location), vars(amplicon))} +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  ##### OUTPUTTING ######
  
  
  ggsave(output_plot_day,by_day,height = 6,width = 9)
  ggsave(output_plot_week,by_week,height = 6,width = 9)
  write.table(dt,output_file, sep = "\t", row.names = FALSE)
  
}
