library(tidyverse)
library(data.table)
library(yaml)
library(pals)
source("scripts/THD_functions.R")
require(optparse)



option_list = list(
  make_option(c("-q", "--seq_file"), action="store", default=NA, type='character',
              help="Path to fasta file with sequences of interest"),
  make_option(c("-m", "--meta_file"), action="store", default=NA, type='character',
              help="Path to metadata file of sequences of interest"),
  make_option(c("-n", "--nextclade_file"), action="store", default=NA, type='character',
              help="Path to nextclade tsv file of sequences of interest"),
  make_option(c("-r", "--mutation_rate"), action="store", default=NA, type="numeric",
              help="mutation rate of given segment"),
  make_option(c("-s", "--timescales"), type = "character", 
              default = "0.15,0.2,0.25,0.3,0.35,0.4", 
              help = "Comma-separated list of timescales",
              metavar = "list"),
  make_option(c("-w", "--timewindows"), type = "character", 
              default = "0.019, 0.038, 0.0833, 0.15, 0.23", 
              help = "Comma-separated list of timewindows",
              metavar = "list"),
  make_option(c("-f", "--output_file"), action="store", 
              default="index_time_parameters.tsv", type='character', 
              help="Path to file used to generate plots"),
  make_option(c("-i", "--output_IDs"), action="store", default="usedSeq_IDs.tsv", 
              type='character', help="Path to IDs of sequences used for plot"),
  make_option(c("-p", "--output_plot"), action="store", 
              default="plot_index_time_parameters.png", type='character',
              help="Path output plot"),
  make_option(c("-c", "--colour_palette"), action="store", default=NA, type='character',
              help="path to json file with colour palette")
  
)



opt = parse_args(OptionParser(option_list=option_list))

######INPUT########

seq_file = opt$seq_file
metadata_file = opt$meta_file
nextclade_file= opt$nextclade_file
palette_file = opt$colour_palette


#parameters, need to be config input
mutation_rate = opt$mutation_rate
timescale_values <- as.numeric(strsplit(opt$timescales, ",")[[1]])
timewindow_values <- as.numeric(strsplit(opt$timewindows, ",")[[1]])



######OUTPUT########
#this would need to be part of rule all
output_file = opt$output_file
output_readIDs = opt$output_IDs


output_plot = opt$output_plot




######CALCULATION########

am_seq = readDNAStringSet(seq_file)

if (file.size(seq_file) == 0L){ 
  
  print("No reads available, empty files outputted.")
  
  write.table(NULL, output_file, row.names = FALSE, col.names = FALSE)
  write.table(NULL,output_readIDs, row.names = FALSE, col.names = FALSE)
  file.create(output_plot)
  
  
}else{
  
  metadata = read.table(metadata_file,header = TRUE, sep = "\t",
                        col.names = c("ID","time")) %>% 
    mutate(time = lubridate::decimal_date(as.Date(time, "%Y-%m-%d")))
  
  nextclade = fread(nextclade_file, header = TRUE, sep = "\t")
  
  
  length_freq_table = table(width(am_seq))
  most_freq_length = as.numeric(names(sort(length_freq_table, 
                                           decreasing = TRUE)[1]))
  
  #strong simplification cause not sure if THD can handle multiple seq lengths 
  #for now
  am_seq_filt = am_seq[width(am_seq) == most_freq_length]
  
  dist_matrix <- pwalign::stringDist(am_seq_filt, method = "hamming")
  distance_mat = as.matrix(dist_matrix)
  
  genome_length = most_freq_length
  
  
  # Create all combinations
  param_grid <- expand.grid(timescale = timescale_values, 
                            time_window = timewindow_values)
  
  df_all_results <- pmap_dfr(param_grid, function(timescale, time_window) {
    message(sprintf("Running: timescale = %.2f, time_window = %.3f", timescale, time_window))
    
    df_temp <- compute_and_merge_results(
      distance_mat = distance_mat,
      time_window = time_window,
      metadata = metadata,
      mutation_rate = mutation_rate,
      timescale = timescale,
      genome_length = genome_length,
      nextclade = nextclade
    )
    
    df_temp$timescale <- timescale
    df_temp$time_window <- time_window
    
    return(df_temp)
  })
  
  # View summary
  message(sprintf("Total combinations: %d", nrow(param_grid)))
  message(sprintf("Total rows in combined dataframe: %d", nrow(df_all_results)))
  
  
  
  ##### PLOTTING ######
  
  if (is.na(palette_file)){
    all_subclades <- unique(df_all_results$subclade)
    random_colors <- as.vector(pals::polychrome(length(all_subclades)))
    custom_palette <- setNames(random_colors, all_subclades)
    
  }else{
    
    config <- read_yaml(palette_file)
    defined_palette <- unlist(config)
    all_subclades <- unique(df_all_results$subclade)
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
  

  df_plot = df_all_results %>% 
    mutate(
      date = lubridate::date_decimal(time),
      week_date = floor_date(as.Date(date), unit = "week"))
  
  p = ggplot(df_plot,
         aes(x = date, y = indicies, colour = subclade)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE)+
    scale_colour_manual(values = custom_palette) +
    facet_grid(cols = vars(time_window),
               rows = vars(timescale),
               labeller = label_both)+
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  ##### OUTPUTTING ######
  
  write.table(df_all_results, output_file, sep = "\t")
  write.table(data.frame(usedIDs = names(am_seq_filt)),
              output_readIDs, sep = "\t", row.names = FALSE)
  
  ggsave(output_plot,p,height = 12,width = 18)
  
}
