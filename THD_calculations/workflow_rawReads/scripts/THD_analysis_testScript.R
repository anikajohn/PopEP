library(tidyverse)
library(data.table)
source("/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/scripts/THD_functions.R")

######INPUT########
result_dir = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_H1_25_26"

#will be handled by wildcards
amp = "amplicon5/updated_ref_dataset"
location = "over_all_locations" #one run also overall 

palette_file = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/colour_palette/H1/colours.json"

# palette_file = NA

seq_file = paste(result_dir,
                  location,
                  amp,
                  "combined_more_downsampled_amp_reads.fasta",
                 sep = "/")

metadata_file = paste(result_dir,
                  location,
                  amp,
                  "combined_metadata.tsv",
                  sep = "/")

nextclade_file= paste(result_dir,
                      location,
                      amp,
                      # "nextclade_results",
                      "nextclade_metadata.tsv",
                      sep = "/")


#parameters, need to be config input
mutation_rate = 0.00329
# timescale_values <- c(0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
timescale_values <- c(0.15, 0.2)
# timewindow_values <- c(0.019, 0.038, 0.0833, 0.15, 0.23)
timewindow_values <- c(0.019)



######OUTPUT########
#this would need to be part of rule all
output_file = paste(result_dir,
                     location,
                     amp,
                     "THD_calculation",
                     "index_time_parameters.tsv",
                    sep = "/")

output_readIDs = paste(result_dir,
                      location,
                      amp,
                      "THD_calculation",
                      "usedSeq_IDs.tsv",
                      sep = "/")


output_plot = paste(result_dir,
                    location,
                    amp,
                    "THD_calculation",
                    "plot_index_time_parameters.png",
                    sep = "/")




######CALCULATION########

am_seq = readDNAStringSet(seq_file)

if (file.size(seq_file) == 0L){ 
  
  print("No reads available, empty files outputted.")
  
  write.table(NULL, output_file, row.names = FALSE, col.names = FALSE)
  write.table(NULL,output_readIDs, row.names = FALSE, col.names = FALSE)
  
  
}else{
  
  
  metadata = read.table(metadata_file,header = TRUE, sep = "\t",
                        col.names = c("ID","time")) %>% 
    mutate(time = lubridate::decimal_date(as.Date(time, "%Y-%m-%d")))
  
  nextclade = fread(nextclade_file, header = TRUE, sep = "\t", 
                    na.strings = c("", "NA", "N/A"))
  
  
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
  p
  
  # write.table(df_all_results, output_file, sep = "\t")
  # write.table(data.frame(usedIDs = names(am_seq_filt)),
  #             output_readIDs, sep = "\t", row.names = FALSE)
  # 
  # ggsave(output_plot,p,height = 12,width = 18)
  # 
  
  }

