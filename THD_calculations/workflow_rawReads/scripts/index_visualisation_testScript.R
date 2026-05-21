library(data.table)
library(tidyverse)
library(ggplot2)
library(pals)

source("/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/scripts/THD_functions.R")

######INPUT########
result_dir = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/test_haplotypes/results_haplo_H1_24_25/"

#will be handled by wildcards
amp = "amplicon4/updated_ref_dataset"
location = "geneve" #one run also overall 

palette_file = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/colour_palette/H1/colours.json"

# palette_file = NA


index_file= paste(result_dir,
                      location,
                      amp,
                      "THD_calculation",
                      "index_time_parameters.tsv",
                      sep = "/")

count_file = paste(result_dir,
                   location,
                   amp,
                   "combined_amp_read_versions_counts.tsv",
                   sep = "/")

amplicon_only = TRUE

######OUTPUT########



output_plot = paste(result_dir,
                    location,
                    amp,
                    "THD_calculation",
                    "THD_plot.png",
                    sep = "/")





######PROCESSING######

df_count = fread(count_file, sep="\t", header=TRUE)
df_count = df_count[df_count$read_name != "read_name", ]
df_count$count = as.numeric(df_count$count)


if (nrow(df_count) == 0L){ 
  
  print("No reads available, empty plot outputted.")
  
  empty_plot <- ggplot() + theme_void()
  ggsave(output_plot, empty_plot, height = 6, width = 9)

}else{

  dt = fread(index_file)
  
  dt = merge(dt, df_count, by.x = c("seqName"), 
             by.y = c("read_name"))
  
  
  ##### PLOTTING ######
  
  palette_colors <- as.vector(pals::cols25(14))
  subclades <- c("C.1.9", "C.1.9.1", "C.1.9.2", "C.1.9.3",
                 "C.1.9.4","D","D.1","D.3.1","D.5",
                 "C.1","C.1.1","C.1.8","D.2","D.3")
  custom_palette <- setNames(palette_colors, subclades)
  
  dt <- dt %>%
      mutate(
        date      = lubridate::date_decimal(time))

    p <- ggplot(dt, aes(x = date, y = indicies, colour = subclade)) +
      geom_point(alpha = 0.4) +
      geom_smooth(aes(x = date, y = indicies, colour = subclade, weight = count),
                  method = "gam", se = FALSE) +
      scale_colour_manual(values = custom_palette) +
      facet_grid(cols = vars(time_window),
                 rows = vars(timescale),
                 labeller = label_both) +
      theme(text = element_text(size = 20),
            axis.text.x = element_text(angle = 45, hjust = 1))

  
  ##### OUTPUTTING ######
  
  
  ggsave(output_plot_day, by_day, height = 6, width = 9)
  ggsave(output_plot_week, by_week, height = 6, width = 9)
  write.table(dt, output_file, sep = "\t", row.names = FALSE)
  
}
