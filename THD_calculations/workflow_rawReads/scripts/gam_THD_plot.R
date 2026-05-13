library(tidyverse)
library(data.table)
library(yaml)

setwd("/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/scripts/")

#H1N1
# amp2 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon2/combined_index_time_parameters.tsv")
# amp3 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon3/combined_index_time_parameters.tsv")
# amp4 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon4/2000n_downsampled/THD_calculation/index_time_parameters.tsv")
# amp5 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon5/combined_index_time_parameters.tsv")
# amp6 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon6/combined_index_time_parameters.tsv")

amp2 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon2/subclade_abundance/")
amp3 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon3/combined_index_time_parameters.tsv")
amp4 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon4/2000n_downsampled/THD_calculation/index_time_parameters.tsv")
amp5 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon5/combined_index_time_parameters.tsv")
amp6 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon6/combined_index_time_parameters.tsv")


#H1N1
# amp2 = fread("../results_H1_25_26_updated_ref_dataset/over_all_locations/amplicon2/combined_index_time_parameters.tsv")
# amp3 = fread("../results_H1_25_26_updated_ref_dataset/over_all_locations/amplicon3/combined_index_time_parameters.tsv")
# amp4 = fread("../results_H1_25_26_updated_ref_dataset/over_all_locations/amplicon4/2000n_downsampled/THD_calculation/index_time_parameters.tsv")
# amp5 = fread("../results_H1_25_26_updated_ref_dataset/over_all_locations/amplicon5/combined_index_time_parameters.tsv")
# amp6 = fread("../results_H1_25_26_updated_ref_dataset/over_all_locations/amplicon6/combined_index_time_parameters.tsv")



######INPUT########

result_dir = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_H1_25_26_updated_ref_dataset"
amp = "amplicon6"
location = "over_all_locations" #one run also overall 

palette_file = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/colour_palette/H1/colours.json"

input_file = paste(result_dir,
                    location,
                    amp,
                    "combined_index_time_parameters.tsv",
                    sep = "/")


######OUTPUT########


output_plot = paste(result_dir,
                    location,
                    amp,
                    "GAMplot_index_time_parameters.png",
                    sep = "/")


######PLOTTING########

if (is.na(palette_file)){
  all_subclades <- unique(dt$subclade)
  random_colors <- as.vector(pals::polychrome(length(all_subclades)))
  custom_palette <- setNames(random_colors, all_subclades)
  
}else{
  
  config <- read_yaml(palette_file)
  defined_palette <- unlist(config)
  all_subclades <- unique(dt$subclade)
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


# dt = fread(input_file)
dt_all = rbind(amp2,amp3,amp4,amp5,amp6)
dt = dt_all[time_window == 0.038 & timescale == 0.35] %>% 
  mutate(amplicon = str_extract(seqName, "amplicon.*"))



df_plot = dt %>% 
  mutate(
    date = lubridate::date_decimal(time),
    week_date = floor_date(as.Date(date), unit = "week"))


p = ggplot(df_plot,
           aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam")+
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon))+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))

p



ggsave(output_plot,p,height = 12,width = 18)
