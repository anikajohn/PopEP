library(tidyverse)
library(data.table)
library(yaml)


setwd("/Users/anjohn/Desktop/THD_calculations_local")

amp2 = fread("nextclade/results_gisaid_2526_AMP2_2000/index_time_parameters_test.tsv")
amp3 = fread("nextclade/results_gisaid_2526_AMP2_2000/index_time_parameters_test.tsv")
amp4 = fread("nextclade/results_gisaid_2526_AMP2_2000/index_time_parameters_test.tsv")
amp5 = fread("nextclade/results_gisaid_2526_AMP2_2000/index_time_parameters_test.tsv")
amp6 = fread("nextclade/results_gisaid_2526_AMP2_2000/index_time_parameters_test.tsv")

palette_file = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/colour_palette/H1/colours.json"

amp2$amplicon = "amplicon2"
amp3$amplicon = "amplicon3"
amp4$amplicon = "amplicon4"
amp5$amplicon = "amplicon5"
amp6$amplicon = "amplicon6"

#want same time span and parameters as ww data
df_all_results = rbind(amp2,amp3,amp4,amp5,amp6) %>% 
  mutate(date = lubridate::date_decimal(time)) %>% 
  dplyr::filter(date >= as.Date("2025-11-01")) %>%
  dplyr::filter(date < as.Date("2026-04-01")) %>% 
  dplyr::filter(time_window == 0.038) %>% 
  dplyr::filter(timescale == 0.35)
  

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
    week_date = floor_date(as.Date(date), unit = "week"))

###INDEX DYNAMIC#####

p = ggplot(df_plot,
           aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam")+
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon))+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))

###SUBCLADE DISTRIBUTION#####

by_day = ggplot(df_plot, aes(fill=subclade, x=date )) + 
  geom_bar(position="stack")+
  scale_fill_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon)) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))

by_week = ggplot(df_plot, aes(fill=subclade, x=week_date )) + 
  geom_bar(position="stack")+
  scale_fill_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon)) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))


##### OUTPUTTING ######
p
by_day 

by_week

