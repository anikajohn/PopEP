library(tidyverse)
library(data.table)
library(yaml)

setwd("/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/scripts/")

#H1N1 24/25

amp2 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon2/subclade_abundance/combined_nextclade_metadata_processed.tsv")
amp3 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon3/subclade_abundance/combined_nextclade_metadata_processed.tsv")
amp4 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon4/2000n_downsampled/subclade_abundance/nextclade_metadata_processed.tsv")
amp5 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon5/subclade_abundance/combined_nextclade_metadata_processed.tsv")
amp6 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon6/subclade_abundance/combined_nextclade_metadata_processed.tsv")


#H1N1 25/26

amp3 = fread("../results_H1_25_26_updated_ref_dataset/over_all_locations/amplicon3/subclade_abundance/combined_nextclade_metadata_processed.tsv")
amp4 = fread("../results_H1_25_26_updated_ref_dataset/over_all_locations/amplicon4/2000n_downsampled/subclade_abundance/nextclade_metadata_processed.tsv")
amp5 = fread("../results_H1_25_26_updated_ref_dataset/over_all_locations/amplicon5/subclade_abundance/combined_nextclade_metadata_processed.tsv")
amp6 = fread("../results_H1_25_26_updated_ref_dataset/over_all_locations/amplicon6/subclade_abundance/combined_nextclade_metadata_processed.tsv")




amplicon_only = TRUE



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
dt_all = rbind(amp3,amp4,amp5,amp6)
dt = dt_all


by_day = ggplot(dt, aes(fill=subclade, x=date )) + 
  geom_bar(position="stack")+
  scale_fill_manual(values = custom_palette) +
  {if(amplicon_only) facet_grid(cols = vars(amplicon)) else facet_grid(vars(location), vars(amplicon))} +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))
by_day

by_week = ggplot(dt, aes(fill=subclade, x=week_date )) + 
  geom_bar(position="stack")+
  scale_fill_manual(values = custom_palette) +
  {if(amplicon_only) facet_grid(cols = vars(amplicon)) else facet_grid(vars(location), vars(amplicon))} +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))

by_week




