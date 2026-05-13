library(data.table)
library(tidyverse)
library(ggplot2)
library(pals)

######INPUT########
result_dir = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_H3_24_25"

#will be handled by wildcards
amp = "amplicon1"
location = "laupen" #one run also overall 

amplicon_only = TRUE

nextclade_meta_file = paste(result_dir,
                      location,
                      amp,
                      "nextclade_metadata.tsv",
                      sep = "/")

readIDs_file = paste(result_dir,
                       location,
                       amp,
                       "THD_claculation",
                       "usedSeq_IDs.tsv",
                       sep = "/")



######OUTPUT########
#this would need to be part of rule all

output_plot_day = paste(result_dir,
                    location,
                    amp,
                    "subclade_abundance",
                    "plot_subclade_frequency_byDay.png",
                    sep = "/")

output_plot_week = paste(result_dir,
                        location,
                        amp,
                        "subclade_abundance",
                        "plot_subclade_frequency_byWeek.png",
                        sep = "/")

output_file = paste(result_dir,
                        location,
                        amp,
                        "subclade_abundance",
                        "nextclade_metadata_processed.tsv",
                        sep = "/")



######PROCESSING######


if (file.size(readIDs_file) == 0L){ 
  
  print("No reads available, empty files outputted.")
  
  write.table(NULL,output_file, row.names = FALSE, col.names = FALSE)
  
  
  
}else{

  readIDs = fread(readIDs_file)
  dt = fread(nextclade_meta_file) %>% 
    filter(seqName %in% as.list(readIDs)$usedIDs) %>% 
    mutate(`legacy-clade` = as.character(`legacy-clade`),
           week_date = floor_date(as.Date(date), unit = "week"),
           month_date = floor_date(as.Date(date), unit = "month"))
  
  
  ##### PLOTTING ######
  
  palette_colors <- as.vector(pals::cols25(14))
  subclades <- c("C.1.9", "C.1.9.1", "C.1.9.2", "C.1.9.3",
                 "C.1.9.4","D","D.1","D.3.1","D.5",
                 "C.1","C.1.1","C.1.8","D.2","D.3")
  custom_palette <- setNames(palette_colors, subclades)
  
  
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
