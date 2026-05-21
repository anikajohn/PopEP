library(data.table)
library(tidyverse)
library(ggplot2)
library(pals)

source("/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/scripts/THD_functions.R")

######INPUT########
result_dir = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25"

#will be handled by wildcards
amp = "amplicon6"
location = "over_all_locations" #one run also overall 

palette_file = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/colour_palette/H1/colours.json"

# palette_file = NA


nextclade_meta_file= paste(result_dir,
                      location,
                      amp,
                      # "nextclade_results",
                      "combined_nextclade_metadata.tsv",
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
                    "subclade_abundance",
                    "subclade_abundance.png",
                    sep = "/")


output_plot = paste(result_dir,
                    location,
                    amp,
                    "subclade_abundance",
                    "subclade_abundance.tsv",
                    sep = "/")



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
  
  palette_colors <- as.vector(pals::cols25(14))
  subclades <- c("C.1.9", "C.1.9.1", "C.1.9.2", "C.1.9.3",
                 "C.1.9.4","D","D.1","D.3.1","D.5",
                 "C.1","C.1.1","C.1.8","D.2","D.3")
  custom_palette <- setNames(palette_colors, subclades)
  
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
  
  
  ggsave(output_plot_day, by_day, height = 6, width = 9)
  ggsave(output_plot_week, by_week, height = 6, width = 9)
  write.table(dt, output_file, sep = "\t", row.names = FALSE)
  
}
