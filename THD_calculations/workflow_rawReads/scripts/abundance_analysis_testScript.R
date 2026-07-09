library(data.table)
library(tidyverse)
library(ggplot2)
library(pals)
library(yaml)

######INPUT########
result_dir = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_SarsCov2_21_22"

#will be handled by wildcards
amp = "amplicon78"
location = "over_all_locations" #one run also overall

palette_file = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/colour_palette/SarsCov2/colours.json"

# palette_file = NA

subclade_col = "partiallyAliased"
clade_col = "clade_who"

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

output_plot_day = paste(result_dir,
                    location,
                    amp,
                    "subclade_abundance",
                    "subclade_abundance_byDay.png",
                    sep = "/")

output_plot_week = paste(result_dir,
                    location,
                    amp,
                    "subclade_abundance",
                    "subclade_abundance_byWeek.png",
                    sep = "/")

output_file = paste(result_dir,
                    location,
                    amp,
                    "subclade_abundance",
                    "nextclade_metadata_processed.tsv",
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
  # write.table(NULL, output_file, row.names = FALSE, col.names = FALSE)

}else{

  dt = fread(nextclade_meta_file,
             select = c("strain", clade_col, subclade_col, "date", "amplicon", "location"))
  if (subclade_col != "subclade"){setnames(dt, subclade_col, "subclade")}
    
  if (clade_col != "legacy-clade"){setnames(dt, clade_col, "legacy-clade")}
    
  dt = dt %>%
    mutate(`legacy-clade` = as.character(`legacy-clade`),
           week_date = floor_date(as.Date(date), unit = "week"),
           month_date = floor_date(as.Date(date), unit = "month"))

  dt = merge(dt, df_count, by.x = c("strain", "amplicon", "location"),
             by.y = c("read_name", "amplicon", "location"))
  
  dt = dt %>% mutate(subclade_simple = case_when(str_starts(subclade,"BA.1") ~ "BA.1",
                                       str_starts(subclade,"BA.2") ~ "BA.2",
                                       str_starts(subclade,"BA.4") ~ "BA.2",
                                       str_starts(subclade,"BA.5") ~ "BA.2",
                                       str_starts(subclade,"XBB") ~ "BA.2",
                                       str_starts(subclade,"XEC") ~ "BA.2",
                                       str_starts(subclade,"XDV") ~ "BA.2",
                                       str_starts(subclade,"B.1.617.2") ~ "B.1.617.2",
                                       .default = subclade))
  


  ##### PLOTTING ######

  make_colors <- function(n) {
    if (n == 0L) return(character(0))
    if (n <= 36L) as.vector(pals::polychrome(max(n, 2L)))[seq_len(n)]
    else colorRampPalette(as.vector(pals::polychrome(36L)))(n)
  }

  if (is.na(palette_file)){
    all_subclades <- unique(dt$subclade)
    custom_palette <- setNames(make_colors(length(all_subclades)), all_subclades)

  }else{

    config <- read_yaml(palette_file)
    defined_palette <- unlist(config)
    all_subclades <- unique(dt$subclade)
    missing_subclades <- setdiff(all_subclades, names(defined_palette))

    if (length(missing_subclades) > 0) {
      random_palette <- setNames(make_colors(length(missing_subclades)), missing_subclades)
      custom_palette<- c(defined_palette, random_palette)

      message("Warning: Assigned random colors to missing subclades: ",
              paste(missing_subclades, collapse = ", "))
    } else {
      custom_palette <- defined_palette
    }
  }

  dt_day <- dt %>% 
    filter(count > 1) %>% 
    group_by(date, subclade_simple,`legacy-clade`, amplicon, location) %>%
    summarise(count_sum = sum(count, na.rm = TRUE), .groups = "drop")

  dt_week <- dt %>%
    filter(count > 1) %>% 
    group_by(week_date, subclade_simple,`legacy-clade`, amplicon, location) %>%
    summarise(count_sum = sum(count, na.rm = TRUE), .groups = "drop")

  by_day = ggplot(dt_day %>% 
                    filter(count_sum > 1) %>%
                    filter(subclade_simple %in% c("BA.1", "BA.2", "B.1.617.2")), aes(fill=subclade_simple, x=date, y=count_sum)) +
    geom_col(position="stack")+
    # scale_fill_manual(values = custom_palette) +
    {if(amplicon_only) facet_grid(cols = vars(amplicon)) else facet_grid(vars(location), vars(amplicon))} +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1))

  by_week = ggplot(dt_week %>% 
                     filter(count_sum > 1) %>%
                     filter(subclade_simple %in% c("BA.1", "BA.2", "B.1.617.2")), aes(fill=subclade_simple, x=week_date, y=count_sum)) +
    geom_col(position="stack")+
    # scale_fill_manual(values = custom_palette) +
    {if(amplicon_only) facet_grid(cols = vars(amplicon)) else facet_grid(vars(location), vars(amplicon))} +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  by_week
  by_day

  ##### OUTPUTTING ######

  # ggsave(output_plot_day, by_day, height = 6, width = 9)
  # ggsave(output_plot_week, by_week, height = 6, width = 9)
  # write.table(dt, output_file, sep = "\t", row.names = FALSE)

}

