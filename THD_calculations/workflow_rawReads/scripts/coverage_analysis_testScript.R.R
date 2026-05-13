library(data.table)
library(tidyverse)
library(ggplot2)
library(pals)

#This coverage contains all reads as ouputted by the extraction filtering pipeline
#does not account for any filtering done within THD pipeline
#for example, downsample again at start of THD, but here original downsampled 

######INPUT########
result_dir = "/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_H1_24_25_updated_ref_dataset"

#will be handled by wildcards
amp = ""
location = "over_all_locations" #one run also overall 

cov_file = paste(result_dir,
                            location,
                            amp,
                            "H1_2425.tsv",
                            sep = "/")



######OUTPUT########
#this would need to be part of rule all

output_plot_numreads = paste(result_dir,
                        location,
                        amp,
                        "plot_numreads.png",
                        sep = "/")

output_plot_covbases = paste(result_dir,
                         location,
                         amp,
                         "plot_covbases.png",
                         sep = "/")

output_file_covbases = paste(result_dir,
                             location,
                             amp,
                             "coverage_processed.tsv",
                             sep = "/")



######PROCESSING######


covdata <- fread(cov_file, col.names = c("sample","amplicon", "name","start","end",
                                         "numreads","covbases", "coverage",
                                         "meandepth","meanbaseq","meanmapq",
                                         "location","amplicon_2")) %>% 
  mutate(date = as.Date(str_extract(sample, "\\d{4}_\\d{2}_\\d{2}"),"%Y_%m_%d")) %>% 
  separate_wider_delim(sample, delim = "/", names = c("sample", "batch"))


##### PLOTTING ######

p_numreads = ggplot(covdata, aes(x=date, y=numreads)) +
  geom_point(alpha = 0.6)+
  geom_line() +
  facet_grid(rows = vars(location), cols = vars(amplicon))+
  theme(text = element_text(size = 20))


p_covbases = ggplot(covdata, aes(x=date, y=covbases)) +
  geom_point(alpha = 0.6)+
  geom_line() +
  facet_grid(rows = vars(location), cols = vars(amplicon))+
  theme(text = element_text(size = 20))


##### OUTPUTTING ######

ggsave(output_plot_numreads,p_numreads,height = 6,width = 9)
ggsave(output_plot_covbases,p_covbases,height = 6,width = 9)
write.table(dt,output_file_covbases, sep = "\t", row.names = FALSE)





