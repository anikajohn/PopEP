library(data.table)
library(tidyverse)
library(ggplot2)
library(pals)
require(optparse)


#This coverage contains all reads as ouputted by the extraction filtering pipeline
#does not account for any filtering done within THD pipeline
#for example, downsample again at start of THD, but here original downsampled 

option_list = list(
  make_option(c("-c", "--coverage_file"), action="store", default=NA, type='character',
              help="Path to coverage summary file (samtools)"),
  make_option(c("-f", "--output_file"), action="store", 
              default="coverage_processed.tsv", type='character', 
              help="Path to file used to generate plots"),
  make_option(c("-b", "--output_plot_covbases"), action="store", 
              default="plot_covbases.png", type='character', 
              help="Path output plot (covbases by date)"),
  make_option(c("-n", "--output_plot_nreads"), action="store", 
              default="plot_numreads.png", type='character',
              help="Path output plot (nreads by date)")
  ,make_option(c("-a", "--amplicon_only"), action="store_true", default=FALSE,
               help="If TRUE, facet only by amplicon (default: FALSE, facet by location and amplicon)")
)



opt = parse_args(OptionParser(option_list=option_list))


######INPUT########

cov_file = opt$coverage_file


######OUTPUT########
#this would need to be part of rule all

output_plot_numreads = opt$output_plot_nreads

output_plot_covbases = opt$output_plot_covbases

output_file_covbases = opt$output_file

amplicon_only = opt$amplicon_only


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
  {if (amplicon_only) facet_grid(cols = vars(amplicon)) else facet_grid(rows = vars(location), cols = vars(amplicon))} +
  theme(text = element_text(size = 20))


p_covbases = ggplot(covdata, aes(x=date, y=covbases)) +
  geom_point(alpha = 0.6)+
  geom_line() +
  {if (amplicon_only) facet_grid(cols = vars(amplicon)) else facet_grid(rows = vars(location), cols = vars(amplicon))} +
  theme(text = element_text(size = 20))


##### OUTPUTTING ######

ggsave(output_plot_numreads,p_numreads,height = 6,width = 9)
ggsave(output_plot_covbases,p_covbases,height = 6,width = 9)
write.table(covdata,output_file_covbases, sep = "\t", row.names = FALSE)





