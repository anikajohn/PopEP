library(tidyverse)
library(data.table)
library(yaml)
library(pals)
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "THD_functions.R"))
require(optparse)



option_list = list(
  make_option(c("-q", "--seq_file"), action="store", default=NA, type='character',
              help="Path to fasta file with sequences of interest"),
  make_option(c("-m", "--meta_file"), action="store", default=NA, type='character',
              help="Path to metadata file of sequences of interest"),
  make_option(c("-n", "--nextclade_file"), action="store", default=NA, type='character',
              help="Path to nextclade tsv file of sequences of interest"),
  make_option(c("-d", "--distance_mat_file"), action="store", default=NA, type='character',
              help="Path to precomputed distance matrix CSV file"),
  make_option(c("-k", "--count_file"), action="store", default=NA, type='character',
              help="Path to combined amp read versions counts TSV file"),
  make_option(c("-r", "--mutation_rate"), action="store", default=NA, type="numeric",
              help="mutation rate of given segment"),
  make_option(c("-s", "--timescales"), type = "character", 
              default = "0.15,0.2,0.25,0.3,0.35,0.4", 
              help = "Comma-separated list of timescales",
              metavar = "list"),
  make_option(c("-w", "--timewindows"), type = "character", 
              default = "0.019, 0.038, 0.0833, 0.15, 0.23", 
              help = "Comma-separated list of timewindows",
              metavar = "list"),
  make_option(c("-f", "--output_file"), action="store", 
              default="index_time_parameters.tsv", type='character', 
              help="Path to file used to generate plots"),
  make_option(c("-c", "--chunk_size"), action="store",
              default=1000L, type='integer',
              help="how many lines of distance matrix read in at once"),
  make_option(c("-b", "--subclade_col"), action="store", default="subclade", type='character',
              help="Column name for subclade in nextclade TSV (e.g. 'subclade' or 'partiallyAliased')"),
  make_option(c("-l", "--clade_col"), action="store", default="legacy-clade", type='character',
              help="Column name for legacy/WHO clade in nextclade TSV (e.g. 'legacy-clade' or 'clade_who')")
)



opt = parse_args(OptionParser(option_list=option_list))

## ---- Cluster parallelism setup ----
# SLURM array task index (1-based); falls back to 1 for local runs
slurm_task_id  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
# Number of cores available on this node for within-job parallelism
mc_cores       <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

######INPUT########

seq_file = opt$seq_file
metadata_file = opt$meta_file
nextclade_file= opt$nextclade_file
distance_mat_file = opt$distance_mat_file
count_file = opt$count_file
palette_file = opt$colour_palette
chunk_size = as.integer(opt$chunk_size)


#parameters, need to be config input
mutation_rate = opt$mutation_rate
timescale_values <- as.numeric(strsplit(opt$timescales, ",")[[1]])
timewindow_values <- as.numeric(strsplit(opt$timewindows, ",")[[1]])



######OUTPUT########
#this would need to be part of rule all
output_file = opt$output_file





######CALCULATION########

am_seq = readDNAStringSet(seq_file, nrec = 1)

if (file.size(seq_file) == 0L){ 
  
  print("No reads available, empty files outputted.")
  
  write.table(NULL, output_file, row.names = FALSE, col.names = FALSE)
  
  
}else{
  
  metadata = read.table(metadata_file,header = TRUE, sep = "\t",
                        col.names = c("ID","time")) %>% 
    mutate(time = lubridate::decimal_date(as.Date(time, "%Y-%m-%d")))
  
  nextclade = fread(nextclade_file, header = TRUE, sep = "\t",
                    na.strings = c("", "NA", "N/A"))
  
  if (opt$subclade_col != "subclade" && opt$subclade_col %in% names(nextclade))
    setnames(nextclade, opt$subclade_col, "subclade")
  if (opt$clade_col != "legacy-clade" && opt$clade_col %in% names(nextclade))
    setnames(nextclade, opt$clade_col, "legacy-clade")

  genome_length = width(am_seq)
  
  
  # Create all combinations
  param_grid <- expand.grid(timescale = timescale_values, 
                            time_window = timewindow_values)
  
  # ## ---- SLURM array dispatch ----
  # ## If running as a SLURM array job (SLURM_ARRAY_TASK_ID is set), each task
  # ## processes exactly one row of param_grid and writes its own output file.
  # ## When running locally (slurm_task_id == 0), all combinations run sequentially.
  # if (slurm_task_id > 0) {
  #   if (slurm_task_id > nrow(param_grid)) {
  #     stop(sprintf("SLURM_ARRAY_TASK_ID (%d) exceeds number of parameter combinations (%d).",
  #                  slurm_task_id, nrow(param_grid)))
  #   }
  #   param_grid <- param_grid[slurm_task_id, , drop = FALSE]
  #   ## Append task id to output paths so array tasks don't clobber each other
  #   output_file   <- sub("(\\.tsv)?$", sprintf("_%02d.tsv", slurm_task_id), output_file, perl = TRUE)
  #   message(sprintf("Array task %d / %d: timescale = %.2f, time_window = %.3f",
  #                   slurm_task_id, nrow(expand.grid(timescale = timescale_values, time_window = timewindow_values)),
  #                   param_grid$timescale, param_grid$time_window))
  # }

  df_all_results <- pmap_dfr(param_grid, function(timescale, time_window) {
    message(sprintf("Running: timescale = %.2f, time_window = %.3f", timescale, time_window))
    
    df_temp <- compute.merge.haplo.results(
      dist_file = distance_mat_file,
      time_window = time_window,
      metadata = metadata,
      mutation_rate = mutation_rate,
      timescale = timescale,
      genome_length = genome_length,
      nextclade = nextclade,
      count_file = count_file,
      chunk_size = chunk_size
    )
    
    df_temp$timescale <- timescale
    df_temp$time_window <- time_window
    
    return(df_temp)
  })
  
  # View summary
  message(sprintf("Total combinations: %d", nrow(param_grid)))
  message(sprintf("Total rows in combined dataframe: %d", nrow(df_all_results)))

  ##### OUTPUTTING ######

  write.table(df_all_results, output_file, sep = "\t")
  
}
