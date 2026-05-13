library(tidyr)
library(Biostrings)
library(msa)
library(ggplot2)
library(stringr)
library(tidyverse)
library(pals)


setwd("/Users/anjohn/Desktop/THD_calculations")

#' @title Compute the index bandwidth from a set timescale, genome length and mutation rate 
#' @param timescale Index timescale (in years)
#' @param mutation_rate Mutation rate (in mutation/site/year)
#' @param genome_length Genome length (in bp)
#' @return Index bandwidth
#' @references Wirth, T, Wong, V, Vandenesch, F, Rasigade, J-P. Applied phyloepidemiology: Detecting drivers of pathogen transmission from genomic signatures using density measures. Evol Appl. 2020; 13: 1513– 1525. https://doi.org/10.1111/eva.12991
#' @export
index.bandwidth = function(timescale, genome_length, mutation_rate) {
  if (timescale == 0) return(NA)
  h = genome_length*(1 - exp(-2 * mutation_rate * timescale))
  objfun = function(b) (1/2 * (1 - b^genome_length) - (1 - b^h))^2
  opt = optimise(objfun, c(1e-09, 1 - 1e-09))
  return(opt$minimum)
}


#' Compute the diversity index from a distance matrix (no tree required)
#'
#' Computes the diversity index directly from a pairwise distance matrix between
#' samples, without requiring a phylogenetic tree. This is useful when users have
#' precomputed genetic distances (e.g., from sequence alignments) without a tree.
#'
#' @param distance_mat Pairwise distance matrix (e.g., Hamming distances) between samples
#' @param time_window Time window for the index computation (positive number, in years)
#' @param metadata Data frame with columns \code{ID} and \code{time}
#' @param mutation_rate Mutation rate (in mutation/site/year)
#' @param timescale Index timescale (in years)
#' @param genome_length Genome length (in bp)
#' @return Named vector of index values for each sample
#' @export
compute.index.from.distance = function(distance_mat, time_window, metadata, mutation_rate, timescale, genome_length){
  ## Input validation
  if (!is.matrix(distance_mat)) {
    stop("'distance_mat' must be a matrix.")
  }
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data.frame.")
  }
  required_cols <- c("ID", "time")
  missing_cols <- setdiff(required_cols, names(metadata))
  if (length(missing_cols) > 0) {
    stop("'metadata' is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(time_window) || length(time_window) != 1 || time_window <= 0) {
    stop("'time_window' must be a single positive number.")
  }
  if (!is.numeric(mutation_rate) || length(mutation_rate) != 1 || mutation_rate <= 0) {
    stop("'mutation_rate' must be a single positive number.")
  }
  if (!is.numeric(timescale) || length(timescale) != 1 || timescale <= 0) {
    stop("'timescale' must be a single positive number.")
  }
  if (!is.numeric(genome_length) || length(genome_length) != 1 || genome_length <= 0) {
    stop("'genome_length' must be a single positive integer.")
  }
  
  ## Match metadata to matrix
  sample_ids <- colnames(distance_mat)
  times <- metadata$time[match(sample_ids, metadata$ID)]
  if (any(is.na(times))) {
    stop("Some sample IDs in 'distance_mat' are not found in 'metadata$ID'.")
  }
  n <- length(sample_ids)
  
  ## Time-window filtering and temporal correction
  hamming_mat <- distance_mat
  window_mat <- matrix(TRUE, nrow = n, ncol = n)
  
  for (j in seq_len(n)) {
    t_j <- times[j]
    for (i in seq_len(n)) {
      if (i == j) next
      t_i <- times[i]
      ## Filter: only keep pairs within the time window
      if (abs(t_i - t_j) > time_window) {
        hamming_mat[i, j] <- NA
        window_mat[i, j] <- FALSE
      } else {
        ## Temporal correction: adjust distance by time difference
        hamming_mat[i, j] <- hamming_mat[i, j] + abs(t_i - t_j)
      }
    }
  }
  
  ## Scale distances
  hamming_mat <- hamming_mat * mutation_rate * genome_length / 2
  
  ## Index computation
  bandwidth <- index.bandwidth(timescale, genome_length, mutation_rate)
  B <- bandwidth^hamming_mat
  
  Bsums <- colSums(B, na.rm = TRUE) - 1 #-1 to exclude self-comparisons
  index <- Bsums / (colSums(window_mat) - 1) #averaging over all allowed times 
  names(index) <- sample_ids
  
  return(index)
}




####Amplicon Length #####

sequence <- readDNAStringSet("nextclade/results_Nexstrain_asOF_Mar26/nextclade.aligned.fasta")
# sequence_amp_len <- subseq(sequence, start = 473, end = 917)

metadata = read.table("../Influenza_gisaid/H1N1/Nexstrain_asOF_Mar26/metadata.tsv", 
                      header = TRUE, 
                      col.names = c("ID", "time"),
                      colClasses = c("character", "Date"))
metadata = drop_na(metadata)

metadata = metadata %>% dplyr::filter(!str_detect(ID, "readVersion"))
# Filter sequences to only include those with valid metadata
valid_ids = intersect(names(sequence), metadata$ID)
sequence = sequence[valid_ids]

sequence_amp_len = subseq(sequence, 473,917)


dist_matrix <- pwalign::stringDist(sequence_amp_len, method = "hamming")
distance_mat = as.matrix(dist_matrix)
time_window = 0.25 #from phylowave paper but for H3N2
mutation_rate = 0.00329 #is it same as clock rate used used for augur refine?
timescale = 0.4 #from phylowave paper but for H3N2
genome_length = 445

# sample_ids <- colnames(distance_mat)
# times <- metadata$time[match(sample_ids, metadata$ID)]

indicies = compute.index.from.distance(distance_mat, 
                                       time_window, 
                                       metadata, 
                                       mutation_rate, 
                                       timescale, 
                                       genome_length)

hist(indicies, breaks = 10, main = "Distribution of Diversity Index", 
     xlab = "Diversity Index", ylab = "Frequency")

indicies_df = as.data.frame(indicies)
indicies_df$seqName <- rownames(indicies_df)

nextclade = read.table("nextclade/results_Nexstrain_asOF_Mar26/nextclade.tsv", header = TRUE, sep = "\t")


gr = nextclade[c("index","seqName","subclade")]
# rownames(gr) <- gr$seqName


ggplot(gr, aes(x = subclade)) + 
  geom_bar()

test = merge(gr, indicies_df, by = 'seqName')

test2 = merge(test, metadata, by.x = 'seqName', by.y = 'ID')

test2$week = strftime(test2$time, format = "%V-%Y")

df_summary = test2 %>%
  # Truncate the full date to the start of its week
  dplyr::filter(time >= as.Date("2024-01-01")) %>% 
  mutate(week_date = floor_date(time, unit = "week")) %>% 
  group_by(week_date, subclade) %>% 
  summarise(median_index = median(indicies),
            mean_index = mean(indicies), .groups = "drop",
            n_sequences = n_distinct(seqName))

amp_length_thd <- ggplot(df_summary, 
                          aes(x = week_date, y = median_index, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE)+
  scale_colour_manual(values = as.vector(pals::cols25(20))) +
  scale_x_date(date_labels = "%V-%Y")

amp_length_thd




df_summary_seq = test2 %>%
  # Truncate the full date to the start of its week
  dplyr::filter(time >= as.Date("2024-01-01")) %>% 
  mutate(week_date = floor_date(time, unit = "week")) %>% 
  group_by(week_date) %>% 
  summarise(n_sequences = n_distinct(seqName))


ggplot(df_summary_seq, aes(x =week_date, y =  n_sequences)) + 
  geom_point()

###ONLY CONSIDER FULL LENGTH SEQUENCES FROM NEXTSTRAIN REF #####

sequence <- readDNAStringSet("nextclade/results_Nexstrain_asOF_Mar26/nextclade.aligned.fasta")
# sequence_amp_len <- subseq(sequence, start = 473, end = 917)

metadata = read.table("../Influenza_gisaid/H1N1/Nexstrain_asOF_Mar26/metadata.tsv", 
                      header = TRUE, 
                      col.names = c("ID", "time"),
                      colClasses = c("character", "Date"))
metadata = drop_na(metadata)

metadata = metadata %>% dplyr::filter(!str_detect(ID, "readVersion"))
# Filter sequences to only include those with valid metadata
valid_ids = intersect(names(sequence), metadata$ID)
sequence = sequence[valid_ids]

dist_matrix <- pwalign::stringDist(sequence, method = "hamming")
distance_mat = as.matrix(dist_matrix)
time_window = 0.25 #from phylowave paper but for H3N2
mutation_rate = 0.00329 #is it same as clock rate used used for augur refine?
timescale = 0.4 #from phylowave paper but for H3N2
genome_length = 1752

# sample_ids <- colnames(distance_mat)
# times <- metadata$time[match(sample_ids, metadata$ID)]

indicies = compute.index.from.distance(distance_mat, 
                                       time_window, 
                                       metadata, 
                                       mutation_rate, 
                                       timescale, 
                                       genome_length)

hist(indicies, breaks = 10, 
     main = "Distribution of Diversity Index", 
     xlab = "Diversity Index", ylab = "Frequency")


#######



###CONSIDER FULL LENGTH SEQUENCES FROM ALL GISAID #####

sequence <- readDNAStringSet("../Influenza_gisaid/filtered_sequences_H1_gisaid_2018TO05Jan25.fasta")

metadata = read.table("../Influenza_gisaid/filtered_sequences_H1_gisaid_2018TO05Jan25_metadata.tsv", 
                      header = TRUE,
                      col.names = c("ID", "time"),
                      colClasses = c("character", "Date"))

metadata = drop_na(metadata)

#Want season 23/24 as background season and 24/25 as current
#Goal to see bump of index in later time --> D.3.1 expansion
gisaid_metadata = metadata  %>% 
  dplyr::filter(time >= as.Date("2023-10-01")) 


gisaid_overview_perDay = gisaid_metadata  %>% 
  group_by(time) %>% 
  summarise(n_sequences = n_distinct(ID))



ggplot(gisaid_overview_perDay , aes(x =time, y =  n_sequences)) + 
  geom_point()+
  geom_smooth(method = "loess")+
  labs(x = "Date", y = "# sequences in reference tree")


# Filter sequences to only include those with valid metadata
valid_ids = intersect(names(sequence), gisaid_metadata$ID)
sequence = sequence[valid_ids]

writeXStringSet(sequence, 
                filepath="nextclade/input_H1_gisaid_2324_2425/sequence.fasta", 
                format="fasta")

write.table(metadata_selected %>% select(c('seq_name', 'Collection_Date')), 
            file = "nextclade/input_H1_gisaid_2324_2425/metadata.tsv",
            sep= "\t",
            row.names = FALSE)

#ALIGNMENT performed with Nextclade, roughly half of the sequences sorted out 
#cause too short?! Left with 17687 sequences over 1.5 seasons 

sequence_aln <- readDNAStringSet("nextclade/results_H1_gisaid_2324_2425/nextclade.aligned.fasta")

gisaid_metadata_aln <- gisaid_metadata %>% 
  dplyr::filter(ID %in% names(sequence_aln))


dist_matrix <- pwalign::stringDist(sequence_aln, method = "hamming")
distance_mat = as.matrix(dist_matrix)
time_window = 0.25 #from phylowave paper but for H3N2
mutation_rate = 0.00329 #is it same as clock rate used used for augur refine?
timescale = 0.4 #from phylowave paper but for H3N2
genome_length = 1752

indicies = compute.index.from.distance(distance_mat, 
                                       time_window, 
                                       metadata, 
                                       mutation_rate, 
                                       timescale, 
                                       genome_length)

hist(indicies, breaks = 10, 
     main = "Distribution of Diversity Index", 
     xlab = "Diversity Index", ylab = "Frequency")


#######

indicies = indicies[is.finite(indicies)]
 
hist(indicies)

indicies_df = as.data.frame(indicies)
indicies_df$seqName <- rownames(indicies_df)

nextclade = read.table("nextclade/results_Nexstrain_asOF_Mar26/nextclade.tsv", header = TRUE, sep = "\t")


gr = nextclade[c("index","seqName","subclade")]
# rownames(gr) <- gr$seqName


ggplot(gr, aes(x = subclade)) + 
  geom_bar()

test = merge(gr, indicies_df, by = 'seqName')

test2 = merge(test, metadata, by.x = 'seqName', by.y = 'ID')

test2$week = strftime(test2$time, format = "%V-%Y")

df_summary = test2 %>%
  # Truncate the full date to the start of its week
  dplyr::filter(time >= as.Date("2024-01-01")) %>% 
  mutate(week_date = floor_date(time, unit = "week")) %>% 
  group_by(week_date, subclade) %>% 
  summarise(median_index = median(indicies),
            mean_index = mean(indicies), .groups = "drop")

full_length_thd <- ggplot(df_summary, 
       aes(x = week_date, y = median_index, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE)+
  scale_colour_manual(values = as.vector(pals::cols25(20))) +
  scale_x_date(date_labels = "%V-%Y")



ggplot(df_summary, 
       aes(x = week_date, y = mean_index, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE)+
  scale_colour_manual(values = as.vector(pals::cols25(20))) +
  scale_x_date(date_labels = "%V-%Y")




amp_length_thd

full_length_thd
