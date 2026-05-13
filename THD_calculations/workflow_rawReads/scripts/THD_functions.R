library(tidyr)
library(Biostrings)
library(msa)
library(ggplot2)
library(stringr)
library(tidyverse)
library(pals)

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
  
  #   # Convert dates to numeric years if needed
  #   if (inherits(times, "Date")) {
  #     times <- as.numeric(times) / 365.25  # Convert days since epoch to years
  #   }
  
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

#' Title
#'
#' @param distance_mat 
#' @param time_window 
#' @param metadata 
#' @param mutation_rate 
#' @param timescale 
#' @param genome_length 
#' @param nextclade nextclade DT
#'
#' @returns
#' @export
#'
#' @examples
compute_and_merge_results <- function(distance_mat, time_window, metadata, 
                                      mutation_rate, timescale, genome_length,
                                      nextclade) {
  indicies <- compute.index.from.distance(distance_mat, 
                                          time_window, 
                                          metadata, 
                                          mutation_rate, 
                                          timescale, 
                                          genome_length)
  indicies = indicies[is.finite(indicies)]
  
  indicies_df = as.data.frame(indicies)
  indicies_df$seqName <- rownames(indicies_df)
  
  df_results <- merge(nextclade[,.(index, seqName, subclade, `legacy-clade`)], 
                      indicies_df, by = 'seqName') %>% 
    merge(., metadata, by.x = 'seqName', by.y = 'ID')
  
  return(df_results)
}


