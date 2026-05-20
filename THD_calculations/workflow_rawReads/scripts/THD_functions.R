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
#' 
compute.index.from.distance = function(distance_mat, time_window, metadata, mutation_rate, timescale, genome_length,
                                        chunk_size = 500L, mc.cores = 1L){
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

  n           <- length(sample_ids)
  scale_fac   <- mutation_rate * genome_length / 2
  bandwidth   <- index.bandwidth(timescale, genome_length, mutation_rate)

  Bsums         <- numeric(n)
  window_counts <- integer(n)

  ## Process columns in chunks to cap peak memory at n * chunk_size entries.
  ## Each chunk builds only an (n x chunk_size) sub-matrix instead of the full
  ## (n x n) matrix, replacing the O(n^2) element-wise R loop with vectorised
  ## C-level arithmetic.
  chunk_breaks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))

  chunk_results <- parallel::mclapply(chunk_breaks, function(cidx) {
    ## Pairwise absolute time differences: (n x |cidx|)
    td <- abs(outer(times, times[cidx], "-"))

    ## Within-window mask; diagonal (self-comparisons) is excluded
    is_diag <- outer(seq_len(n), cidx, `==`)
    in_win <- (td <= time_window) & !is_diag

    ## Distance sub-matrix with temporal correction applied only within window,
    ## then scaled; out-of-window entries are set to NA so colSums ignores them
    d_scaled <- (distance_mat[, cidx, drop = FALSE] + td) * scale_fac
    d_scaled[!in_win] <- NA_real_

    B <- bandwidth ^ d_scaled

    list(
      cidx          = cidx,
      Bsums_chunk   = colSums(B, na.rm = TRUE),
      window_chunk  = colSums(in_win)
    )
  }, mc.cores = mc.cores)

  for (res in chunk_results) {
    Bsums[res$cidx]         <- res$Bsums_chunk
    window_counts[res$cidx] <- res$window_chunk
  }

  index <- Bsums / window_counts
  names(index) <- sample_ids

  return(index)
}




#' Compute the diversity index from a haplotype distance matrix and haplotype counts
#'
#' Computes the diversity index directly from a pairwise distance matrix between
#' samples, without requiring a phylogenetic tree. This is useful when users have
#' precomputed genetic distances (e.g., from sequence alignments) without a tree.
#'
#' @param distance_mat Pairwise distance matrix (e.g., Hamming distances) between samples
#' @param count_file Tab-delimited file with columns \code{haplo_name} and \code{count}
#' @param time_window Time window for the index computation (positive number, in years)
#' @param metadata Data frame with columns \code{ID} and \code{time}
#' @param mutation_rate Mutation rate (in mutation/site/year)
#' @param timescale Index timescale (in years)
#' @param genome_length Genome length (in bp)
#' @return Named vector of index values for each sample
#' @export
#' 
compute.index.from.haplo.distance = function(distance_mat, count_file,
time_window, metadata, mutation_rate, timescale, 
genome_length,chunk_size = 500L, mc.cores = 1L){

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

  n           <- length(sample_ids)
  scale_fac   <- mutation_rate * genome_length / 2
  bandwidth   <- index.bandwidth(timescale, genome_length, mutation_rate)

  Bsums         <- numeric(n)
  window_counts <- integer(n)

  df_count = read.table(count_file, sep="\t", header=TRUE)
  count <- df_count$count
  names(count) <- df_count$read_name
  count <- count[sample_ids] # align to matrix row/col order by name

  baseline = (count-1)

  ## Process columns in chunks to cap peak memory at n * chunk_size entries.
  ## Each chunk builds only an (n x chunk_size) sub-matrix instead of the full
  ## (n x n) matrix, replacing the O(n^2) element-wise R loop with vectorised
  ## C-level arithmetic.
  chunk_breaks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))

  chunk_results <- parallel::mclapply(chunk_breaks, function(cidx) {
    ## Pairwise absolute time differences: (n x |cidx|)
    td <- abs(outer(times, times[cidx], "-"))

    ## Within-window mask; diagonal (self-comparisons) is excluded
    is_diag <- outer(seq_len(n), cidx, `==`)
    in_win <- (td <= time_window) & !is_diag

    ## Distance sub-matrix with temporal correction applied only within window,
    ## then scaled; out-of-window entries are set to NA so colSums ignores them
    d_scaled <- (distance_mat[, cidx, drop = FALSE] + td) * scale_fac
    d_scaled[!in_win] <- NA_real_

    B <- bandwidth ^ d_scaled

    #accounts for the fact that if if have the same sequence n times will have 
    #n-1 time B =1 and the window counts will have n-1 additional counts

    list(
      cidx          = cidx,
      Bsums_chunk   = baseline[cidx] + colSums(B * count, na.rm = TRUE),
      window_chunk  = baseline[cidx] + colSums(in_win * count)
    )
  }, mc.cores = mc.cores)

  for (res in chunk_results) {
    Bsums[res$cidx]         <- res$Bsums_chunk
    window_counts[res$cidx] <- res$window_chunk
  }

  index <- Bsums / window_counts
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
                                      nextclade, mc.cores = 1L) {
  
  indicies <- compute.index.from.distance(distance_mat, 
                                          time_window, 
                                          metadata, 
                                          mutation_rate, 
                                          timescale, 
                                          genome_length,
                                          mc.cores = mc.cores)
  indicies = indicies[is.finite(indicies)]
  
  indicies_df = as.data.frame(indicies)
  indicies_df$seqName <- rownames(indicies_df)
  
  df_results <- merge(nextclade[,.(index, seqName, subclade, `legacy-clade`)], 
                      indicies_df, by = 'seqName') %>% 
    merge(., metadata, by.x = 'seqName', by.y = 'ID')
  
  return(df_results)
}


