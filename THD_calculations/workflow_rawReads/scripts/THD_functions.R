library(tidyr)
library(Biostrings)
# library(msa)
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
  
  #2. Nexclade filters out low quality sequences, need to account for this

  matched  <- sample_ids %in% metadata$ID
  n_lost   <- sum(!matched)
  pct_lost <- 100 * n_lost / length(sample_ids)

  if (n_lost > 0) {
    message(sprintf("%d of %d sample IDs (%.1f%%) not found in metadata$ID.",
                     n_lost, length(sample_ids), pct_lost))
  }

  # #50% is a bit abitrary here
  # if (pct_lost > max_loss_pct) {
  #   stop(sprintf("More than %d%% of sample IDs (%.1f%%) missing from metadata$ID.",
  #                 50, pct_lost))
  # }

  metadata <- metadata[match(sample_ids[matched], metadata$ID), ]

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
#' @param dist_file Pairwise distance matrix (e.g., Hamming distances) between samples
#' @param count_file Tab-delimited file with columns \code{haplo_name} and \code{count}
#' @param time_window Time window for the index computation (positive number, in years)
#' @param metadata Data frame with columns \code{ID} and \code{time}
#' @param mutation_rate Mutation rate (in mutation/site/year)
#' @param timescale Index timescale (in years)
#' @param genome_length Genome length (in bp)
#' @return Named vector of index values for each sample
#' @export
#' 
compute.index.from.haplo.distance = function(dist_file, count_file,
time_window, metadata, mutation_rate, timescale, 
genome_length,chunk_size = 500L, mc.cores = 1L){

  ## Input validation
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

  

  # --- Chunked distance-matrix computation via fread + sequential connection ---
  
  # 1. Extract sample IDs cheaply: shell `cut` reads only bytes up to the first
  #    tab per line, so it touches a tiny fraction of a large file.
  sample_ids <- fread(
    cmd    = if (grepl("\\.gz$", dist_file))
               sprintf("gunzip -c '%s' | cut -f1", dist_file)
             else
               sprintf("cut -f1 '%s'", dist_file),
    header = FALSE,
    sep    = "\n"
  )[[1]]
  
  #some trailing empty lines in distance mat probs.
  sample_ids <- sample_ids[nzchar(sample_ids)]
  
  n <- length(sample_ids)

  

  #2. Nexclade filters out low quality sequences, need to account for this

  matched  <- sample_ids %in% metadata$ID
  n_lost   <- sum(!matched)
  pct_lost <- 100 * n_lost / n

  if (n_lost > 0) {
    message(sprintf("%d of %d sample IDs (%.1f%%) not found in metadata$ID:",
                     n_lost, n, pct_lost))
    message(paste(sample_ids[!matched], collapse = "\n"))
  }

  # #50% is a bit abitrary here
  # if (pct_lost > 50) {
  #   stop(sprintf("More than %d%% of sample IDs (%.1f%%) missing from metadata$ID.",
  #                 50, pct_lost))
  # }

  metadata <- metadata[match(sample_ids[matched], metadata$ID), ]
  
  
  times <- metadata$time[match(sample_ids, metadata$ID)]
  
 
  
  
  message(paste("not matching IDs:", metadata[is.na(metadata$time), ], collapse = "\n"))

  
  if (any(is.na(times))) {
    stop("Some sample IDs in 'distance_mat' are not found in 'metadata$ID'.")
  }

  
  df_count  <- read.table(count_file, sep = "\t", header = TRUE)
  df_count  <- df_count[df_count$read_name != "read_name", ]
  df_count$count <- as.numeric(df_count$count)
  count_vec <- setNames(df_count$count, df_count$read_name)[sample_ids]
  baseline  <- count_vec - 1L
  scale_fac <- mutation_rate * genome_length / 2
  bandwidth <- index.bandwidth(timescale, genome_length, mutation_rate)
  
  # 3. Accumulators (length-N vectors — the only large objects kept in memory)
  Bsums         <- numeric(n)
  window_counts <- numeric(n)
  
  # 4. Open a sequential connection and process chunk_size rows at a time.
  #    readLines(con, n=chunk_size) advances the file pointer without re-scanning,
  #    so the file is read exactly once from start to finish.
  con        <- if (grepl("\\.gz$", dist_file)) gzfile(dist_file, open = "r") else file(dist_file, open = "r")
  row_offset <- 0L
  
  repeat {
    lines <- readLines(con, n = chunk_size)
    if (length(lines) == 0L) break
    
    # fread parses the raw text lines at C speed; first col is seq name
    chunk   <- fread(paste(lines, collapse = "\n"), header = FALSE, sep = "\t")
    chunk_n <- nrow(chunk)
    cidx    <- (row_offset + 1L):(row_offset + chunk_n)
    row_offset <- row_offset + chunk_n
    
    # Transpose: read chunk_n rows → n x chunk_n column matrix.
    # Valid because the distance matrix is symmetric: row j == col j.
    # After transpose the chunk has the same shape as the original column-chunk
    # (n x chunk_size), so all downstream logic is unchanged.
    dist_chunk <- t(as.matrix(chunk[, -1L]))
    storage.mode(dist_chunk) <- "double"
    
    # Identical to original column-chunk logic with cidx = current row indices
    td      <- abs(outer(times, times[cidx], "-"))
    is_diag <- outer(seq_len(n), cidx, `==`)
    in_win  <- (td <= time_window) & !is_diag
    
    d_scaled <- (dist_chunk + td) * scale_fac
    d_scaled[!in_win] <- NA_real_
    B <- bandwidth ^ d_scaled
    
    Bsums[cidx]         <- baseline[cidx] + colSums(B * count_vec, na.rm = TRUE)
    window_counts[cidx] <- baseline[cidx] + colSums(in_win * count_vec)
    print(sprintf("Processed cols %d-%d of %d", cidx[1], cidx[chunk_n], n))
    gc()
    
  }
  close(con)
  
  # 5. Compute the final index (baseline already applied per-chunk above)
  index_chunked <- Bsums / window_counts
  names(index_chunked) <- sample_ids
  

  return(index_chunked)
}




#'
#' @returns
#' @export
#'
#' @examples
compute.merge.haplo.results <- function(dist_file,  time_window, metadata, 
                                      mutation_rate, timescale, genome_length,
                                      nextclade, count_file, mc.cores = 1L, 
                                      chunk_size = 500L) {
  
  indicies <- compute.index.from.haplo.distance(dist_file, 
                                          count_file,
                                          time_window, 
                                          metadata, 
                                          mutation_rate, 
                                          timescale, 
                                          genome_length,
                                          chunk_size,
                                          mc.cores = mc.cores)
  indicies = indicies[is.finite(indicies)]
  
  indicies_df <- data.frame(
  seqName = names(indicies),
  indicies   = as.numeric(indicies))

  message(sprintf("Total number rows before merge: %d", nrow(indicies_df)))

  
  df_results <- merge(nextclade[,.(seqName, subclade, `legacy-clade`)], 
                      indicies_df, by = 'seqName') %>% 
    merge(., metadata, by.x = 'seqName', by.y = 'ID')
  
  message(sprintf("Total number rows after merge: %d", nrow(df_results)))
  
  return(df_results)
  # return(indicies_df)
  # return(indicies_df)

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
compute_and_merge_results <- function(distance_mat, count_file, time_window, metadata, 
                                      mutation_rate, timescale, genome_length,
                                      nextclade, mc.cores = 1L) {
  
  indicies <- compute.index.from.haplo.distance(distance_mat,
                                                count_file,
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


