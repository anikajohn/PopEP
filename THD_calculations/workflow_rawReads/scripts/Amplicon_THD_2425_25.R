library(tidyverse)
library(data.table)
source("THD_functions.R")

####AMP SEQ PREPERATION#######

# 1. amplicon sclicing, followed by nextclade classification

full_length = readDNAStringSet("nextclade/input_H1_gisaid_2425_25/downsampled_sequence.fasta")
metadata = read.table("nextclade/input_H1_gisaid_2425_25/downsampled_metadata.tsv")

amp_all = read.table("nextclade/input_H1_gisaid_2425_25/amplicon_map/H1.amp_all.bed",
                  header = TRUE)


amp2 = subseq(full_length,amp_all[,"start"][1],amp_all[,"end"][1])
amp3 = subseq(full_length,amp_all[,"start"][2],amp_all[,"end"][2])
amp4 = subseq(full_length,amp_all[,"start"][3],amp_all[,"end"][3])
amp5 = subseq(full_length,amp_all[,"start"][4],amp_all[,"end"][4])
amp6 = subseq(full_length,amp_all[,"start"][5],amp_all[,"end"][5])


writeXStringSet(amp2,"nextclade/input_H1_gisaid_2425_25/AMP2/sequence.fasta")
writeXStringSet(amp3,"nextclade/input_H1_gisaid_2425_25/AMP3/sequence.fasta")
writeXStringSet(amp4,"nextclade/input_H1_gisaid_2425_25/AMP4/sequence.fasta")
writeXStringSet(amp5,"nextclade/input_H1_gisaid_2425_25/AMP5/sequence.fasta")
writeXStringSet(amp6,"nextclade/input_H1_gisaid_2425_25/AMP6/sequence.fasta")


##### THD CALCULATION BY AMPLICON ######
amp = "AMP2"
# seq_file = paste0("nextclade/input_H1_gisaid_2425_25/",
                  # amp,
                  # "/sequence.fasta")
# nextclade_file = paste0("nextclade/results_H1_gisaid_2425_25_",
#                         amp,
#                         "/nextclade.tsv")
output_file = paste0("nextclade/results_H1_gisaid_2425_25_",
                     amp,
                     "/index_time_parameters_test.tsv")
seq_file ="nextclade/input_H1_gisaid_2425_25/downsampled_sequence.fasta"
nextclade_file="nextclade/results_H1_gisaid_2425_25/nextclade.tsv"



metadata = read.table("nextclade/input_H1_gisaid_2425_25/downsampled_metadata.tsv",
                      header = TRUE, sep = "\t")
nextclade = read.table(nextclade_file, header = TRUE, sep = "\t")

am_seq = readDNAStringSet(seq_file)
dist_matrix <- pwalign::stringDist(am_seq, method = "hamming")
distance_mat = as.matrix(dist_matrix)
mutation_rate = 0.00329 #is it same as clock rate used used for augur refine?
genome_length = width(am_seq) %>% unique()



# Parameter grid search
timescale_values <- c(0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
timewindow_values <- c(0.019, 0.038, 0.0833, 0.15, 0.23)

# Create all combinations
param_grid <- expand.grid(timescale = timescale_values, 
                          time_window = timewindow_values)

df_all_results <- pmap_dfr(param_grid, function(timescale, time_window) {
  message(sprintf("Running: timescale = %.2f, time_window = %.3f", timescale, time_window))
  
  df_temp <- compute_and_merge_results(
    distance_mat = distance_mat,
    time_window = time_window,
    metadata = metadata,
    mutation_rate = mutation_rate,
    timescale = timescale,
    genome_length = genome_length,
    nextclade = nextclade
  )
  
  df_temp$timescale <- timescale
  df_temp$time_window <- time_window
  
  return(df_temp)
})

# View summary
message(sprintf("Total combinations: %d", nrow(param_grid)))
message(sprintf("Total rows in combined dataframe: %d", nrow(df_all_results)))

write.table(df_all_results,
            output_file,
            sep = "\t")



##### PLOTTING ######

palette_colors <- as.vector(pals::cols25(14))
subclades <- c("C.1.9", "C.1.9.1", "C.1.9.2", "C.1.9.3",
               "C.1.9.4","D","D.1","D.3.1","D.5",
               "C.1","C.1.1","C.1.8","D.2","D.3")
custom_palette <- setNames(palette_colors, subclades)

df_plot = df_all_results %>% 
  mutate(
    date = lubridate::date_decimal(time),
    week_date = floor_date(as.Date(date), unit = "week"))

ggplot(df_plot,
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE)+
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(time_window),
             rows = vars(timescale),
             labeller = label_both)

#TODO get rid of this filtering in commandline script
ggplot(df_plot %>% dplyr::filter(time_window >= 0.0833),
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE)+
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(time_window),
             rows = vars(timescale),
             labeller = label_both)



##### PLOTTING DISTANCE HEATMAPS######

# Get subclade info for each sample and order by subclade
sample_subclade <- nextclade %>%
  dplyr::select(seqName, subclade) %>%
  dplyr::filter(seqName %in% colnames(distance_mat)) %>%
  dplyr::arrange(subclade, seqName)

# Reorder matrix by subclade, with hierarchical clustering within each subclade
ordered_ids <- c()  # Initialize empty vector

for (sc in unique(sample_subclade$subclade)) {
  # Get samples in this subclade
  subclade_samples <- sample_subclade %>%
    dplyr::filter(subclade == sc) %>%
    dplyr::pull(seqName)
  
  if (length(subclade_samples) > 2) {
    # Extract submatrix and perform hierarchical clustering
    sub_dist <- as.dist(distance_mat[subclade_samples, subclade_samples])
    hc <- hclust(sub_dist, method = "average")
    # Get ordered sample names from dendrogram
    ordered_subclade_ids <- subclade_samples[hc$order]
  } else {
    # If only 1-2 samples, keep original order
    ordered_subclade_ids <- subclade_samples
  }
  
  ordered_ids <- c(ordered_ids, ordered_subclade_ids)
}

# # Update sample_subclade to match new ordering (for downstream plotting)
# sample_subclade <- sample_subclade %>%
#   dplyr::slice(match(ordered_ids, seqName))

distance_mat_ordered <- distance_mat[ordered_ids, ordered_ids]

# Convert to long format with subclade info
distance_heatmap_df <- as.data.frame(as.table(distance_mat_ordered)) %>%
  dplyr::rename(sample_1 = Var1,
                sample_2 = Var2,
                distance = Freq) %>%
  dplyr::left_join(sample_subclade, by = c("sample_1" = "seqName")) %>%
  dplyr::rename(subclade_1 = subclade) %>%
  dplyr::left_join(sample_subclade, by = c("sample_2" = "seqName")) %>%
  dplyr::rename(subclade_2 = subclade)

# Preserve ordering as factor levels
distance_heatmap_df$sample_1 <- factor(distance_heatmap_df$sample_1, levels = ordered_ids)
distance_heatmap_df$sample_2 <- factor(distance_heatmap_df$sample_2, levels = ordered_ids)

subclade_colors <- custom_palette

# Calculate subclade positions for axis labels (midpoint of each group)
subclade_positions <- sample_subclade %>%
  mutate(position = match(seqName, ordered_ids)) %>%
  group_by(subclade) %>%
  summarise(mid_pos = mean(position), 
            start_pos = min(position) - 0.5,
            end_pos = max(position) + 0.5,
            .groups = "drop")

# Get boundary positions between subclades (for separator lines)
subclade_boundaries <- subclade_positions$end_pos[-nrow(subclade_positions)]

# Create annotation bars data for x and y axes
n_samples <- length(ordered_ids)
annotation_width <- n_samples * 0.03  # Width of annotation bar

x_annotation <- sample_subclade %>%
  mutate(position = match(seqName, ordered_ids),
         ymin = -annotation_width,
         ymax = 0)

y_annotation <- sample_subclade %>%
  mutate(position = match(seqName, ordered_ids),
         xmin = -annotation_width,
         xmax = 0)

ggplot(distance_heatmap_df, aes(x = as.numeric(sample_1), y = as.numeric(sample_2))) +
  # Main heatmap
  geom_tile(aes(fill = distance)) +
  # X-axis annotation bar (bottom)
  geom_rect(data = x_annotation,
            aes(xmin = position - 0.5, xmax = position + 0.5,
                ymin = ymin, ymax = ymax),
            fill = subclade_colors[x_annotation$subclade],
            inherit.aes = FALSE) +
  # Y-axis annotation bar (left)
  geom_rect(data = y_annotation,
            aes(xmin = xmin, xmax = xmax,
                ymin = position - 0.5, ymax = position + 0.5),
            fill = subclade_colors[y_annotation$subclade],
            inherit.aes = FALSE) +
  # Subclade separator lines (vertical)
  geom_vline(xintercept = subclade_boundaries, color = "black", linewidth = 0.3) +
  # Subclade separator lines (horizontal)
  geom_hline(yintercept = subclade_boundaries, color = "black", linewidth = 0.3) +
  coord_equal(clip = "off") +
  scale_fill_gradient(low = "white", high = "firebrick", na.value = "grey90") +
  scale_x_continuous(breaks = subclade_positions$mid_pos,
                     labels = subclade_positions$subclade,
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = subclade_positions$mid_pos,
                     labels = subclade_positions$subclade,
                     expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = "Hamming\ndistance",
       title = "Pairwise distance matrix (ordered by subclade)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 10, 10, 30)) 

