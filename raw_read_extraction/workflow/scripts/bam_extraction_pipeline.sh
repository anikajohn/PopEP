#!/bin/bash
#SBATCH --job-name="bam_extraction"
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=64
####### --partition=light
#SBATCH --time=72:00:00

#TODO: intermediate outputs should be skipped or compressed

# 1. filter out amplicon 3 reads 
find ../results_Aire/ -type f -name "REF_aln_trim.bam" | while read bam_file; do
    dir=$(dirname "$bam_file")
    barcode=$(echo "$bam_file" | grep -o 'barcode[0-9]\+')
    tmpfile=$(mktemp)
    # base=$(basename "$dir")
    
    # 1. filter out amplicon 3 reads
    python filtering_bam.py -i "$bam_file" -o "$dir/REF_aln_amp3.bam" -s 473 -e 917 

    samtools sort "$dir/REF_aln_amp3.bam" -o "$dir/REF_aln_amp3_sorted.bam" 

    samtools index "$dir/REF_aln_amp3_sorted.bam" 

    rm "$dir/REF_aln_amp3.bam"

    # 2. Convert to fasta
    python BamToFasta.py -i "$dir/REF_aln_amp3_sorted.bam" -o "$dir/WWAmplicon.fasta" 

    # 4. Removing identical sequences #PROBLEM FOR FREQUENCY TREES/ DIVERSITY QUANTIFICATION
    # seqkit rmdup -s < "$dir/WWAmplicon.fasta" > "$dir/WWAmplicon_clean.fasta"
    seqkit rmdup -s < "$dir/WWAmplicon.fasta" | gzip > "$dir/WWAmplicon_clean.fasta.gz"

    gzip "$dir/WWAmplicon.fasta" 

    #sample identifier to header to make metadata later
    zcat "$dir/WWAmplicon.fasta.gz" | awk -v bc="$barcode" '/^>/ {print ">" bc "_" substr($0, 2)} !/^>/ {print}' | gzip > "$tmpfile"
    mv "$tmpfile" "$dir/WWAmplicon.fasta.gz"

    # 5. Sub-sample fasta randomly (e.g., 10% of sequences)
    # seqkit sample -p 0.1 "$dir/WWAmplicon_clean.fasta" -o "$dir/WWAmplicon_subsampled.fasta"
done


# Combine like all clinical data is combined over the season, sub-sampled for easier processing
find ../results_Aire/ -type f -name "WWAmplicon.fasta.gz" -exec zcat {} + | seqkit sample -p 0.1 -o WWAmplicon_subsampled.fasta


