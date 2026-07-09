#!/usr/bin/env python3
"""
Command-line tool to convert basecnt files to mutation TSV format.

Usage example:
    python basecnt_to_mutation_tsv.py \
        --fasta /path/to/H3.fasta \
        --reference CY163680.1 \
        --timeline /path/to/timeline.tsv \
        --basecnt '/path/to/basecnt/*/*/alignments/basecnt.tsv.gz' \
        --output /path/to/output/mutations.tsv \
        --location "Zürich (ZH)" \
        --location_code 10 \
        --time_start 2025-10-16 \
        --time_end 2026-07-11
"""

import argparse
import numpy as np
import pandas as pd
import json
import glob
import re
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description='Convert basecnt files to mutation TSV format for Influenza analysis.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--fasta', '-f', required=True,
                        help='Path to the FASTA file containing reference sequences')
    parser.add_argument('--reference', '-r', required=True,
                        help='Reference sequence ID in the FASTA file (e.g., CY163680.1)')
    parser.add_argument('--timeline', '-t', required=True,
                        help='Path to the timeline TSV file')
    parser.add_argument('--basecnt', '-b', required=True,
                        help='Glob pattern for basecnt.tsv.gz files (use quotes)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file for the results')
    parser.add_argument('--location', '-l', default='Zürich (ZH)',
                        help='Location filter for samples')
    parser.add_argument('--location_code', '-c', default='10',
                        help='Location code for output filename')
    parser.add_argument('--time_start', default='2025-10-16',
                        help='Start date for filtering samples (YYYY-MM-DD)')
    parser.add_argument('--time_end', default='2026-07-11',
                        help='End date for filtering samples (YYYY-MM-DD)')
    parser.add_argument('--threshold', default=0.005, type=float,
                        help='Frequency threshold value for mutations')
    parser.add_argument('--rare_mutation_days', default=2, type=int,
                        help='Minimum number of days a mutation must appear above threshold')
    parser.add_argument('--coverage_threshold', default=10, type=int,
                        help='Minimum coverage threshold')

    return parser.parse_args()


def read_fasta(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def extract_sample_ID(timeline_file_dir, time_start, time_end, location):
    """Extract sample IDs from selected time period, location, and protocol."""
    timeline_file = pd.read_csv(timeline_file_dir, sep='\t', usecols=['sample', 'proto', 'date', 'location'],
                                encoding='utf-8')
    
    # some batch set ups had bugs in the name
    timeline_file = timeline_file[~timeline_file['proto'].astype(str).str.contains('<built-in', na=False)]

    timeline_file["date"] = pd.to_datetime(timeline_file["date"], errors='coerce', format='mixed')


    selected_rows = timeline_file[
        (timeline_file["date"] > time_start) &
        (timeline_file["date"] < time_end) &
        (timeline_file["location"].isin([location]))
    ]
    samples_ID = selected_rows[['sample', 'date']]
    return samples_ID


def extract_cov(basecnt_tsv_file, sample_name, coverage_threshold):
    """Extract coverage for each nucleotide at each genome position and make a frequency np array."""
    sample_batch = pd.read_csv(
        filepath_or_buffer=basecnt_tsv_file,
        sep='\t',
        header=0,
        usecols=[2]).columns[0]
    
    #different naming conventions in bascnt headers
    if '/' in sample_batch:
        sample = sample_batch.split('/')[0]
    else:
        if sample_name in sample_batch:
            sample = sample_name
        else:
            sample = sample_batch


    if sample_name != sample:
        print(f"Sample name mismatch! Expected '{sample_name}', found '{sample}'.")

    test_basecnt = pd.read_csv(
        filepath_or_buffer=basecnt_tsv_file,
        sep='\t',
        header=2,
        usecols=[0, 1, 2, 3, 4, 5])
    test_basecnt.columns = ['ref', 'pos', 'A', 'C', 'G', 'T']

    test_basecnt['total'] = test_basecnt[['A', 'C', 'G', 'T']].sum(axis=1)
    test_basecnt.loc[test_basecnt['total'] < coverage_threshold, ['A', 'C', 'G', 'T', 'total']] = np.nan

    freq_basecnt = test_basecnt[['A', 'C', 'G', 'T']].div(test_basecnt['total'], axis=0)
    freq_basecnt['pos'] = test_basecnt['pos']
    freq_basecnt['date'] = sample_name
    freq_basecnt['location'] = sample_name.split('_')[1] if len(sample_name.split('_')) > 1 else 'unknown'
    freq_basecnt['total_cov'] = test_basecnt['total']

    return freq_basecnt


def process_multiple_basecnt_files(input_dir, sample_IDs, coverage_threshold):
    """Process multiple basecnt files and return combined dataframe."""
    basecnt_files = glob.glob(input_dir, recursive=True)
    rows = []

    for basecnt_file in basecnt_files:
        sample_name = basecnt_file.split('/')[-4]

        if sample_name in sample_IDs.iloc[:, 0].values:
            df = extract_cov(basecnt_file, sample_name, coverage_threshold)
            rows.append(df)

    basecnt_out = pd.concat(rows, axis=0, ignore_index=False)
    return basecnt_out


def create_mut_freq_dict(col):
    """Create a dictionary of mutation proportions for a sample."""
    mutation_proportions = {}
    for mutation, frequency in col.items():
        if pd.isna(frequency):
            mutation_proportions[mutation] = None
        else:
            mutation_proportions[mutation] = frequency

    mutation_json = json.dumps(mutation_proportions, ensure_ascii=False)
    return mutation_json


def main():
    args = parse_args()

    ref_seq = read_fasta(args.fasta)[args.reference]
    sample_IDs = extract_sample_ID(args.timeline, args.time_start, args.time_end, args.location)

    collected_basecnt_freq = process_multiple_basecnt_files(
        args.basecnt, sample_IDs, args.coverage_threshold
    )[['date', 'location', 'pos', 'A', 'C', 'G', 'T']]

    collected_coverage_A = collected_basecnt_freq.pivot_table(values='A', index='pos', columns='date', sort=False)
    collected_coverage_A.index = [f'{ref_seq[int(pos)-1]}{int(pos)}A' for pos in collected_coverage_A.index]
    collected_coverage_A = collected_coverage_A[~collected_coverage_A.index.str.startswith('A')]

    collected_coverage_C = collected_basecnt_freq.pivot_table(values='C', index='pos', columns='date', sort=False)
    collected_coverage_C.index = [f'{ref_seq[int(pos)-1]}{int(pos)}C' for pos in collected_coverage_C.index]
    collected_coverage_C = collected_coverage_C[~collected_coverage_C.index.str.startswith('C')]

    collected_coverage_G = collected_basecnt_freq.pivot_table(values='G', index='pos', columns='date', sort=False)
    collected_coverage_G.index = [f'{ref_seq[int(pos)-1]}{int(pos)}G' for pos in collected_coverage_G.index]
    collected_coverage_G = collected_coverage_G[~collected_coverage_G.index.str.startswith('G')]

    collected_coverage_T = collected_basecnt_freq.pivot_table(values='T', index='pos', columns='date', sort=False)
    collected_coverage_T.index = [f'{ref_seq[int(pos)-1]}{int(pos)}T' for pos in collected_coverage_T.index]
    collected_coverage_T = collected_coverage_T[~collected_coverage_T.index.str.startswith('T')]

    collected_frequency = pd.concat(
        [collected_coverage_A, collected_coverage_C, collected_coverage_G, collected_coverage_T],
        axis=0, ignore_index=False
    )

    num_greater_than_threshold = (collected_frequency >= args.threshold).sum(axis=1)
    collected_frequency = collected_frequency[num_greater_than_threshold >= args.rare_mutation_days]

    nucleotide_mut_freq = collected_frequency.apply(create_mut_freq_dict, axis=0)

    timeline_file = pd.read_csv(args.timeline, sep='\t', usecols=['sample', 'proto', 'date', 'location', 'batch'],
                                encoding='utf-8')
    
    
    #have 09-31 somewhere
    timeline_file["date"] = pd.to_datetime(timeline_file["date"], errors='coerce')

    # Filter out the rows with the broken Python regex object string (Sars-CoV-2) in the 'batch' column
    timeline_file = timeline_file[~timeline_file['batch'].str.contains('method|<', na=False)]

    #work-around cause origianlly batch column was not needed and (back to back compatibality with prep samples)
    timeline_file = timeline_file.drop(columns=['batch'])


    selected_rows = timeline_file[
        (timeline_file["date"] > args.time_start) &
        (timeline_file["date"] < args.time_end) &
        (timeline_file["location"].isin([args.location]))
    ]
    tsv_samples_locations = selected_rows.copy()

    valid_samples = set(tsv_samples_locations['sample'].values)
    filtered_nucleotide_mut_freq = {
        sample: freq for sample, freq in nucleotide_mut_freq.items() if sample in valid_samples
    }
    tsv_samples_locations['nucleotideMutationFrequency'] = tsv_samples_locations["sample"].map(
        filtered_nucleotide_mut_freq)

    tsv_samples_locations.to_csv(args.output, sep='\t', index=False, quoting=3)
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
