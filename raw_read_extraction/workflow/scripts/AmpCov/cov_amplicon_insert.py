import argparse
import numpy as np
import pandas as pd
import os
import re
import gzip
import subprocess


def parse_args():
    '''Parsing of command line args'''
    parser = argparse.ArgumentParser(
            description="Script to calculate coverage of amplicon-inserts on sequence of interest.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-b", required=True, metavar='PROTOCOL=BED',
                        dest='bed_files', action='append', type=str,
                        help="Mapping of primer_protocol to insert bed file, e.g. v3=path/to/v3.insert.bed. Repeat for each protocol.")
    parser.add_argument("-s", required=True, metavar='TSV',
                        dest='samp_file', help="tsv file like samples.tsv.")
    parser.add_argument("-p", required=True, metavar='PATH',
                        dest='samp_path', help="main path to samples")
    parser.add_argument("-m", required=False, dest='paird_reads', help="Specify reads paird or not.", action='store_true')
    parser.add_argument("-o", required=False, default=os.getcwd(),
                        metavar='PATH', dest='outdir',
                        help="Output directory")

    return parser.parse_args()


def parse_bed_map(bed_args):
    '''Parse list of 'protocol=path' strings into a dict'''
    bed_map = {}
    for entry in bed_args:
        protocol, path = entry.split("=", 1)
        bed_map[protocol] = path
    return bed_map


def get_samples_paths(main_samples_path, samplestsv):
    '''make list of (coverage_path, primer_protocol) tuples from samples.tsv. Credit: https://github.com/dr-david/amplicon_cov'''
    with open(samplestsv, 'r') as f:
        first_field = f.readline().split("\t")[0]
    has_header = first_field == "sample_id"

    df = pd.read_csv(samplestsv, sep="\t", header=0 if has_header else None,
                     names=["sample_id", "batch", "length", "primer_protocol"])

    sam_paths_list = []
    for _, row in df.iterrows():
        path = main_samples_path + "/" + row["sample_id"] + "/" + row["batch"] + "/alignments/coverage.tsv.gz"
        sam_paths_list.append((path, row["primer_protocol"]))

    return sam_paths_list


def read_bed(bed_file):

    df = pd.read_table(bed_file, header=None)
    try:
        df["amplicon"] = df[3].astype(int)
    except (ValueError, TypeError):
        df["amplicon"] = [int(re.search("_([0-9]+)", i).group(1)) for i in df[3]]

    return df

def read_cov(cov_file):

    df = pd.read_csv(cov_file, sep="\t", compression="gzip",names=["Reference","Position","Coverage"], header=0)

    return df
    

def amplicon_median(cov_df,amp_start,amp_end):
    
    amp_median = cov_df[(cov_df.Position >= amp_start) & (cov_df.Position <= amp_end)]["Coverage"].median()

    return amp_median


def make_amplicon_df(bed_file,coverage_file):

    bed_df = read_bed(bed_file)
    cov_df = read_cov(coverage_file)

    s = coverage_file.split("/")[-4]
    
    amplicon_df = pd.DataFrame({"sample": s, "amplicon": bed_df.amplicon, "start_amp":bed_df[1], "end_amp":bed_df[2]})
    amplicon_df['median_coverage']= amplicon_df.apply(lambda row: amplicon_median(cov_df, row['start_amp'], row['end_amp']), axis=1)

    amplicon_df = pd.pivot_table(amplicon_df[['sample','amplicon','median_coverage']], 
               values='median_coverage', index='sample', columns='amplicon')

    amplicon_df.columns.name = None
    amplicon_df = amplicon_df.reset_index()

    
    return amplicon_df

def main():

    args = parse_args()
    samp_file = args.samp_file
    samp_path = args.samp_path
    bed_map = parse_bed_map(args.bed_files)
    seq = args.paird_reads
    outdir = args.outdir

    cov_files = get_samples_paths(samp_path, samp_file)
    amp_dfs = []
    for cov_file, protocol in cov_files:
        if protocol not in bed_map:
            raise ValueError(f"No bed file specified for protocol '{protocol}'. Provide it with -b {protocol}=<path>")
        amp_dfs.append(make_amplicon_df(bed_map[protocol], cov_file))
    temp_df = pd.concat(amp_dfs)
    temp_df = temp_df.set_index('sample')

    if seq == True:
        #for paired read sequencing divide by two for read number per amplicon
        temp_df = temp_df.div(2)
    else:
        pass
    
    amplicon_df = temp_df.reset_index()
    amplicon_norm_df = temp_df.div(temp_df.sum(axis=1), axis=0).reset_index()

    os.makedirs(outdir, exist_ok=True)
    print("\nOutputting .csv's")
    amplicon_df.to_csv(outdir + "/amplicons_coverages.csv", index = False)
    amplicon_norm_df.to_csv(outdir + "/amplicons_coverages_norm.csv", index = False)

if __name__ == '__main__':
    main()

