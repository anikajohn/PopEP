import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Combine files across locations")
parser.add_argument("--fasta-files", nargs="+", required=True, help="Input FASTA files")
parser.add_argument("--metadata-files", nargs="+", required=True, help="Input metadata TSV files")
parser.add_argument("--nextclade-meta-files", nargs="+", required=True, help="Input nextclade metadata TSV files")
parser.add_argument("--coverage-files", nargs="+", required=True, help="Input coverage TSV files")
parser.add_argument("--out-fasta", required=True, help="Output combined FASTA file")
parser.add_argument("--out-metadata", required=True, help="Output combined metadata TSV")
parser.add_argument("--out-nextclade-meta", required=True, help="Output combined nextclade metadata TSV")
parser.add_argument("--out-coverage", required=True, help="Output combined coverage TSV")

args = parser.parse_args()

# Combine FASTA files
with open(args.out_fasta, 'w') as outf:
    for f in args.fasta_files:
        with open(f, 'r') as inf:
            outf.write(inf.read())

# Combine TSV files (with headers)
def combine_tsvs(file_list, output_path):
    dfs = [pd.read_csv(f, sep="\t") for f in file_list]
    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(output_path, sep="\t", index=False)

combine_tsvs(args.metadata_files, args.out_metadata)
combine_tsvs(args.nextclade_meta_files, args.out_nextclade_meta)
combine_tsvs(args.coverage_files, args.out_coverage)
