#!/usr/bin/env python3
"""
Script to add meta info to nextclade classification.
"""

import pandas as pd
import argparse

def add_metadata_to_nextclade(metadata_file, nextclade_file, location, amplicon, output_file):

    metadata = pd.read_csv(metadata_file, sep="\t")
    nextclade = pd.read_csv(nextclade_file, sep="\t")
    
    new_nextclade = nextclade.merge(metadata, left_on="seqName", right_on="strain")
    new_nextclade['amplicon'] = amplicon
    new_nextclade['location'] = location
    
    new_nextclade.to_csv(output_file, sep="\t", index=False)

    return new_nextclade

def main():
    parser = argparse.ArgumentParser(description="Add metadata to nextclade results.")
    parser.add_argument("-m", "--metadata_file", help="Path to metadata file.")
    parser.add_argument("-n", "--nextclade_file", help="Path to nextclade results file.")
    parser.add_argument("-l", "--location", help="Location of the samples.")
    parser.add_argument("-a", "--amplicon", help="Amplicon name.")
    parser.add_argument("-o", "--output_file", help="Path to output file.")

    args = parser.parse_args()

    add_metadata_to_nextclade(
        metadata_file=args.metadata_file,
        nextclade_file=args.nextclade_file,
        location=args.location,
        amplicon=args.amplicon,
        output_file=args.output_file
    )

if __name__ == "__main__":
    main()
