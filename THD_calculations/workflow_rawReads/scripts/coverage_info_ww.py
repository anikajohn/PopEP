import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Add location and amplicon columns to coverage file")
parser.add_argument("-i", "--input", required=True, help="Input coverage file")
parser.add_argument("-o", "--output", required=True, help="Output TSV file")
parser.add_argument("-l", "--location", required=True, help="Location name")
parser.add_argument("-a", "--amplicon", required=True, help="Amplicon name")

args = parser.parse_args()

colnames = ["sample_name", "amplicon_name", "ref_name", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"] 
df = pd.read_csv(args.input, sep="\t", names=colnames, header=0)
df["location"] = args.location
df["amplicon"] = args.amplicon
df.to_csv(args.output, sep="\t", index=False)
