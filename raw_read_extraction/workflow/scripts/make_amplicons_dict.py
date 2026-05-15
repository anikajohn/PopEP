import pandas as pd
import json
import yaml

amp_bed_file = snakemake.input.bedfile
mut_dict_path = snakemake.input.mut_dict
output_file = snakemake.output[0]

with open(mut_dict_path, "r") as f:
    mut_dict = yaml.safe_load(f)

amp_bed = pd.read_csv(amp_bed_file, sep="\t")

amplicons_dict = {}

for i in range(amp_bed.shape[0]):
    amp_dict = {}
    for pos_str, muts in mut_dict.items():
        pos = int(pos_str)
        if amp_bed.loc[i, "start"] <= pos <= amp_bed.loc[i, "stop"]:
            if not isinstance(muts, list):
                muts = [muts]
            amp_dict[pos] = muts
    amplicons_dict[i + 1] = amp_dict if amp_dict else None

with open(output_file, "w") as f:
    json.dump(amplicons_dict, f, indent=2)
