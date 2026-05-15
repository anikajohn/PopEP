import pandas as pd
import pysam
import json
import yaml
import argparse
import os



def extract_sample_ID(timeline_file_dir, datefrom, dateto, location):

    timeline_file = pd.read_csv(timeline_file_dir, sep = '\t', usecols = ['sample', 'batch', 'proto', 'date', 'location'],
                                encoding='utf-8')
    # convert the "date" column to datetime type:
    timeline_file["date"] = pd.to_datetime(timeline_file["date"])

    selected_rows = timeline_file[
        # select the rows with date from date 1 to date 2 (according to .tsv file)
        (timeline_file["date"] > datefrom) &
        (timeline_file["date"] < dateto) &
        (timeline_file["location"].isin([location]))
        # filtering condition to take only Artic v4.1 protocol:
       # (timeline_file["proto"] == "v41")
    ]
    samples_ID = selected_rows[['sample', 'batch']]
    samples_ID_batch = pd.DataFrame()
    samples_ID_batch['sample_batch'] = samples_ID['sample'] + '/' + samples_ID['batch']

    return samples_ID_batch


def main(input_timeline_tsv, datefrom, dateto, location_list, output_yaml):

    final_samples = list()
    for location in location_list:

        samples = extract_sample_ID(input_timeline_tsv, datefrom, dateto, location)
        sample_list = samples['sample_batch'].tolist()
        final_samples.extend(sample_list)

#    samples = extract_sample_ID(input_timeline_tsv, datefrom, dateto, location)
#    sample_list = samples['sample_batch'].tolist()

    with open(output_yaml, 'w') as file:
        yaml.dump({'sample_batch': final_samples}, file, default_flow_style=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='make config/samples.yaml file')
    parser.add_argument("input_timeline_tsv")
    parser.add_argument("datefrom")
    parser.add_argument("dateto")
    parser.add_argument("location")
    parser.add_argument("output_yaml")

    args = parser.parse_args()

    output_path = os.path.expanduser(args.output_yaml)  
    main(args.input_timeline_tsv, args.datefrom, args.dateto, [args.location], output_path)
    print(f"YAML file created at: {output_path}")
