import yaml
import csv
import ast
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

csv.field_size_limit(100000000)


def read_amp_dict_yaml(file_path):
    """Read the amp_dict file as a dictionary using YAML/JSON parser."""
    with open(file_path, 'r') as file:
        amp_dict = yaml.safe_load(file)
    return amp_dict


def simulate_mutated_reads(csv_file_path, fasta_file, amp_file, bed_file, sample_name):
    """
    Process cooccurrences CSV and reconstruct unique read sequences with their raw counts.

    Returns:
        list of (amplicon_sequence, raw_count) tuples — one entry per unique co-occurrence pattern
    """
    from load_bed import load_bed

    amp_dict_full = read_amp_dict_yaml(amp_file)
    # amp_dict is keyed by amplicon number (1-based); take the first key since
    # this script is invoked per-amplicon and the BAM is already scoped to one amplicon
    first_key = str(list(amp_dict_full.keys())[0])
    amp_dict = amp_dict_full[first_key]

    # amp_dict keys may be int (from JSON) or str (from YAML) — normalise to str
    amp_dict = {str(k): v for k, v in amp_dict.items()}

    bed = load_bed(bed_file)
    records = list(SeqIO.parse(fasta_file, "fasta"))
    pos_mutated = list(amp_dict.keys())

    start = int(bed['start'].iloc[0] - 1)  # 0-based
    stop = int(bed['stop'].iloc[0])

    possible_reads = []

    with open(csv_file_path, 'r', encoding='utf-8') as file:
        csv_reader = csv.reader(file)

        for row in csv_reader:
            overall_dict = ast.literal_eval(row[1])

            for tuple_key, value_dict in overall_dict.items():
                if not (isinstance(value_dict, dict) and
                        all(k in value_dict for k in ['amplicon', 'covered', 'counts', 'mutated'])):
                    print(f"Error parsing row {row}")
                    continue

                read_seq = str(records[0].seq)
                read_list = list(read_seq)
                covered = value_dict['covered']
                count = value_dict['counts']   # raw integer count from scan_bam.py
                mutated = value_dict['mutated']

                for i in range(len(mutated)):
                    pos_0 = int(pos_mutated[i]) - 1  # 0-based

                    if covered[i] == 1 and mutated[i] == 1:
                        read_list[pos_0] = amp_dict[pos_mutated[i]][0]

                    elif covered[i] == 1 and mutated[i] == 0:
                        pass  # keep reference base

                    elif covered[i] == 0 and mutated[i] == 0:
                        read_list[pos_0] = "N"

                    elif covered[i] == 0 and mutated[i] == 1:
                        raise ValueError(
                            f"covered=0 but mutated=1 at position {pos_mutated[i]}.")

                    elif covered[i] == 1 and isinstance(mutated[i], list):
                        if 1 in mutated[i]:
                            j = mutated[i].index(1)
                            read_list[pos_0] = amp_dict[pos_mutated[i]][j]
                        # else: not mutated, keep reference base

                    elif covered[i] == 0 and isinstance(mutated[i], list):
                        if 1 in mutated[i]:
                            j = mutated[i].index(1)
                            raise ValueError(
                                f"covered=0 but mutated=1 at position {pos_mutated[i]} mut index {j}.")
                        read_list[pos_0] = "N"

                    else:
                        raise ValueError(
                            f"Unexpected case at position {pos_mutated[i]}: "
                            f"covered={covered[i]}, mutated={mutated[i]}.")

                amplicon_seq = ''.join(read_list)[start:stop]
                possible_reads.append((amplicon_seq, count))

    return possible_reads


def write_fasta_and_counts(possible_reads, sample_name, amplicon, location, output_fasta, output_tsv):
    """
    Write each unique read ONCE to a FASTA and record its raw count in a TSV.

    Output FASTA header:  {sample_name}_readVersion_{i}
    Output TSV columns:   read_name, count, amplicon, location
    """
    seq_records = []
    tsv_rows = []

    for i, (seq_str, count) in enumerate(possible_reads):
        header = f"{sample_name}_readVersion_{i}"
        seq_records.append(SeqRecord(Seq(seq_str), id=header, description=""))
        tsv_rows.append((header, count))

    SeqIO.write(seq_records, output_fasta, "fasta")
    print(f"FASTA written: {output_fasta} ({len(seq_records)} unique reads)")

    with open(output_tsv, 'w') as f:
        f.write("read_name\tcount\tamplicon\tlocation\n")
        for header, count in tsv_rows:
            f.write(f"{header}\t{count}\t{amplicon}\t{location}\n")
    print(f"Counts TSV written: {output_tsv}")


def main(csv_file_path, fasta_file, amp_file, bed_file, sample_name, amplicon, location, output_fasta, output_tsv):

    print("Step 0: Checking if cooccurrences file is empty...")
    with open(csv_file_path, 'r', encoding='utf-8') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            overall_dict = ast.literal_eval(row[1])

            if not overall_dict:
                print("Empty file: no co-occurrences found.")
                SeqIO.write([], output_fasta, "fasta")
                with open(output_tsv, 'w') as f:
                    f.write("read_name\tcount\tamplicon\tlocation\n")
            else:
                print("Step 1: Simulating mutated reads...")
                possible_reads = simulate_mutated_reads(
                    csv_file_path, fasta_file, amp_file, bed_file, sample_name)

                print("Step 2: Writing FASTA and counts TSV...")
                write_fasta_and_counts(possible_reads, sample_name, amplicon, location, output_fasta, output_tsv)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate unique amplicon reads and output raw counts')
    parser.add_argument("csv_file_path", help="Path to cooccurrences_summarized.csv")
    parser.add_argument("fasta_file",    help="Path to reference FASTA")
    parser.add_argument("amp_file",      help="Path to amp_dict YAML")
    parser.add_argument("bed_file",      help="Path to BED file with amplicon coordinates")
    parser.add_argument("sample_name",   help="Sample name for FASTA headers")
    parser.add_argument("amplicon",       help="Amplicon name (wildcard)")
    parser.add_argument("location",       help="Location name (wildcard)")
    parser.add_argument("output_fasta",  help="Output FASTA file path")
    parser.add_argument("output_tsv",    help="Output TSV file path (read name → count)")
    args = parser.parse_args()

    main(args.csv_file_path, args.fasta_file, args.amp_file, args.bed_file,
         args.sample_name, args.amplicon, args.location, args.output_fasta, args.output_tsv)
