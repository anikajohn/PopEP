import pysam
import argparse

parser = argparse.ArgumentParser(description="Convert Bam to Fasta, extracting only aligned sequences (excluding soft-clipped bases)")


parser.add_argument("-i","--input_bam", 
                    required=True,type=str,dest="i",
                    help="Input BAM filename")
parser.add_argument("-o","--output_fasta", 
                    required=True,type=str,dest="o",
                    help="Output FASTA filename")
parser.add_argument("-s","--sample_name",
                    required=False,type=str,dest="s",
                    help="Sample name to include in FASTA headers (optional)")
parser.add_argument("-a","--amplicon_name",
                    required=False,type=str,dest="a",
                    help="Amplicon name to include in FASTA headers (optional)")



args = parser.parse_args()


def extract_aligned_sequences(bam_file_path, output_fasta_path, sample_name=None, amplicon_name=None):
    """
    Extracts aligned sequences (excluding soft-clipped bases) from a BAM file and writes them to a FASTA file.

    Args:
        bam_file_path: Path to the input BAM file.
        output_fasta_path: Path to the output FASTA file.
        sample_name: Optional sample name to include in FASTA headers.
        amplicon_name: Optional amplicon name to include in FASTA headers.
    """
    with pysam.AlignmentFile(bam_file_path, "rb") as bam, \
        open(output_fasta_path, "w") as fasta_file:
        for read in bam.fetch():
                
            if read.is_unmapped:
                continue
            
            aligned_sequence = ""
            qpos = 0
            
            for operation, length in read.cigartuples:
                if operation == 0:  # Alignment match
                    aligned_sequence += read.query_sequence[qpos:qpos + length]
                    qpos += length
                elif operation == 7:  # Sequence match
                    aligned_sequence += read.query_sequence[qpos:qpos + length]
                    qpos += length
                elif operation == 8: # Sequence mismatch
                    aligned_sequence += read.query_sequence[qpos:qpos + length]
                    qpos += length
                elif operation == 2: # Deletion from reference (no read bases)
                    continue
                elif operation == 1:  # Insertion to reference skip, but advance qpos
                    qpos += length
                elif operation == 4:  # Soft clip
                    continue
                elif operation == 5: # Hard clip
                    continue
                
            if aligned_sequence:
                header = f">{read.query_name}"
                if sample_name:
                    header += f"_{sample_name}"
                if amplicon_name:
                    header += f"_{amplicon_name}"
                fasta_file.write(f"{header}\n{aligned_sequence}\n")

def main():
    extract_aligned_sequences(args.i, args.o, args.s, args.a)

if __name__ == "__main__":
    main()



