import pysam
import argparse

parser = argparse.ArgumentParser(description="Filter BAM by start/end positions.")


parser.add_argument("-i","--input_bam", 
                    required=True,type=str,dest="i",
                    help="Input BAM filename")
parser.add_argument("-o","--output_bam", 
                    required=True,type=str,dest="o",
                    help="Output BAM filename")
parser.add_argument("-s","--target_start", 
                    required=True,type=int,dest="s",
                    help="Target start position (1-based)")
parser.add_argument("-e","--target_end", 
                    required=True,type=int,dest="e",
                    help="Target end position (1-based)")
# parser.add_argument("-b","--buffer", 
#                     required=False,type=int, dest="b",default=5, 
#                     help="Buffer size around target positions, defines a range of possible start and end positions")

args = parser.parse_args()


# def filter_bam_by_position(input_bam, output_bam, target_start, target_end, buffer=5):
    
#     start_range = range(target_start - buffer, target_start + buffer)
#     end_range = range(target_end - buffer, target_end + buffer)

#     with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
#          pysam.AlignmentFile(output_bam, "wb", template=in_bam) as out_bam:

#         for read in in_bam:
#             if read.is_unmapped:
#                 continue
            
#             start = read.reference_start
#             end = read.reference_end

#             spans_region = start <= target_start and end >= target_end

#             # TODO: don't need buffer anymore I guess
#             if (start in start_range) and (end in end_range) and spans_region:
#                 out_bam.write(read)

def filter_bam_by_position(input_bam, output_bam, target_start, target_end):
    
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", template=in_bam) as out_bam:

        for read in in_bam:
            if read.is_unmapped:
                continue
            
            start = read.reference_start
            end = read.reference_end

            spans_region = start <= target_start and end >= target_end

            if spans_region:
                out_bam.write(read)
            

def main():
    # filter_bam_by_position(args.i, args.o, args.s, args.e, args.b)
    filter_bam_by_position(args.i, args.o, args.s, args.e)

if __name__ == "__main__":
    main()

