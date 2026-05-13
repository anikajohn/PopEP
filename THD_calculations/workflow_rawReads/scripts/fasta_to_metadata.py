#!/usr/bin/env python3
"""
Script to extract headers from a FASTA file and generate a metadata file.
Output: tab-separated file with strain and date columns.
"""

import sys
import re
import argparse

def extract_date_from_header(header):
    """
    Extract date from header format like: AV233803:AV083:2414523030:1:10801:0024:2668_A1_05_2025_01_29/20250221_2414523030_amplicon4
    Returns date in YYYY-MM-DD format.
    """
    # Pattern to match YYYY_MM_DD in the header
    match = re.search(r'(\d{4})_(\d{2})_(\d{2})', header)
    if match:
        year, month, day = match.groups()
        return f"{year}-{month}-{day}"
    return None

def parse_fasta_headers(fasta_file):
    """
    Parse FASTA file and extract headers (lines starting with >).
    """
    headers = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Remove the '>' and any trailing whitespace
                header = line[1:].strip()
                headers.append(header)
    return headers

def generate_metadata(fasta_file, output_file):
    """
    Generate metadata file from FASTA headers.
    """
    headers = parse_fasta_headers(fasta_file)
    
    with open(output_file, 'w') as f:
        # Write header row
        f.write("strain\tdate\n")
        
        for header in headers:
            date = extract_date_from_header(header)
            if date:
                f.write(f"{header}\t{date}\n")
            else:
                print(f"Warning: Could not extract date from header: {header}", file=sys.stderr)
                f.write(f"{header}\t\n")
    
    print(f"Metadata file generated: {output_file}")
    print(f"Total entries: {len(headers)}")

def main():
    parser = argparse.ArgumentParser(description='Generate metadata file from FASTA headers')
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('-o', '--output', help='Output metadata file (default: metadata.tsv)', 
                        default='metadata.tsv')
    
    args = parser.parse_args()
    generate_metadata(args.fasta_file, args.output)

if __name__ == '__main__':
    main()
