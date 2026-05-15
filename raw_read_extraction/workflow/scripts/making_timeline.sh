#!/bin/bash

# Usage: ./making_timeline.sh <BASE_DIR> <OUTPUT_FILE>

BASE_DIR="$1"
OUTPUT_FILE="$2"

# Location mapping
declare -A LOCATIONS
LOCATIONS[05]="Lugano (TI)"
LOCATIONS[10]="Zürich (ZH)"
LOCATIONS[15]="Basel (BS)"
LOCATIONS[16]="Genève (GE)"
LOCATIONS[17]="Chur (GR)"
LOCATIONS[18]="Luzern (LU)"
LOCATIONS[19]="Altenrhein (SG)"
LOCATIONS[25]="Laupen (BE)"

# Create header
echo -e "sample\tbatch\treads\tproto\tlocation_code\tdate\tlocation" > "$OUTPUT_FILE"

# Process each directory
for dir in "$BASE_DIR"/*; do
    if [[ -d "$dir" ]]; then
        sample_name=$(basename "$dir")
        
        # Skip directories ending with _3x or similar patterns
        if [[ "$sample_name" =~ _[0-9]+x$ ]]; then
            continue
        fi
        
        # Extract components from sample name (e.g., E2_05_2025_04_03)
        if [[ "$sample_name" =~ ^([A-Z][0-9]+)_([0-9]+)_([0-9]{4})_([0-9]{2})_([0-9]{2})$ ]]; then
            location_code="${BASH_REMATCH[2]}"
            year="${BASH_REMATCH[3]}"
            month="${BASH_REMATCH[4]}"
            day="${BASH_REMATCH[5]}"
            date="$year-$month-$day"
            
            location="${LOCATIONS[$location_code]:-Unknown}"
            
            batch_dir=$(find "$dir" -maxdepth 1 -type d -name "2*" | head -1)
            if [[ -n "$batch_dir" ]]; then
                batch=$(basename "$batch_dir")
            else
                batch="unknown"
            fi
            
            reads="250"
            proto="dummy"
            
            echo -e "$sample_name\t$batch\t$reads\t$proto\t$location_code\t$date\t$location" >> "$OUTPUT_FILE"
        fi
    fi
done

echo "Timeline file generated: $OUTPUT_FILE"